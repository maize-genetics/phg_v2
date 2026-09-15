package net.maizegenetics.phgv2.pathing.ropebwt

import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.io.TempDir
import java.io.File
import java.nio.ByteBuffer
import java.nio.ByteOrder
import java.util.zip.GZIPOutputStream

class BuildPresenceTableTest {

    /** One lift point: a carrier position anchored to a reference position. */
    private data class Point(val carrierPos: Long, val refPos: Long, val carrierSeq: Int, val refSeq: Int)

    /**
     * Writes a `ropebwt3 lift` file. [pointsBySeq] is indexed by carrier sequence id, matching the
     * offset table the real format carries, and the points are written in carrier order.
     */
    private fun writeLift(file: File, pointsBySeq: List<List<Point>>) {
        val total = pointsBySeq.sumOf { it.size }
        val offsets = LongArray(pointsBySeq.size + 1)
        for (i in pointsBySeq.indices) offsets[i + 1] = offsets[i] + pointsBySeq[i].size
        val buffer = ByteBuffer.allocate(21 + 8 * offsets.size + 24 * total).order(ByteOrder.LITTLE_ENDIAN)
        buffer.put("LIFT".toByteArray(Charsets.ISO_8859_1))
        buffer.put(1)
        buffer.putLong(pointsBySeq.size.toLong())
        buffer.putLong(total.toLong())
        offsets.forEach { buffer.putLong(it) }
        for (points in pointsBySeq) {
            for (p in points) {
                buffer.putLong(p.carrierPos); buffer.putLong(p.refPos)
                buffer.putInt(p.carrierSeq); buffer.putInt(p.refSeq)
            }
        }
        file.writeBytes(buffer.array())
    }

    private fun writeLengths(file: File, entries: List<Pair<String, Long>>) {
        GZIPOutputStream(file.outputStream()).bufferedWriter().use { writer ->
            entries.forEach { (name, length) -> writer.write("$name\t$length\n") }
        }
    }

    /**
     * Three taxa over one reference chromosome of 300 bp with 100 bp windows, so three windows.
     *
     * Reference anchor grid, per window:
     *   window 0 (0-99):    10, 20, 30   -- founderA has all three, founderB has one
     *   window 1 (100-199): 150          -- founderA only
     *   window 2 (200-299): (none)       -- no anchors at all
     */
    private fun simpleCase(dir: File): Triple<File, File, BuildPresenceTable.Blocks> {
        val lift = File(dir, "test.lift")
        val lengths = File(dir, "test.len.gz")
        writeLengths(lengths, listOf(
            "B73_chr1" to 300L,        // seq 0: the reference chromosome
            "founderA_ctg1" to 250L,   // seq 1
            "founderB_ctg1" to 250L    // seq 2
        ))
        writeLift(lift, listOf(
            emptyList(),
            listOf(Point(0, 10, 1, 0), Point(5, 20, 1, 0), Point(9, 30, 1, 0), Point(40, 150, 1, 0)),
            listOf(Point(0, 20, 2, 0))
        ))
        return Triple(lift, lengths, BuildPresenceTable().buildBlocks(lift, lengths, "B73", 100))
    }

    @Test
    fun presenceIsTheFractionOfTheWindowAnchorGridATaxonHolds(@TempDir dir: File) {
        val (_, _, blocks) = simpleCase(dir)
        assertEquals(listOf("B73", "founderA", "founderB"), blocks.taxa)
        assertEquals(1, blocks.chroms.size)

        val chrom = blocks.chroms.single()
        assertEquals("chr1", chrom.name, "the reference prefix is stripped from the contig name")
        assertEquals(4, chrom.nWindows, "300 / 100 + 1")

        val n = blocks.taxa.size
        fun presence(window: Int, taxon: String) = chrom.presence!![window * n + blocks.taxa.indexOf(taxon)]

        // Window 0's grid is {10, 20, 30}: founderA holds all of it, founderB holds only 20.
        assertEquals(1.0f, presence(0, "founderA"))
        assertEquals(1.0f / 3.0f, presence(0, "founderB"))
        // Window 1's grid is {150}, which only founderA anchors to.
        assertEquals(1.0f, presence(1, "founderA"))
        assertEquals(0.0f, presence(1, "founderB"))
        // Window 2 has no anchors, so nobody is present rather than everybody.
        assertEquals(0.0f, presence(2, "founderA"))
        assertEquals(0.0f, presence(2, "founderB"))
    }

    @Test
    fun theReferenceHoldsEveryAnchorItHasNoPointsOfItsOwn(@TempDir dir: File) {
        val (_, _, blocks) = simpleCase(dir)
        val chrom = blocks.chroms.single()
        val n = blocks.taxa.size
        val slot = blocks.taxa.indexOf("B73")
        assertEquals(1.0f, chrom.presence!![0 * n + slot], "window with anchors")
        assertEquals(1.0f, chrom.presence!![1 * n + slot], "window with anchors")
        assertEquals(0.0f, chrom.presence!![2 * n + slot], "window with no anchors at all")
    }

    @Test
    fun duplicatePositionsFromSeveralContigsOfOneFounderCountOnce(@TempDir dir: File) {
        val lift = File(dir, "dup.lift")
        val lengths = File(dir, "dup.len.gz")
        writeLengths(lengths, listOf(
            "B73_chr1" to 100L,
            "founderA_ctg1" to 50L,
            "founderA_ctg2" to 50L
        ))
        // Both of founderA's contigs anchor to reference position 10; that is one anchor, not two,
        // so presence must stay at 1.0 rather than exceeding it.
        writeLift(lift, listOf(
            emptyList(),
            listOf(Point(0, 10, 1, 0)),
            listOf(Point(0, 10, 2, 0))
        ))
        val blocks = BuildPresenceTable().buildBlocks(lift, lengths, "B73", 100)
        val n = blocks.taxa.size
        val chrom = blocks.chroms.single()
        assertEquals(1.0f, chrom.presence!![0 * n + blocks.taxa.indexOf("founderA")])
    }

    @Test
    fun aReferenceChromosomeWithNoAnchorsGetsNoBlock(@TempDir dir: File) {
        val lift = File(dir, "empty.lift")
        val lengths = File(dir, "empty.len.gz")
        writeLengths(lengths, listOf(
            "B73_chr1" to 100L,
            "B73_chr2" to 100L,
            "founderA_ctg1" to 50L
        ))
        writeLift(lift, listOf(emptyList(), emptyList(), listOf(Point(0, 10, 2, 0))))
        val blocks = BuildPresenceTable().buildBlocks(lift, lengths, "B73", 100)
        assertEquals(listOf("chr1", "chr2"), blocks.chroms.map { it.name })
        assertEquals(2, blocks.chroms[0].nWindows)
        assertEquals(0, blocks.chroms[1].nWindows, "chr2 has no anchors")
        assertNull(blocks.chroms[1].presence)
    }

    @Test
    fun chromosomesAreOrderedNumericallyNotLexically(@TempDir dir: File) {
        val lift = File(dir, "order.lift")
        val lengths = File(dir, "order.len.gz")
        val names = listOf("B73_chr1", "B73_chr2", "B73_chr10", "founderA_ctg1")
        writeLengths(lengths, names.map { it to 100L })
        writeLift(lift, listOf(emptyList(), emptyList(), emptyList(), listOf(Point(0, 10, 3, 0))))
        val blocks = BuildPresenceTable().buildBlocks(lift, lengths, "B73", 100)
        assertEquals(listOf("chr1", "chr2", "chr10"), blocks.chroms.map { it.name },
            "chr10 must sort after chr2, not between chr1 and chr2")
    }

    @Test
    fun theWrittenFileIsReadBackByPresenceTable(@TempDir dir: File) {
        val (_, _, blocks) = simpleCase(dir)
        val out = File(dir, "presence.bin")
        BuildPresenceTable().writeTable(out, 100, blocks.taxa, blocks.chroms)

        val table = PresenceTable(out)
        assertEquals(100, table.windowSize)
        assertEquals(3, table.nTaxa)
        assertTrue(table.hasContig("chr1"))
        assertEquals(4, table.windowCount("chr1"))
        assertTrue(table.missingTaxa(listOf("B73", "founderA", "founderB")).isEmpty())
        assertEquals(listOf("nope"), table.missingTaxa(listOf("B73", "nope")))

        assertEquals(1.0, table.presence("chr1", 0, "founderA"), 1e-6)
        assertEquals(1.0 / 3.0, table.presence("chr1", 0, "founderB"), 1e-6)
        assertEquals(0.0, table.presence("chr1", 2, "founderB"), 1e-6)
        assertEquals(1.0, table.presence("chr1", 1, "B73"), 1e-6)
        assertTrue(table.presence("chr1", 0, "unknown").isNaN())
        assertTrue(table.presence("chrX", 0, "founderA").isNaN())
    }

    @Test
    fun aLiftFileWithTheWrongMagicIsRejected(@TempDir dir: File) {
        val lift = File(dir, "bad.lift")
        lift.writeBytes(ByteArray(64) { 0 })
        val lengths = File(dir, "bad.len.gz")
        writeLengths(lengths, listOf("B73_chr1" to 100L))
        val error = assertThrows(IllegalArgumentException::class.java) {
            BuildPresenceTable().buildBlocks(lift, lengths, "B73", 100)
        }
        assertTrue(error.message!!.contains("not a ropebwt3 lift file"))
    }

    @Test
    fun aLengthFileThatDisagreesWithTheLiftIsRejected(@TempDir dir: File) {
        val lift = File(dir, "mismatch.lift")
        val lengths = File(dir, "mismatch.len.gz")
        writeLift(lift, listOf(emptyList(), emptyList()))
        writeLengths(lengths, listOf("B73_chr1" to 100L))   // one name, but the lift holds two
        val error = assertThrows(IllegalArgumentException::class.java) {
            BuildPresenceTable().buildBlocks(lift, lengths, "B73", 100)
        }
        assertTrue(error.message!!.contains("lift has 2 sequences"))
    }

    @Test
    fun pointsAreAttributedToTheCarrierSequenceTheOffsetTableNames(@TempDir dir: File) {
        // Every point carries a refPos in the same window, so the only thing that can distinguish
        // the two founders is whether the offset walk assigns each point to the right carrier.
        val lift = File(dir, "carrier.lift")
        val lengths = File(dir, "carrier.len.gz")
        writeLengths(lengths, listOf(
            "B73_chr1" to 100L,
            "founderA_ctg1" to 50L,
            "founderB_ctg1" to 50L
        ))
        writeLift(lift, listOf(
            emptyList(),
            listOf(Point(0, 10, 1, 0), Point(1, 20, 1, 0)),
            listOf(Point(0, 30, 2, 0))
        ))
        val blocks = BuildPresenceTable().buildBlocks(lift, lengths, "B73", 100)
        val n = blocks.taxa.size
        val presence = blocks.chroms.single().presence!!
        // Grid is {10, 20, 30}: A holds two thirds of it, B one third.
        assertEquals(2.0f / 3.0f, presence[blocks.taxa.indexOf("founderA")], 1e-6f)
        assertEquals(1.0f / 3.0f, presence[blocks.taxa.indexOf("founderB")], 1e-6f)
    }
}
