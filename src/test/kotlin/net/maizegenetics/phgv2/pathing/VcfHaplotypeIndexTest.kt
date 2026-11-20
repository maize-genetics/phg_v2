package net.maizegenetics.phgv2.pathing

import htsjdk.variant.variantcontext.Allele
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.Position
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNull
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
class VcfHaplotypeIndexTest {

    companion object {
        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.testOutputDir).deleteRecursively()
            File(TestExtension.testOutputDir).mkdirs()
        }

        @JvmStatic
        @AfterAll
        fun tearDown() {
            //comment out the following line to inspect the test results after the tests have been run
            File(TestExtension.testOutputDir).deleteRecursively()
        }
    }

    @Test
    fun testAlleleKeyPlainAllele() {
        assertEquals("A", VcfHaplotypeIndexUtils.alleleKey(Allele.create("A", true)))
        assertEquals("T", VcfHaplotypeIndexUtils.alleleKey(Allele.create("T", false)))
    }

    @Test
    fun testAlleleKeySymbolicAllele() {
        val delAllele = Allele.create("<DEL>", false)

        //Pin the htsjdk trap: baseString silently returns "" for symbolic alleles because
        //isSymbolic() alleles have no actual bases. If a future refactor swaps alleleKey to use
        //baseString instead of displayString, this assertion documents why that is wrong.
        assertEquals("", delAllele.baseString)

        assertEquals("<DEL>", VcfHaplotypeIndexUtils.alleleKey(delAllele))
    }

    @Test
    fun testAlleleKeyRejectsNoCall() {
        assertFailsWith<IllegalArgumentException> { VcfHaplotypeIndexUtils.alleleKey(Allele.NO_CALL) }
    }

    @Test
    fun testFormatAndParseEntryRoundTrip() {
        val entry = VcfHaplotypeIndexEntry(
            contig = "1",
            position = 1001,
            allele = "A",
            refRange = ReferenceRange("1", 1001, 5500),
            hapIds = listOf("hap1", "hap2", "hap3")
        )

        val formatted = VcfHaplotypeIndexUtils.formatEntry(entry)
        assertEquals(4, formatted.count { it == '\t' })
        assertEquals("1\t1001\tA\t1:1001-5500\thap1,hap2,hap3", formatted)

        val parsed = VcfHaplotypeIndexUtils.parseEntry(formatted)
        assertEquals(entry, parsed)
    }

    @Test
    fun testParseEntryRejectsMalformedRow() {
        assertFailsWith<IllegalArgumentException> { VcfHaplotypeIndexUtils.parseEntry("1\t1001\tA\t1:1001-5500") }
    }

    @Test
    fun testParseRefRangeWithDashInContig() {
        val range = VcfHaplotypeIndexUtils.parseRefRange("scaffold-1:100-200")
        assertEquals(ReferenceRange("scaffold-1", 100, 200), range)

        //ReferenceRange.parse splits on the FIRST "-", so on this same string it tries to parse
        //"1:100-200" as an Int and throws - this is why parseRefRange cannot just delegate to
        //ReferenceRange.parse.
        assertFailsWith<NumberFormatException> { ReferenceRange.parse("scaffold-1:100-200") }
    }

    @Test
    fun testParseRefRangeSimpleContig() {
        assertEquals(ReferenceRange("1", 1, 100), VcfHaplotypeIndexUtils.parseRefRange("1:1-100"))
    }

    @Test
    fun testWriteAndReadEntriesRoundTrip() {
        val entries = listOf(
            VcfHaplotypeIndexEntry("1", 10, "A", ReferenceRange("1", 1, 100), listOf("hapA1", "hapA2")),
            VcfHaplotypeIndexEntry("1", 10, "T", ReferenceRange("1", 1, 100), listOf("hapB1")),
            VcfHaplotypeIndexEntry("1", 150, "C", ReferenceRange("1", 101, 200), listOf("hapC1"))
        )

        val outputFile = "${TestExtension.testOutputDir}roundTrip.txt"
        VcfHaplotypeIndexUtils.writeIndex(outputFile, entries, commandLine = "phg build-vcf-haplotype-index ...",
            hvcfDir = "hvcfDir", panelVcf = "panel.vcf", graphChecksum = "checksum123")

        val readBack = VcfHaplotypeIndexUtils.readEntries(outputFile)
        assertEquals(entries, readBack)
    }

    @Test
    fun testReadIndexGroupsByPosition() {
        val entries = listOf(
            VcfHaplotypeIndexEntry("1", 10, "A", ReferenceRange("1", 1, 100), listOf("hapA1", "hapA2")),
            VcfHaplotypeIndexEntry("1", 10, "T", ReferenceRange("1", 1, 100), listOf("hapB1")),
            VcfHaplotypeIndexEntry("1", 150, "C", ReferenceRange("1", 101, 200), listOf("hapC1"))
        )
        val outputFile = "${TestExtension.testOutputDir}groupedIndex.txt"
        VcfHaplotypeIndexUtils.writeIndex(outputFile, entries)

        val index = VcfHaplotypeIndexUtils.readIndex(outputFile)
        assertEquals(2, index.size)

        val pos10 = index[Position("1", 10)]!!
        assertEquals(ReferenceRange("1", 1, 100), pos10.refRange)
        assertEquals(mapOf("A" to listOf("hapA1", "hapA2"), "T" to listOf("hapB1")), pos10.alleleToHapIds)

        val pos150 = index[Position("1", 150)]!!
        assertEquals(ReferenceRange("1", 101, 200), pos150.refRange)
        assertEquals(mapOf("C" to listOf("hapC1")), pos150.alleleToHapIds)
    }

    @Test
    fun testReadIndexMergesDuplicatePositions() {
        val entries = listOf(
            VcfHaplotypeIndexEntry("1", 10, "A", ReferenceRange("1", 1, 100), listOf("hap2", "hap1")),
            VcfHaplotypeIndexEntry("1", 10, "A", ReferenceRange("1", 1, 100), listOf("hap3", "hap1"))
        )
        val outputFile = "${TestExtension.testOutputDir}dupIndex.txt"
        VcfHaplotypeIndexUtils.writeIndex(outputFile, entries)

        val index = VcfHaplotypeIndexUtils.readIndex(outputFile)
        assertEquals(1, index.size)
        assertEquals(listOf("hap1", "hap2", "hap3"), index[Position("1", 10)]!!.alleleToHapIds["A"])
    }

    @Test
    fun testReadEntriesRejectsFileWithoutColumnHeader() {
        val outputFile = "${TestExtension.testOutputDir}noHeader.txt"
        File(outputFile).writeText("1\t10\tA\t1:1-100\thap1\n")
        assertFailsWith<IllegalArgumentException> { VcfHaplotypeIndexUtils.readEntries(outputFile) }
    }

    @Test
    fun testReadEntriesSkipsCommentLines() {
        val outputFile = "${TestExtension.testOutputDir}withComments.txt"
        File(outputFile).writeText(
            "${VcfHaplotypeIndexUtils.MAGIC}\n" +
                "#version=${VcfHaplotypeIndexUtils.VERSION}\n" +
                "#PHG Command: some command\n" +
                "${VcfHaplotypeIndexUtils.COLUMN_HEADER}\n" +
                "1\t10\tA\t1:1-100\thap1\n" +
                "#a stray trailing comment\n"
        )
        val entries = VcfHaplotypeIndexUtils.readEntries(outputFile)
        assertEquals(1, entries.size)
        assertEquals("A", entries[0].allele)
    }

    @Test
    fun testReadEntriesRejectsMalformedRow() {
        val outputFile = "${TestExtension.testOutputDir}malformed.txt"
        File(outputFile).writeText("${VcfHaplotypeIndexUtils.COLUMN_HEADER}\n1\t10\tA\t1:1-100\n")
        assertFailsWith<IllegalArgumentException> { VcfHaplotypeIndexUtils.readEntries(outputFile) }
    }

    @Test
    fun testReadEmptyIndex() {
        val outputFile = "${TestExtension.testOutputDir}emptyIndex.txt"
        VcfHaplotypeIndexUtils.writeIndex(outputFile, emptyList())
        assertTrue(VcfHaplotypeIndexUtils.readEntries(outputFile).isEmpty())
        assertTrue(VcfHaplotypeIndexUtils.readIndex(outputFile).isEmpty())
    }

    @Test
    fun testReadGraphChecksum() {
        val withChecksum = "${TestExtension.testOutputDir}withChecksum.txt"
        VcfHaplotypeIndexUtils.writeIndex(withChecksum, emptyList(), graphChecksum = "abc123")
        assertEquals("abc123", VcfHaplotypeIndexUtils.readGraphChecksum(withChecksum))

        val withoutChecksum = "${TestExtension.testOutputDir}withoutChecksum.txt"
        VcfHaplotypeIndexUtils.writeIndex(withoutChecksum, emptyList())
        assertNull(VcfHaplotypeIndexUtils.readGraphChecksum(withoutChecksum))
    }

    @Test
    fun testGzipRoundTrip() {
        val entries = listOf(
            VcfHaplotypeIndexEntry("1", 10, "A", ReferenceRange("1", 1, 100), listOf("hapA1", "hapA2"))
        )
        val outputFile = "${TestExtension.testOutputDir}gzippedIndex.txt.gz"
        VcfHaplotypeIndexUtils.writeIndex(outputFile, entries)

        val readBack = VcfHaplotypeIndexUtils.readEntries(outputFile)
        assertEquals(entries, readBack)
    }

    //region Streaming reader and contig-order scanning

    private fun writeIndexText(name: String, rows: String): String {
        val file = File("${TestExtension.testOutputDir}stream_$name")
        file.writeText(
            "${VcfHaplotypeIndexUtils.MAGIC}\n#version=1.0\n${VcfHaplotypeIndexUtils.COLUMN_HEADER}\n" +
                rows.trimIndent() + "\n"
        )
        return file.path
    }

    private fun readAllPositions(file: String): List<VcfHaplotypeIndexPosition> {
        val result = mutableListOf<VcfHaplotypeIndexPosition>()
        VcfHaplotypeIndexStreamReader(file).use { reader ->
            while (true) result.add(reader.next() ?: break)
        }
        return result
    }

    @Test
    fun testStreamReaderMatchesReadIndex() {
        val file = writeIndexText("matches.txt", """
            1	10	A	1:1-200	h1
            1	10	T	1:1-200	h1,h2
            1	10	T	1:1-200	h3
            1	20	C	1:1-200	h3
            1	250	G	1:201-400	h9,h8
            2	5	A	2:1-100	h4
        """)

        val streamed = readAllPositions(file)
        val loaded = VcfHaplotypeIndexUtils.readIndex(file)

        //same positions, same order, same merged/sorted allele -> hapIds lookups
        assertEquals(loaded.values.toList(), streamed)
        //duplicate rows for 1:10 T are unioned and sorted
        assertEquals(listOf("h1", "h2", "h3"), streamed[0].alleleToHapIds["T"])
        //an unsorted hapId list is sorted, so it compares equal to the same set elsewhere
        assertEquals(listOf("h8", "h9"), streamed[2].alleleToHapIds["G"])
    }

    @Test
    fun testStreamReaderPoolsHapIdStringsWithinRange() {
        val file = writeIndexText("pool.txt", """
            1	10	A	1:1-200	h1,h2
            1	20	A	1:1-200	h1,h2
        """)
        val (first, second) = readAllPositions(file)
        //equal hapIds from different rows of the same range share one String instance
        assertTrue(first.alleleToHapIds.getValue("A")[0] === second.alleleToHapIds.getValue("A")[0])
    }

    @Test
    fun testStreamReaderEmptyIndex() {
        val file = writeIndexText("empty.txt", "")
        VcfHaplotypeIndexStreamReader(file).use { assertNull(it.next()) }
    }

    @Test
    fun testStreamReaderGzip() {
        val rows = "1\t10\tA\t1:1-200\th1\n"
        val gz = "${TestExtension.testOutputDir}stream_index.txt.gz"
        java.util.zip.GZIPOutputStream(File(gz).outputStream()).bufferedWriter().use {
            it.write("${VcfHaplotypeIndexUtils.COLUMN_HEADER}\n$rows")
        }
        val positions = readAllPositions(gz)
        assertEquals(listOf(Position("1", 10)), positions.map { it.position })
    }

    @Test
    fun testStreamReaderRejectsUnsortedPositions() {
        val file = writeIndexText("unsortedPos.txt", """
            1	20	A	1:1-200	h1
            1	10	A	1:1-200	h2
        """)
        val error = assertFailsWith<IllegalArgumentException> { readAllPositions(file) }
        assertTrue(error.message!!.contains("not sorted by position"), error.message)
    }

    @Test
    fun testStreamReaderRejectsNonContiguousContig() {
        val file = writeIndexText("nonContiguous.txt", """
            1	10	A	1:1-200	h1
            2	10	A	2:1-200	h2
            1	30	A	1:1-200	h3
        """)
        val error = assertFailsWith<IllegalArgumentException> { readAllPositions(file) }
        assertTrue(error.message!!.contains("non-contiguous"), error.message)
    }

    @Test
    fun testStreamReaderRejectsMissingOrLateColumnHeader() {
        val noHeader = File("${TestExtension.testOutputDir}stream_noHeader.txt").apply {
            writeText("${VcfHaplotypeIndexUtils.MAGIC}\n#version=1.0\n")
        }.path
        assertFailsWith<IllegalArgumentException> { VcfHaplotypeIndexStreamReader(noHeader) }

        val dataFirst = File("${TestExtension.testOutputDir}stream_dataFirst.txt").apply {
            writeText("1\t10\tA\t1:1-200\th1\n${VcfHaplotypeIndexUtils.COLUMN_HEADER}\n")
        }.path
        assertFailsWith<IllegalArgumentException> { VcfHaplotypeIndexStreamReader(dataFirst) }
    }

    @Test
    fun testScanContigOrderHandlesCrlfLongNamesAndHeaderOnlyFiles() {
        val longContig = "scaffold_" + "x".repeat(1000) //longer than the scanner's initial 256-byte buffer
        val crlf = File("${TestExtension.testOutputDir}scan_crlf.txt").apply {
            writeText("#header\r\nchrA\t1\r\nchrA\t2\r\n$longContig\t1\r\nchrB\r\n")
        }.path
        //"chrB" has no tab at all (malformed line) - its whole line is still read as the contig
        assertEquals(listOf("chrA", longContig, "chrB"), VcfHaplotypeIndexUtils.scanContigOrder(crlf))

        val headerOnly = File("${TestExtension.testOutputDir}scan_headerOnly.txt").apply { writeText("#only\n#header\n") }.path
        assertEquals(emptyList(), VcfHaplotypeIndexUtils.scanContigOrder(headerOnly))

        //no trailing newline on the last line
        val noTrailingNewline = File("${TestExtension.testOutputDir}scan_noNewline.txt").apply { writeText("chrA\t1\nchrB") }.path
        assertEquals(listOf("chrA", "chrB"), VcfHaplotypeIndexUtils.scanContigOrder(noTrailingNewline))
    }

    @Test
    fun testScanContigOrderRejectsNonContiguousContig() {
        val file = File("${TestExtension.testOutputDir}scan_nonContiguous.txt").apply { writeText("a\t1\nb\t1\na\t2\n") }.path
        val error = assertFailsWith<IllegalArgumentException> { VcfHaplotypeIndexUtils.scanContigOrder(file) }
        assertTrue(error.message!!.contains("'a'"), error.message)
    }

    @Test
    fun testScanContigOrderBcf() {
        val header = htsjdk.variant.vcf.VCFHeader(emptySet(), listOf("S1")).apply {
            sequenceDictionary = htsjdk.samtools.SAMSequenceDictionary(
                listOf(htsjdk.samtools.SAMSequenceRecord("chr2", 1000), htsjdk.samtools.SAMSequenceRecord("chr10", 1000))
            )
            addMetaDataLine(htsjdk.variant.vcf.VCFFormatHeaderLine("GT", 1, htsjdk.variant.vcf.VCFHeaderLineType.String, "Genotype"))
        }
        val bcf = File("${TestExtension.testOutputDir}scan.bcf")
        htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder()
            .setOutputFile(bcf)
            .setOutputFileType(htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType.BCF)
            .unsetOption(htsjdk.variant.variantcontext.writer.Options.INDEX_ON_THE_FLY)
            .build().use { writer ->
                writer.writeHeader(header)
                listOf("chr2" to 5, "chr2" to 9, "chr10" to 3).forEach { (contig, pos) ->
                    val alleles = listOf(Allele.create("A", true), Allele.create("T", false))
                    writer.add(
                        htsjdk.variant.variantcontext.VariantContextBuilder("test", contig, pos.toLong(), pos.toLong(), alleles)
                            .genotypes(htsjdk.variant.variantcontext.GenotypeBuilder("S1", alleles).make())
                            .make()
                    )
                }
            }
        assertEquals(listOf("chr2", "chr10"), VcfHaplotypeIndexUtils.scanContigOrder(bcf.path))
    }

    //endregion
}
