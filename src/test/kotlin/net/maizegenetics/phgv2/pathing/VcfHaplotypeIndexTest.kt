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
}
