package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.testing.test
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertFailsWith
import kotlin.test.assertFalse
import kotlin.test.assertNull
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
class BuildVcfHaplotypeIndexTest {

    companion object {
        const val hvcfDir = "data/test/vcfHaplotypeIndex/hvcfs"
        const val panelVcf = "data/test/vcfHaplotypeIndex/testPanel.vcf"
        const val expectedIndexFile = "data/test/vcfHaplotypeIndex/expectedIndex.txt"

        const val realHvcfDir = "data/test/hvcf2vcf/asmHvcfs/hvcf_files"
        const val realPanelVcf = "data/test/hvcf2vcf/asmPangenomeVcf/merged.vcf"

        val testGraph by lazy { HaplotypeGraph(hvcfDir) }

        val rangeA = ReferenceRange("1", 1, 100)
        val rangeB = ReferenceRange("1", 101, 200)

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

        //Keyed by (contig, start) - not just start - since the fixture panel VCF deliberately
        //has records on two different contigs sharing the same position (1:10 and 9:10).
        fun panelRecords(): Map<Pair<String, Int>, VariantContext> {
            return VCFFileReader(File(panelVcf), false).use { reader ->
                reader.associateBy { it.contig to it.start }
            }
        }
    }

    //region Fixture sanity checks (verifying the hand-written hvcf fixtures parse as intended)

    @Test
    fun testFixtureGraphShape() {
        assertEquals(listOf(rangeA, rangeB), testGraph.ranges())
        assertEquals(listOf("TestA", "TestB"), testGraph.samples())
    }

    @Test
    fun testFixtureGraphHaplotypesRangeA() {
        val hapIds = testGraph.sampleGameteToHaplotypeId(rangeA)
        assertEquals("hapA1g0", hapIds[SampleGamete("TestA", 0)])
        assertEquals("hapA1g1", hapIds[SampleGamete("TestA", 1)])
        assertEquals("hapB1g0", hapIds[SampleGamete("TestB", 0)])
        assertEquals("hapA1g0", hapIds[SampleGamete("TestB", 1)])
    }

    @Test
    fun testFixtureGraphHaplotypesRangeBOmitsTestB() {
        val hapIds = testGraph.sampleGameteToHaplotypeId(rangeB)
        assertEquals(setOf(SampleGamete("TestA", 0), SampleGamete("TestA", 1)), hapIds.keys)
        assertEquals("hapA2g0", hapIds[SampleGamete("TestA", 0)])
        assertEquals("hapA2g1", hapIds[SampleGamete("TestA", 1)])
    }

    //endregion

    //region CLI validation

    @Test
    fun testCliktParams() {
        val cmd = BuildVcfHaplotypeIndex()

        val noHvcfDir = cmd.test("--panel-vcf testDir.vcf --output-file testDir.txt")
        assertEquals(1, noHvcfDir.statusCode)
        assertEquals(
            "Usage: build-vcf-haplotype-index [<options>]\n\nError: missing option --hvcf-dir\n",
            noHvcfDir.stderr
        )

        val noPanelVcf = cmd.test("--hvcf-dir $hvcfDir --output-file testDir.txt")
        assertEquals(1, noPanelVcf.statusCode)
        assertEquals(
            "Usage: build-vcf-haplotype-index [<options>]\n\nError: missing option --panel-vcf\n",
            noPanelVcf.stderr
        )

        val noOutputFile = cmd.test("--hvcf-dir $hvcfDir --panel-vcf $panelVcf")
        assertEquals(1, noOutputFile.statusCode)
        assertEquals(
            "Usage: build-vcf-haplotype-index [<options>]\n\nError: missing option --output-file\n",
            noOutputFile.stderr
        )
    }

    @Test
    fun testCliktParamsInvalidPaths() {
        val cmd = BuildVcfHaplotypeIndex()

        val badHvcfDir = cmd.test("--hvcf-dir notADir --panel-vcf $panelVcf --output-file testDir.txt")
        assertEquals(1, badHvcfDir.statusCode)
        assertTrue(badHvcfDir.stderr.contains("is not a valid directory"))

        val badPanelVcf = cmd.test("--hvcf-dir $hvcfDir --panel-vcf notAFile.vcf --output-file testDir.txt")
        assertEquals(1, badPanelVcf.statusCode)
        assertTrue(badPanelVcf.stderr.contains("is not a valid file"))
    }

    //endregion

    //region RefRangeLocator

    @Test
    fun testRefRangeLocatorFindRange() {
        val ranges = mapOf(
            "1" to listOf(ReferenceRange("1", 1, 100), ReferenceRange("1", 101, 200)),
            "2" to listOf(ReferenceRange("2", 1, 50))
        )
        val locator = RefRangeLocator(ranges)

        assertEquals(ReferenceRange("1", 1, 100), locator.findRange("1", 1))
        assertEquals(ReferenceRange("1", 1, 100), locator.findRange("1", 100))
        assertEquals(ReferenceRange("1", 101, 200), locator.findRange("1", 101))
        assertEquals(ReferenceRange("1", 101, 200), locator.findRange("1", 200))
        assertNull(locator.findRange("1", 201))
        assertNull(locator.findRange("1", 0))
        assertEquals(ReferenceRange("2", 1, 50), locator.findRange("2", 25))
        assertNull(locator.findRange("2", 60))
        assertNull(locator.findRange("3", 10))
    }

    @Test
    fun testRefRangeLocatorSingleAndEmptyContigs() {
        val ranges = mapOf(
            "1" to listOf(ReferenceRange("1", 1, 100)),
            "2" to emptyList<ReferenceRange>()
        )
        val locator = RefRangeLocator(ranges)
        assertEquals(ReferenceRange("1", 1, 100), locator.findRange("1", 50))
        assertNull(locator.findRange("2", 1))

        val emptyLocator = RefRangeLocator(emptyMap())
        assertNull(emptyLocator.findRange("1", 1))
    }

    @Test
    fun testRefRangeLocatorMixedContigNames() {
        val ranges = mapOf(
            "1" to listOf(ReferenceRange("1", 1, 100)),
            "2" to listOf(ReferenceRange("2", 1, 100)),
            "10" to listOf(ReferenceRange("10", 1, 100)),
            "chrUn" to listOf(ReferenceRange("chrUn", 1, 100))
        )
        val locator = RefRangeLocator(ranges)

        assertEquals(ReferenceRange("1", 1, 100), locator.findRange("1", 50))
        assertEquals(ReferenceRange("2", 1, 100), locator.findRange("2", 50))
        assertEquals(ReferenceRange("10", 1, 100), locator.findRange("10", 50))
        assertEquals(ReferenceRange("chrUn", 1, 100), locator.findRange("chrUn", 50))
    }

    @Test
    fun testRefRangeLocatorFromGraphAgreesWithGraphRanges() {
        val locator = RefRangeLocator(testGraph)
        for (range in testGraph.ranges()) {
            assertEquals(range, locator.findRange(range.contig, range.start))
            assertEquals(range, locator.findRange(range.contig, range.end))
            assertEquals(range, locator.findRange(range.contig, (range.start + range.end) / 2))
        }
    }

    @Test
    fun testRefRangeLocatorDetectOverlapsNoneOnFixtureGraph() {
        assertTrue(RefRangeLocator(testGraph).detectOverlaps().isEmpty())
    }

    @Test
    fun testRefRangeLocatorDetectOverlapsFound() {
        val ranges = mapOf("1" to listOf(ReferenceRange("1", 1, 100), ReferenceRange("1", 50, 200)))
        val locator = RefRangeLocator(ranges)
        val overlaps = locator.detectOverlaps()
        assertEquals(1, overlaps.size)

        //a position within the overlap resolves to the range with the greater start
        assertEquals(ReferenceRange("1", 50, 200), locator.findRange("1", 75))
    }

    //endregion

    //region HapIdsByRangeCache

    @Test
    fun testHapIdsByRangeCache() {
        val cache = HapIdsByRangeCache(testGraph)

        val first = cache.hapIds(rangeA)
        val second = cache.hapIds(rangeA)
        val other = cache.hapIds(rangeB)
        val backToFirst = cache.hapIds(rangeA)

        assertEquals(testGraph.sampleGameteToHaplotypeId(rangeA), first)
        assertEquals(first, second)
        assertEquals(testGraph.sampleGameteToHaplotypeId(rangeB), other)
        assertEquals(first, backToFirst)
    }

    //endregion

    //region processRecord

    @Test
    fun testProcessRecordBiallelic() {
        val record = panelRecords().getValue("1" to 10)
        val locator = RefRangeLocator(testGraph)
        val cmd = BuildVcfHaplotypeIndex()
        val stats = VcfHaplotypeIndexStats()

        val entries = cmd.processRecord(record, locator, testGraph.samples().toSet(), testGraph::sampleGameteToHaplotypeId, stats)
            .associateBy { it.allele }

        assertEquals(2, entries.size)
        assertEquals(listOf("hapA1g0", "hapA1g1"), entries.getValue("A").hapIds)
        assertEquals(listOf("hapA1g0", "hapB1g0"), entries.getValue("T").hapIds)
        assertEquals(1, stats.positionsIndexed)
    }

    @Test
    fun testProcessRecordMultiAllelic() {
        val record = panelRecords().getValue("1" to 50)
        val locator = RefRangeLocator(testGraph)
        val cmd = BuildVcfHaplotypeIndex()
        val stats = VcfHaplotypeIndexStats()

        val entries = cmd.processRecord(record, locator, testGraph.samples().toSet(), testGraph::sampleGameteToHaplotypeId, stats)

        assertEquals(listOf("C", "G", "T"), entries.map { it.allele })
        val byAllele = entries.associateBy { it.allele }
        assertEquals(listOf("hapA1g0"), byAllele.getValue("C").hapIds)
        assertEquals(listOf("hapA1g1", "hapB1g0"), byAllele.getValue("G").hapIds)
        assertEquals(listOf("hapA1g0"), byAllele.getValue("T").hapIds)
    }

    @Test
    fun testProcessRecordNoCall() {
        val record = panelRecords().getValue("1" to 60)
        val locator = RefRangeLocator(testGraph)
        val cmd = BuildVcfHaplotypeIndex()
        val stats = VcfHaplotypeIndexStats()

        val entries = cmd.processRecord(record, locator, testGraph.samples().toSet(), testGraph::sampleGameteToHaplotypeId, stats)
            .associateBy { it.allele }

        assertEquals(2, stats.noCallAlleles)
        assertEquals(listOf("hapB1g0"), entries.getValue("A").hapIds)
        assertEquals(listOf("hapA1g0"), entries.getValue("T").hapIds)
    }

    @Test
    fun testProcessRecordSymbolicAllele() {
        val record = panelRecords().getValue("1" to 100)
        val locator = RefRangeLocator(testGraph)
        val cmd = BuildVcfHaplotypeIndex()
        val stats = VcfHaplotypeIndexStats()

        val entries = cmd.processRecord(record, locator, testGraph.samples().toSet(), testGraph::sampleGameteToHaplotypeId, stats)
            .associateBy { it.allele }

        assertEquals(2, entries.size)
        assertTrue(entries.containsKey("<DEL>"))
        assertEquals(listOf("hapA1g1"), entries.getValue("<DEL>").hapIds)
        assertEquals(listOf("hapA1g0", "hapB1g0"), entries.getValue("G").hapIds)
    }

    @Test
    fun testProcessRecordSampleMissingHaplotypeInRange() {
        val record = panelRecords().getValue("1" to 101)
        val locator = RefRangeLocator(testGraph)
        val cmd = BuildVcfHaplotypeIndex()
        val stats = VcfHaplotypeIndexStats()

        val entries = cmd.processRecord(record, locator, testGraph.samples().toSet(), testGraph::sampleGameteToHaplotypeId, stats)
            .associateBy { it.allele }

        assertEquals(2, stats.gametesWithNoHaplotypeInRange)
        assertEquals(listOf("hapA2g0"), entries.getValue("A").hapIds)
        assertEquals(listOf("hapA2g1"), entries.getValue("T").hapIds)
    }

    @Test
    fun testProcessRecordAlleleWithNoHapIdsIsSkipped() {
        val record = panelRecords().getValue("1" to 150)
        val locator = RefRangeLocator(testGraph)
        val cmd = BuildVcfHaplotypeIndex()
        val stats = VcfHaplotypeIndexStats()

        val entries = cmd.processRecord(record, locator, testGraph.samples().toSet(), testGraph::sampleGameteToHaplotypeId, stats)

        assertEquals(1, entries.size)
        assertEquals("A", entries[0].allele)
        assertEquals(1, stats.allelesWithNoHapIds)
    }

    @Test
    fun testProcessRecordUncoveredPositionAndUnknownContig() {
        val records = panelRecords()
        val locator = RefRangeLocator(testGraph)
        val cmd = BuildVcfHaplotypeIndex()
        val stats = VcfHaplotypeIndexStats()

        val uncoveredEntries = cmd.processRecord(records.getValue("1" to 201), locator, testGraph.samples().toSet(), testGraph::sampleGameteToHaplotypeId, stats)
        assertTrue(uncoveredEntries.isEmpty())

        val unknownContigEntries = cmd.processRecord(records.getValue("9" to 10), locator, testGraph.samples().toSet(), testGraph::sampleGameteToHaplotypeId, stats)
        assertTrue(unknownContigEntries.isEmpty())

        assertEquals(2, stats.positionsWithNoRefRange)
        assertEquals(0, stats.positionsIndexed)
    }

    @Test
    fun testProcessRecordSampleNotInGraphDoesNotThrow() {
        val locator = RefRangeLocator(testGraph)
        val cmd = BuildVcfHaplotypeIndex()
        val stats = VcfHaplotypeIndexStats()
        val graphSampleNames = testGraph.samples().toSet()

        panelRecords().values.forEach { record ->
            cmd.processRecord(record, locator, graphSampleNames, testGraph::sampleGameteToHaplotypeId, stats)
        }

        assertEquals(setOf("Ghost"), stats.samplesNotInGraph)
        //Ghost is diploid and appears in every record, but only the 6 records whose position
        //falls inside a reference range (1:10, 1:50, 1:60, 1:100, 1:101, 1:150) ever reach the
        //genotype loop - the 2 uncovered records (1:201, 9:10) return before Ghost is examined.
        assertEquals(12, stats.gametesFromSamplesNotInGraph)

        //Documents WHY processRecord routes through sampleGameteToHaplotypeId(range) instead of
        //graph.sampleToHapId(range, sample): the latter throws for an unknown sample name, which
        //would abort the whole run the first time a panel VCF contained an extra sample.
        assertFailsWith<IllegalArgumentException> {
            testGraph.sampleToHapId(rangeA, SampleGamete("Ghost", 0))
        }
    }

    //endregion

    //region End-to-end

    @Test
    fun testEndToEndSmallFixture() {
        val outputFile = "${TestExtension.testOutputDir}smallIndex.txt"
        val result = BuildVcfHaplotypeIndex().test("--hvcf-dir $hvcfDir --panel-vcf $panelVcf --output-file $outputFile")
        assertEquals(0, result.statusCode)

        val dataLines = File(outputFile).readLines().filterNot { it.startsWith("#") }
        val expectedLines = File(expectedIndexFile).readLines().filter { it.isNotBlank() }
        assertEquals(expectedLines, dataLines)

        assertTrue(File(outputFile).readLines().contains(VcfHaplotypeIndexUtils.COLUMN_HEADER))
    }

    @Test
    fun testEndToEndGzipOutput() {
        val outputFile = "${TestExtension.testOutputDir}smallIndex.txt.gz"
        val result = BuildVcfHaplotypeIndex().test("--hvcf-dir $hvcfDir --panel-vcf $panelVcf --output-file $outputFile")
        assertEquals(0, result.statusCode)

        val entries = VcfHaplotypeIndexUtils.readEntries(outputFile)
        val expectedLines = File(expectedIndexFile).readLines().filter { it.isNotBlank() }
        assertEquals(expectedLines.size, entries.size)
    }

    @Test
    fun testEndToEndReaderConsumesWriterOutput() {
        val outputFile = "${TestExtension.testOutputDir}handoffIndex.txt"
        val result = BuildVcfHaplotypeIndex().test("--hvcf-dir $hvcfDir --panel-vcf $panelVcf --output-file $outputFile")
        assertEquals(0, result.statusCode)

        val index = VcfHaplotypeIndexUtils.readIndex(outputFile)
        val pos50 = index[net.maizegenetics.phgv2.utils.Position("1", 50)]!!
        assertEquals(setOf("C", "G", "T"), pos50.alleleToHapIds.keys)
        assertEquals(listOf("hapA1g0"), pos50.alleleToHapIds["C"])
        assertEquals(listOf("hapA1g1", "hapB1g0"), pos50.alleleToHapIds["G"])
        assertEquals(listOf("hapA1g0"), pos50.alleleToHapIds["T"])
    }

    @Test
    fun testEndToEndRealFixtures() {
        val outputFile = "${TestExtension.testOutputDir}realIndex.txt"
        val result = BuildVcfHaplotypeIndex().test("--hvcf-dir $realHvcfDir --panel-vcf $realPanelVcf --output-file $outputFile")
        assertEquals(0, result.statusCode)

        val entries = VcfHaplotypeIndexUtils.readEntries(outputFile)

        val distinctPositions = entries.map { it.positionKey }.toSet()
        assertEquals(1737, distinctPositions.size)

        entries.forEach { entry ->
            assertTrue(entry.contig == "1" || entry.contig == "2")
            assertTrue(entry.position <= 16500)
            assertEquals(entry.contig, entry.refRange.contig)
            assertTrue(entry.position in entry.refRange.start..entry.refRange.end)
            assertTrue(entry.hapIds.isNotEmpty())
            assertEquals(entry.hapIds.distinct().sorted(), entry.hapIds)
        }

        val byPosAllele = entries.associateBy { Triple(it.contig, it.position, it.allele) }

        //A -> T,<DEL>: LineA and LineA2 share the same haplotype checksum in this range, so the
        //"T" row must contain exactly one hapId even though two samples carry it.
        val delRow = byPosAllele.getValue(Triple("1", 3131, "<DEL>"))
        assertEquals(listOf("8967fabf10e55d881caa6fe192e7d4ca"), delRow.hapIds)
        val tRow = byPosAllele.getValue(Triple("1", 3131, "T"))
        assertEquals(listOf("3149b3144f93134eb29661bade697fc6"), tRow.hapIds)

        //T -> A,<DEL>: allele "A" is carried only by LineA/LineA2, neither of which has a
        //haplotype in 2:12001-16500 (only LineB does) - so allele "A" must be entirely absent.
        assertFalse(byPosAllele.containsKey(Triple("2", 12019, "A")))
        val del2Row = byPosAllele.getValue(Triple("2", 12019, "<DEL>"))
        assertEquals(listOf("6fb2de47c835bd9ab026c02d62f49807"), del2Row.hapIds)
    }

    //endregion
}
