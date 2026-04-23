package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import java.nio.file.Files
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class MappingCountTableQcTest {

    companion object {
        val tempTestDir = "${TestExtension.tempDir}mappingCountTableQcTest/"
        val hvcfDir = "${tempTestDir}hvcfDir/"
        val readMappingDir = "${tempTestDir}readMappingDir/"
        val outputDir = "${tempTestDir}outputDir/"

        @JvmStatic
        @BeforeAll
        fun setup() {
            resetDirs()
            setupDebugLogging()

            // Copy smallseq hvcf files into the test hvcf directory
            Files.copy(File(TestExtension.smallseqLineAHvcfFile).toPath(), File("${hvcfDir}LineA.h.vcf").toPath())
            Files.copy(File(TestExtension.smallseqLineBHvcfFile).toPath(), File("${hvcfDir}LineB.h.vcf").toPath())
            Files.copy(File(TestExtension.smallseqRefHvcfFile).toPath(), File("${hvcfDir}Ref.h.vcf").toPath())

            // Copy the ropebwt read mapping files into the test read mapping directory
            Files.copy(
                File("data/test/ropebwt/LineA_1_readMapping.txt").toPath(),
                File("${readMappingDir}LineA_1_readMapping.txt").toPath()
            )
            Files.copy(
                File("data/test/ropebwt/LineA_2_readMapping.txt").toPath(),
                File("${readMappingDir}LineA_2_readMapping.txt").toPath()
            )
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            resetDirs()
        }

        private fun resetDirs() {
            File(tempTestDir).deleteRecursively()
            File(tempTestDir).mkdirs()
            File(hvcfDir).mkdirs()
            File(readMappingDir).mkdirs()
            File(outputDir).mkdirs()
        }
    }

    // -------------------------------------------------------------------------
    // CLI parameter validation
    // -------------------------------------------------------------------------

    @Test
    fun testCliktMissingHvcfDir() {
        val result = MappingCountTableQc().test(
            "--read-mapping-dir $readMappingDir --output-dir $outputDir"
        )
        assertEquals(1, result.statusCode)
        assertTrue(result.stderr.contains("missing option --hvcf-dir"))
    }

    @Test
    fun testCliktMissingReadMappingDir() {
        val result = MappingCountTableQc().test(
            "--hvcf-dir $hvcfDir --output-dir $outputDir"
        )
        assertEquals(1, result.statusCode)
        assertTrue(result.stderr.contains("missing option --read-mapping-dir"))
    }

    @Test
    fun testCliktMissingOutputDir() {
        val result = MappingCountTableQc().test(
            "--hvcf-dir $hvcfDir --read-mapping-dir $readMappingDir"
        )
        assertEquals(1, result.statusCode)
        assertTrue(result.stderr.contains("missing option --output-dir"))
    }

    // -------------------------------------------------------------------------
    // buildCountOutputFile
    // -------------------------------------------------------------------------

    @Test
    fun testBuildCountOutputFile() {
        val qc = MappingCountTableQc()
        val outDir = "/some/output/dir"
        val inputFile = File("/some/input/dir/MySample_readMapping.txt")

        val result = qc.buildCountOutputFile(outDir, inputFile)

        assertEquals(File(outDir, "MySample_mappingCounts.txt"), result)
    }

    @Test
    fun testBuildCountOutputFilePreservesPrefix() {
        val qc = MappingCountTableQc()
        val outDir = "/out"
        val inputFile = File("/reads/LineA_1_readMapping.txt")

        val result = qc.buildCountOutputFile(outDir, inputFile)

        assertEquals(File("/out", "LineA_1_mappingCounts.txt"), result)
    }

    // -------------------------------------------------------------------------
    // findRefRangeForHapIdSet
    // -------------------------------------------------------------------------

    @Test
    fun testFindRefRangeForHapIdSetSingleHap() {
        val qc = MappingCountTableQc()
        val range1 = ReferenceRange("1", 1, 1000)
        val range2 = ReferenceRange("1", 1001, 5500)
        val hapIdToRefRange = mapOf(
            "hapA" to listOf(range1),
            "hapB" to listOf(range2)
        )

        val result = qc.findRefRangeForHapIdSet(listOf("hapA"), hapIdToRefRange)
        assertEquals(range1, result)
    }

    @Test
    fun testFindRefRangeForHapIdSetMultipleHapsAgreeOnRange() {
        val qc = MappingCountTableQc()
        val range1 = ReferenceRange("1", 1, 1000)
        val hapIdToRefRange = mapOf(
            "hapA" to listOf(range1),
            "hapB" to listOf(range1)
        )

        val result = qc.findRefRangeForHapIdSet(listOf("hapA", "hapB"), hapIdToRefRange)
        assertEquals(range1, result)
    }

    @Test
    fun testFindRefRangeForHapIdSetDisagreementReturnsUnknown() {
        val qc = MappingCountTableQc()
        val range1 = ReferenceRange("1", 1, 1000)
        val range2 = ReferenceRange("1", 1001, 5500)
        val hapIdToRefRange = mapOf(
            "hapA" to listOf(range1),
            "hapB" to listOf(range2)
        )

        val result = qc.findRefRangeForHapIdSet(listOf("hapA", "hapB"), hapIdToRefRange)
        assertEquals("UNKNOWN", result.contig)
    }

    // -------------------------------------------------------------------------
    // convertReadMappingToHapIdCounts
    // -------------------------------------------------------------------------

    @Test
    fun testConvertReadMappingToHapIdCountsSingleMappings() {
        val qc = MappingCountTableQc()
        val graph = HaplotypeGraph(hvcfDir)
        val hapIdAndRangeToSampleGametes = graph.hapIdAndRefRangeToSampleGametes()
        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()

        // Use a single unambiguous LineA hap id from range 1:1-1000
        // hapId from LineA.h.vcf: 12f0cec9102e84a161866e37072443b7 -> RefRange=1:1-1000
        val hapId = "12f0cec9102e84a161866e37072443b7"
        val readMappings = mapOf(listOf(hapId) to 100)

        val result = qc.convertReadMappingToHapIdCounts(readMappings, hapIdAndRangeToSampleGametes, hapIdToRefRangeMap)

        // Should have an entry for LineA:0 at range 1:1-1000 with count 100
        val range = ReferenceRange("1", 1, 1000)
        val lineAGamete = SampleGamete("LineA", 0)
        assertTrue(result.containsKey(Pair(range, lineAGamete)),
            "Expected count for LineA at 1:1-1000")
        assertEquals(100, result[Pair(range, lineAGamete)])
    }

    @Test
    fun testConvertReadMappingToHapIdCountsAccumulatesAcrossEntries() {
        val qc = MappingCountTableQc()
        val graph = HaplotypeGraph(hvcfDir)
        val hapIdAndRangeToSampleGametes = graph.hapIdAndRefRangeToSampleGametes()
        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()

        // Two entries both pointing to LineA's hap at 1:1-1000
        val hapId = "12f0cec9102e84a161866e37072443b7"
        val readMappings = mapOf(
            listOf(hapId) to 200,
            listOf(hapId) to 300   // duplicate key — map collapses to last value (300)
        )

        val result = qc.convertReadMappingToHapIdCounts(readMappings, hapIdAndRangeToSampleGametes, hapIdToRefRangeMap)
        val range = ReferenceRange("1", 1, 1000)
        val lineAGamete = SampleGamete("LineA", 0)
        // Only last entry survives Map deduplication, so count = 300
        assertEquals(300, result[Pair(range, lineAGamete)])
    }

    @Test
    fun testConvertReadMappingSkipsUnknownRefRange() {
        val qc = MappingCountTableQc()
        val graph = HaplotypeGraph(hvcfDir)
        val hapIdAndRangeToSampleGametes = graph.hapIdAndRefRangeToSampleGametes()
        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()

        // Two hapIds from different ranges — findRefRangeForHapIdSet should return UNKNOWN
        val hapIdRange1 = "12f0cec9102e84a161866e37072443b7" // 1:1-1000
        val hapIdRange2 = "3149b3144f93134eb29661bade697fc6" // 1:1001-5500
        val readMappings = mapOf(listOf(hapIdRange1, hapIdRange2) to 50)

        val result = qc.convertReadMappingToHapIdCounts(readMappings, hapIdAndRangeToSampleGametes, hapIdToRefRangeMap)

        // The ambiguous set should be skipped; result should be empty
        assertTrue(result.isEmpty(), "Ambiguous hapId set should produce no counts, got: $result")
    }

    // -------------------------------------------------------------------------
    // writeOutCounts
    // -------------------------------------------------------------------------

    @Test
    fun testWriteOutCountsHeaderAndStructure() {
        val qc = MappingCountTableQc()
        val graph = HaplotypeGraph(hvcfDir)
        val outputFile = File("${outputDir}writeOutCountsTest.txt")

        // Build a minimal counts map: give LineA:0 a count of 42 at 1:1-1000
        val range = ReferenceRange("1", 1, 1000)
        val lineAGamete = SampleGamete("LineA", 0)
        val hapIdCounts = mapOf(Pair(range, lineAGamete) to 42)

        qc.writeOutCounts(hapIdCounts, graph, outputFile)

        assertTrue(outputFile.exists())
        val lines = bufferedReader(outputFile.absolutePath).readLines()

        // Header line
        val header = lines[0]
        assertTrue(header.startsWith("RefRange\t"), "Header should start with RefRange, got: $header")
        // All sample gametes from the graph should appear in the header
        val sampleGametes = graph.sampleGametesInGraph()
        sampleGametes.forEach { sg ->
            assertTrue(header.contains(sg.toString()), "Header should contain $sg")
        }

        // Every range in the graph should have an output row
        val ranges = graph.ranges()
        assertEquals(ranges.size + 1, lines.size, "Expected header + one row per range")
    }

    @Test
    fun testWriteOutCountsValuesCorrect() {
        val qc = MappingCountTableQc()
        val graph = HaplotypeGraph(hvcfDir)
        val outputFile = File("${outputDir}writeOutCountsValues.txt")

        val range = ReferenceRange("1", 1, 1000)
        val sampleGametes = graph.sampleGametesInGraph().toList()
        // Give the first sample gamete count=7 at range 1:1-1000
        val sg0 = sampleGametes[0]
        val hapIdCounts = mapOf(Pair(range, sg0) to 7)

        qc.writeOutCounts(hapIdCounts, graph, outputFile)

        val lines = bufferedReader(outputFile.absolutePath).readLines()
        val headerCols = lines[0].split("\t")

        // Find the row for range 1:1-1000
        val rangeRow = lines.drop(1).firstOrNull { it.startsWith("1:1-1000\t") }
            ?: error("Row for 1:1-1000 not found")
        val rowCols = rangeRow.split("\t")

        val sg0ColIndex = headerCols.indexOf(sg0.toString())
        assertTrue(sg0ColIndex >= 1, "Column for $sg0 not found in header")
        assertEquals("7", rowCols[sg0ColIndex], "Count for $sg0 at 1:1-1000 should be 7")

        // All other sample gametes for this range should have count 0
        sampleGametes.filter { it != sg0 }.forEach { sg ->
            val colIdx = headerCols.indexOf(sg.toString())
            if (colIdx >= 1) {
                assertEquals("0", rowCols[colIdx], "Count for $sg at 1:1-1000 should be 0")
            }
        }
    }

    // -------------------------------------------------------------------------
    // buildMappingTables / processSingleReadMappingFile (integration)
    // -------------------------------------------------------------------------

    @Test
    fun testBuildMappingTablesProducesOutputFiles() {
        val qc = MappingCountTableQc()
        val singleReadMappingDir = File("${tempTestDir}singleReadMappingDir/")
        val singleOutputDir = File("${tempTestDir}singleOutput/")
        singleReadMappingDir.mkdirs()
        singleOutputDir.mkdirs()

        Files.copy(
            File("data/test/ropebwt/LineA_1_readMapping.txt").toPath(),
            File(singleReadMappingDir, "LineA_1_readMapping.txt").toPath()
        )

        qc.buildMappingTables(hvcfDir, singleReadMappingDir, singleOutputDir)

        val outputFiles = singleOutputDir.listFiles()?.toList() ?: emptyList()
        assertEquals(1, outputFiles.size, "Expected exactly one output file")
        assertEquals("LineA_1_mappingCounts.txt", outputFiles[0].name)
    }

    @Test
    fun testBuildMappingTablesProcessesAllReadMappingFiles() {
        val qc = MappingCountTableQc()
        val multiOutputDir = File("${tempTestDir}multiOutput/")
        multiOutputDir.mkdirs()

        // readMappingDir contains both LineA_1 and LineA_2 files (set up in @BeforeAll)
        qc.buildMappingTables(hvcfDir, File(readMappingDir), multiOutputDir)

        val outputFiles = multiOutputDir.listFiles()?.map { it.name }?.sorted() ?: emptyList()
        assertEquals(2, outputFiles.size, "Expected two output files for two read mapping files")
        assertTrue(outputFiles.contains("LineA_1_mappingCounts.txt"))
        assertTrue(outputFiles.contains("LineA_2_mappingCounts.txt"))
    }

    @Test
    fun testBuildMappingTablesOutputContainsExpectedStructure() {
        val qc = MappingCountTableQc()
        val structOutputDir = File("${tempTestDir}structOutput/")
        structOutputDir.mkdirs()

        Files.copy(
            File("data/test/ropebwt/LineA_1_readMapping.txt").toPath(),
            File("${tempTestDir}structReadMapping/LineA_1_readMapping.txt").also { it.parentFile.mkdirs() }.toPath(),
            java.nio.file.StandardCopyOption.REPLACE_EXISTING
        )

        qc.buildMappingTables(hvcfDir, File("${tempTestDir}structReadMapping"), structOutputDir)

        val outputFile = File(structOutputDir, "LineA_1_mappingCounts.txt")
        assertTrue(outputFile.exists())

        val lines = bufferedReader(outputFile.absolutePath).readLines()
        assertTrue(lines.isNotEmpty())

        // Header should start with "RefRange"
        assertTrue(lines[0].startsWith("RefRange\t"), "First line should be header, got: ${lines[0]}")

        // Each subsequent line should parse as <range>\t<tab-separated integers>
        val graph = HaplotypeGraph(hvcfDir)
        val numSampleGametes = graph.sampleGametesInGraph().size
        lines.drop(1).forEach { line ->
            val cols = line.split("\t")
            assertEquals(numSampleGametes + 1, cols.size,
                "Row should have RefRange + $numSampleGametes count columns: $line")
            // All count columns should be parseable as integers
            cols.drop(1).forEach { countStr ->
                assertTrue(countStr.toIntOrNull() != null,
                    "Count column should be integer, got: '$countStr' in line: $line")
            }
        }
    }

    @Test
    fun testBuildMappingTablesIgnoresNonReadMappingFiles() {
        val qc = MappingCountTableQc()
        val filteredReadMappingDir = File("${tempTestDir}filteredReadMappingDir/")
        val filteredOutputDir = File("${tempTestDir}filteredOutput/")
        filteredReadMappingDir.mkdirs()
        filteredOutputDir.mkdirs()

        // Copy a valid read mapping file
        Files.copy(
            File("data/test/ropebwt/LineA_1_readMapping.txt").toPath(),
            File(filteredReadMappingDir, "LineA_1_readMapping.txt").toPath()
        )
        // Create a file that does NOT end in _readMapping.txt — should be ignored
        File(filteredReadMappingDir, "someOtherFile.txt").writeText("not a read mapping file")

        qc.buildMappingTables(hvcfDir, filteredReadMappingDir, filteredOutputDir)

        val outputFiles = filteredOutputDir.listFiles()?.toList() ?: emptyList()
        assertEquals(1, outputFiles.size, "Only the valid read mapping file should produce output")
        assertEquals("LineA_1_mappingCounts.txt", outputFiles[0].name)
    }
}
