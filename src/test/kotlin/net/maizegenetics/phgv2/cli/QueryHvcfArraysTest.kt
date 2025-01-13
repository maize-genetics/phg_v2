package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class QueryHvcfArraysTest {
    companion object {

        private val multiInputDir = "${TestExtension.tempDir}/multi-input/"
        private val outputHvcfDir = "${TestExtension.tempDir}/output/"
        private val outputQueryResultsDir = "${TestExtension.tempDir}/queryResults/"
        val lineAhvcf = "data/test/tiledbCoreHvcf/LineA.h.vcf"
        val lineBhvcf = "data/test/tiledbCoreHvcf/LineB.h.vcf"
        val imputedHvcf = "data/test/resequenceHaplotypeVCF/Imputation.h.vcf"
        val dbPath = TestExtension.testTileDBURI
        val altHeaderArray = dbPath + "/alt_header_array"
        val variantsArray = dbPath + "/hvcf_variants_array"

        @BeforeAll
        @JvmStatic
        fun setup() {
            File(multiInputDir).mkdirs()
            File(outputHvcfDir).mkdirs()
            File(dbPath).mkdirs()
            File(outputQueryResultsDir).mkdirs()

            File(lineAhvcf).copyTo(File(multiInputDir + File(lineAhvcf).name))
            File(lineBhvcf).copyTo(File(multiInputDir + File(lineBhvcf).name))
            File(imputedHvcf).copyTo(File(multiInputDir + File(imputedHvcf).name))
        }

        @AfterAll
        @JvmStatic
        fun teardown() {
            File(outputHvcfDir).deleteRecursively()
            File(multiInputDir).deleteRecursively()
            File(outputQueryResultsDir).deleteRecursively()
            // delete the tempDir
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testCliktParams() {
        val queryHvcfArrays = QueryHvcfArrays()

        // Test missing output-file parameter
        val resultMissingOutputFile = queryHvcfArrays.test("--db-path ${TestExtension.testTileDBURI} --query-type DISTINCT_SAMPLES")
        assertEquals(resultMissingOutputFile.statusCode, 1)
        assertEquals("Usage: query-hvcf-arrays [<options>]\n" +
                "\n" +
                "Error: missing option --output-file\n",resultMissingOutputFile.output)

        // Test missing query-type parameter
        val resultMissingQueryType = queryHvcfArrays.test("--db-path ${TestExtension.testTileDBURI} --output-file dummyOutFile")
        assertEquals(resultMissingQueryType.statusCode, 1)
        assertEquals("Usage: query-hvcf-arrays [<options>]\n" +
                "\n" +
                "Error: missing option --query-type\n",resultMissingQueryType.output)

    }

    @Test
    fun testQueryHvcfArays() {
        // First create the arrays:
        val initHvcfArray = InitHvcfArray()
        val resultInit = initHvcfArray.test("--db-path ${TestExtension.testTileDBURI}")
        assertEquals(resultInit.statusCode, 0)

        // Now load the hvcf files
        val loadHVCF = LoadHvcf()
        val resultLoad = loadHVCF.test("--hvcf-dir $multiInputDir --db-path ${TestExtension.testTileDBURI}")
        assertEquals(resultLoad.statusCode, 0)

        // Now query the arrays
        val queryHvcfArrays = QueryHvcfArrays()
        var resultQuery = queryHvcfArrays.test("--db-path ${TestExtension.testTileDBURI} --query-type DISTINCT_SAMPLES --output-file $outputQueryResultsDir/sampleNames.txt")
        assertEquals(resultQuery.statusCode, 0)

        resultQuery = queryHvcfArrays.test("--db-path ${TestExtension.testTileDBURI} --query-type DISTINCT_RANGES --output-file $outputQueryResultsDir/ranges.txt")
        assertEquals(resultQuery.statusCode, 0)

        // Read the results files and compare to expected

        val actualSampleNames = File("$outputQueryResultsDir/sampleNames.txt").readText()
        // print the sample names
        println("SampleNames file contents")
        println(actualSampleNames)
        assertEquals("LineA\nLineB\nTestLine2", actualSampleNames.trim())
        assertEquals(3, actualSampleNames.trim().split("\n").size)

        // query the sample names from the altHeaderArray - should be 1 less
        resultQuery = queryHvcfArrays.test("--db-path ${TestExtension.testTileDBURI} --query-type DISTINCT_SAMPLES --output-file $outputQueryResultsDir/sampleNamesAltHeader.txt --array-type ALT_HEADER")
        assertEquals(resultQuery.statusCode, 0)
        val actualSampleNamesAltHeader = File("$outputQueryResultsDir/sampleNamesAltHeader.txt").readText()
        // print the sample names
        println("SampleNamesAltHeader file contents")
        println(actualSampleNamesAltHeader)
        assertEquals("LineA\nLineB", actualSampleNamesAltHeader.trim())
        assertEquals(2, actualSampleNamesAltHeader.trim().split("\n").size)

        // Read the actual refRanges
        val actualRefRanges = File("$outputQueryResultsDir/ranges.txt").readText()
        // print the refRanges
        println("RefRanges file contents:")
        println(actualRefRanges)
        assertEquals(38, actualRefRanges.trim().split("\n").size)

    }

}