package net.maizegenetics.phgv2.utils

import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.cli.TestExtension
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File

@ExtendWith(TestExtension::class)
class TiledbAltHeadersQueriesTest {
    private val myLogger = LogManager.getLogger(TiledbAltHeadersQueriesTest::class.java)

    companion object {

        private val multiInputDir = "${TestExtension.tempDir}/multi-input/"
        private val outputHvcfDir = "${TestExtension.tempDir}/output/"
        private val imputedHvcfDir = "${TestExtension.tempDir}/imputed/"
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
            File(imputedHvcfDir).mkdirs()

            File(lineAhvcf).copyTo(File(multiInputDir + File(lineAhvcf).name))
            File(lineBhvcf).copyTo(File(multiInputDir + File(lineBhvcf).name))
            File(imputedHvcf).copyTo(File(imputedHvcfDir + File(imputedHvcf).name))
        }

        @AfterAll
        @JvmStatic
        fun teardown() {
            File(outputHvcfDir).deleteRecursively()
            File(multiInputDir).deleteRecursively()
            File(imputedHvcfDir).deleteRecursively()
            // delete the tempDir
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testQuerySampleNameID_byRefRange() {
        // testing output from parseTiledbAltHeaders
        var vcfReader = VCFFileReader(File(lineAhvcf), false)
        val altHeadersLineA = parseTiledbAltHeaders(vcfReader)
        println("Finished parsing lineAhvcf")

        vcfReader = VCFFileReader(File(lineBhvcf), false)
        val altHeadersLineB = parseTiledbAltHeaders(vcfReader)
        println("Finished parsing lineBhvcf")

        // Create the tiledb array by calling the function createTileDBArray
        println("Creating tiledb array")
        // ensure the top folder exists
        File(dbPath).mkdirs()
        createTileDBCoreArrays(dbPath)

        // Write the altHeaders to the tiledb array
        println("Writing LineA altHeaders to tiledb array")
        writeAltDataToTileDB(altHeaderArray, altHeadersLineA)
        println("Writing LineB altHeaders to tiledb array")
        writeAltDataToTileDB(altHeaderArray, altHeadersLineB)

        // To test what was written, we must extract data
        println("TestCase: Extracting data from tiledb array")

        val rangesToQuery = listOf("1:1-1000", "1:6501-11000", "2:1-1000", "2:5501-6500")

        // Try with dynamic buffers
        val streamingResults = queryWithStreaming_sampleNameIdByRefRange(TiledbCoreHvcfUtilsTest.altHeaderArray, rangesToQuery)

        println("\nStreaming buffer results:")
        streamingResults.forEach { println(it) }

        Assertions.assertTrue(streamingResults.size == 4)
        val firstRefRange = streamingResults.get("1:1-1000")
        Assertions.assertEquals(2, firstRefRange?.size) // entries for both LineA and LineB
        // verify the firstRefRange contains an entry for SampleName=LineA, ID=12f0cec9102e84a161866e37072443b7
        // Define the expected map entry
        val expectedEntryLineA = mapOf("SampleName" to "LineA", "ID" to "12f0cec9102e84a161866e37072443b7")
        val expectedEntryLineB = mapOf("SampleName" to "LineB", "ID" to "4fc7b8af32ddd74e07cb49d147ef1938")

        // Assert that firstRefRange contains the expected entry
        Assertions.assertTrue(firstRefRange?.any { it == expectedEntryLineA } ?: false,
            "Expected entry LineA not found in firstRefRange")
        Assertions.assertTrue(firstRefRange?.any { it == expectedEntryLineB } ?: false,
            "Expected entry LineB not found in firstRefRange")

        // delete the tiledbArray so next tests can recreate what they need
        File(TiledbCoreHvcfUtilsTest.dbPath).deleteRecursively()

    }

    @Test
    fun testQueryIdsByRefRange() {

        // testing output from parseTiledbAltHeaders
        var vcfReader = VCFFileReader(File(lineAhvcf), false)
        val altHeadersLineA = parseTiledbAltHeaders(vcfReader)
        println("Finished parsing LineA alt headers ")
        vcfReader = VCFFileReader(File(lineBhvcf), false)
        val altHeadersLineB = parseTiledbAltHeaders(vcfReader)
        println("Finished parsing lineBhvcf alt headers")


        // ensure the top folder exists
        File(dbPath).mkdirs()
        // Create the tiledb array by calling the function createTileDBArray
        println("Creating tiledb array")
        //createTileDBArraySingleDimensionID(tiledbArrayName)
        createTileDBCoreArrays(dbPath)

        // Write the altHeaders to the tiledb array
        println("Writing altHeaders to tiledb array")
        writeAltDataToTileDB(altHeaderArray, altHeadersLineA)
        writeAltDataToTileDB(altHeaderArray, altHeadersLineB)
        // Try to query by refRange
        val resultsByRefRange = queryIDsByRefRange(TiledbCoreHvcfUtilsTest.altHeaderArray, listOf("1:1-1000"))
        println("Results with 1 refRange: $resultsByRefRange")
        val resultsByRefRange2 = queryIDsByRefRange(TiledbCoreHvcfUtilsTest.altHeaderArray, listOf("1:1-1000","2:12001-16500"))
        println("Results with 2 refRanges: $resultsByRefRange2")

        val resultsByRefRange5Ranges = queryIDsByRefRange(TiledbCoreHvcfUtilsTest.altHeaderArray, listOf("1:1-1000","1:12001-16500","1:33001-34000","2:1-1000","2:5501-6500"))
        println("Results with 5 refRanges: $resultsByRefRange5Ranges")

        // Need some asserts here !!

        // delete the tiledbArray so next tests can recreate what they need
        File(TiledbCoreHvcfUtilsTest.dbPath).deleteRecursively()

    }

}