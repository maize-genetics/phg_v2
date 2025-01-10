package net.maizegenetics.phgv2.utils

import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.cli.TestExtension
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertTrue

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
        assertEquals(2, firstRefRange?.size) // entries for both LineA and LineB
        // verify the firstRefRange contains an entry for SampleName=LineA, ID=12f0cec9102e84a161866e37072443b7
        // Define the expected map entry
        val expectedEntryLineA = mapOf("SampleName" to "LineA", "ID" to "12f0cec9102e84a161866e37072443b7")
        val expectedEntryLineB = mapOf("SampleName" to "LineB", "ID" to "4fc7b8af32ddd74e07cb49d147ef1938")

        // Assert that firstRefRange contains the expected entry
        assertTrue(firstRefRange?.any { it == expectedEntryLineA } ?: false,
            "Expected entry LineA not found in firstRefRange")
        assertTrue(firstRefRange?.any { it == expectedEntryLineB } ?: false,
            "Expected entry LineB not found in firstRefRange")

        // delete the tiledbArray so next tests can recreate what they need
        File(TiledbCoreHvcfUtilsTest.dbPath).deleteRecursively()

    }

    @Test
    fun testQueryIdsByRefRange() {

        // testing output from parseTiledbAltHeaders
        // The values used for the asserts were copied from the values in the hvcf files
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
        assertEquals(1, resultsByRefRange.keys.size)
        val rangeResults = resultsByRefRange.get("1:1-1000")
        assertEquals(2, rangeResults?.size)
        // verify rangeResults contains both 12f0cec9102e84a161866e37072443b7 and 4fc7b8af32ddd74e07cb49d147ef1938
        assertTrue(rangeResults?.contains("12f0cec9102e84a161866e37072443b7") ?: false)
        assertTrue(rangeResults?.contains("4fc7b8af32ddd74e07cb49d147ef1938") ?: false)

        val resultsByRefRange2 = queryIDsByRefRange(TiledbCoreHvcfUtilsTest.altHeaderArray, listOf("1:1-1000","2:12001-16500"))
        println("Results with 2 refRanges: $resultsByRefRange2")
        assertEquals(2, resultsByRefRange2.keys.size)
        var firstRangeResults = resultsByRefRange2.get("1:1-1000")
        assertEquals(2, firstRangeResults?.size)
        // verify firstRangeResults contains both 12f0cec9102e84a161866e37072443b7 and 4fc7b8af32ddd74e07cb49d147ef1938
        assertTrue(firstRangeResults?.contains("12f0cec9102e84a161866e37072443b7") ?: false)
        assertTrue(firstRangeResults?.contains("4fc7b8af32ddd74e07cb49d147ef1938") ?: false)

        var secondRangeResults = resultsByRefRange2.get("2:12001-16500")
        assertEquals(2, secondRangeResults?.size)
        // verify secondRangeResults contains both 184a72815a2ba5949635cc38769cedd0 and 6fb2de47c835bd9ab026c02d62f49807
        assertTrue(secondRangeResults?.contains("184a72815a2ba5949635cc38769cedd0") ?: false)
        assertTrue(secondRangeResults?.contains("6fb2de47c835bd9ab026c02d62f49807") ?: false)


        val resultsByRefRange5Ranges = queryIDsByRefRange(TiledbCoreHvcfUtilsTest.altHeaderArray, listOf("1:1-1000","1:12001-16500","1:33001-34000","2:1-1000","2:5501-6500"))
        println("Results with 5 refRanges: $resultsByRefRange5Ranges")

        assertEquals(5, resultsByRefRange5Ranges.keys.size)
        firstRangeResults = resultsByRefRange5Ranges.get("1:1-1000")
        assertEquals(2, firstRangeResults?.size)
        // verify firstRangeResults contains both 12f0cec9102e84a161866e37072443b7 and 4fc7b8af32ddd74e07cb49d147ef1938
        assertTrue(firstRangeResults?.contains("12f0cec9102e84a161866e37072443b7") ?: false)
        assertTrue(firstRangeResults?.contains("4fc7b8af32ddd74e07cb49d147ef1938") ?: false)

        secondRangeResults = resultsByRefRange5Ranges.get("1:12001-16500")
        assertEquals(2, secondRangeResults?.size)
        // verify secondRangeResults contains both d4c8b5505d7046b41d7f69b246063ebb and aff71f19de448514a6d9208b1fcb4e8a
        assertTrue(secondRangeResults?.contains("d4c8b5505d7046b41d7f69b246063ebb") ?: false)
        assertTrue(secondRangeResults?.contains("aff71f19de448514a6d9208b1fcb4e8a") ?: false)

        val thirdRangeResults = resultsByRefRange5Ranges.get("1:33001-34000")
        assertEquals(2, thirdRangeResults?.size)
        // verify thirdRangeResults contains both 5a2ff3fb844d5647987da5c194d1c728 and 79146831745c85d26117f13b5873935f
        assertTrue(thirdRangeResults?.contains("5a2ff3fb844d5647987da5c194d1c728") ?: false)
        assertTrue(thirdRangeResults?.contains("79146831745c85d26117f13b5873935f") ?: false)

        val fourthRangeResults = resultsByRefRange5Ranges.get("2:1-1000")
        assertEquals(2, fourthRangeResults?.size)
        // verify fourthRangeResults contains both 13417ecbb38b9a159e3ca8c9dade7088 and 180417a01edbfed525d7c238910e0ff4
        assertTrue(fourthRangeResults?.contains("13417ecbb38b9a159e3ca8c9dade7088") ?: false)
        assertTrue(fourthRangeResults?.contains("180417a01edbfed525d7c238910e0ff4") ?: false)

        val fifthRangeResults = resultsByRefRange5Ranges.get("2:5501-6500")
        assertEquals(2, fifthRangeResults?.size)
        // verify fifthRangeResults contains both 50044914d5111c5b5ec58c9d720e3b2d and 45b121547c7ae517a181fdd2621495c4
        assertTrue(fifthRangeResults?.contains("50044914d5111c5b5ec58c9d720e3b2d") ?: false)
        assertTrue(fifthRangeResults?.contains("45b121547c7ae517a181fdd2621495c4") ?: false)

        // delete the tiledbArray so next tests can recreate what they need
        File(TiledbCoreHvcfUtilsTest.dbPath).deleteRecursively()

    }

}