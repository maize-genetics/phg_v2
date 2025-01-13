package net.maizegenetics.phgv2.utils

import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.cli.TestExtension
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.assertTrue
import org.jetbrains.kotlinx.dataframe.api.*

/**
 * This Class contains tests for the queries in TiledbCoreVariantQueries.kt
 */
class TileDBCoreVariantQueriesTest {
    private val myLogger = LogManager.getLogger(TileDBCoreVariantQueriesTest::class.java)

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
    fun testQueryByRefRange_variantArray() {

        // filter on specific refRanges and get it correct for the variants array,
        // Make sure the top folder exists
        File(dbPath).mkdirs()

        // testing output from parseTiledbAltHeaders
        var vcfReader = VCFFileReader(File(lineAhvcf), false)

        val variantDataLineA = TiledbCoreHvcfUtils.parseTiledbVariantData(vcfReader)
        println("Finished parsing lineAhvcf variant data")

        vcfReader = VCFFileReader(File(lineBhvcf), false)
        val variantDataLineB = TiledbCoreHvcfUtils.parseTiledbVariantData(vcfReader)
        println("Finished parsing lineBhvcf variant data")

        // Create the tiledb array by calling the function createTileDBArray
        println("Creating tiledb array")
        //createTileDBArraySingleDimensionID(tiledbArrayName)
        TiledbCoreHvcfUtils.createTileDBCoreArrays(dbPath)

        //Write to the variants array
        TiledbCoreHvcfUtils.writeVariantsDataToTileDB(variantsArray, variantDataLineA)
        TiledbCoreHvcfUtils.writeVariantsDataToTileDB(variantsArray, variantDataLineB)

        // Load the imputed hvcf file

        vcfReader = VCFFileReader(File(imputedHvcf), false)
        val altHeadersImputed = TiledbCoreHvcfUtils.parseTiledbAltHeaders(vcfReader)
        println("Finished parsing imputedHvcf alt headers ")
        val variantData = TiledbCoreHvcfUtils.parseTiledbVariantData(vcfReader)
        println("Finished parsing imputedHvcf variant data")

        // Now load the data
        val arrayName = "${dbPath}/alt_header_array"
        TiledbCoreHvcfUtils.writeAltDataToTileDB(arrayName, altHeadersImputed)
        val variantArrayName = "${dbPath}/hvcf_variants_array"
        TiledbCoreHvcfUtils.writeVariantsDataToTileDB(variantArrayName, variantData)

        // Query the variants array
        val results = queryVariantArray_byRefRange(variantsArray, listOf("1:1-1000","1:12001-16500","1:33001-34000","2:1-1000","2:5501-6500"))
        println("\nResults on variantsArray:")
        results.forEach { println(it) }

        // Manually verifying - need asserts here !!
        assertEquals(5, results.size)
        val firstRefRange = results.get("1:1-1000")
        assertEquals(3, firstRefRange?.size) // entries for LineA, LineB, and TestLine2
        // verify the firstRefRange contains an entry for SampleName=LineA, ID=12f0cec9102e84a161866e37072443b7
        // Define the expected map entry
        var expectedEntryLineA = mapOf("SampleName" to "LineA", "ID1" to "12f0cec9102e84a161866e37072443b7","ID2" to "12f0cec9102e84a161866e37072443b7")
        var expectedEntryLineB = mapOf("SampleName" to "LineB", "ID1" to "4fc7b8af32ddd74e07cb49d147ef1938", "ID2" to "4fc7b8af32ddd74e07cb49d147ef1938")
        // Many of the imputed valuesare the same as LineA
        var expectedEntryTestLine2 = mapOf("SampleName" to "TestLine2", "ID1" to "12f0cec9102e84a161866e37072443b7", "ID2" to "12f0cec9102e84a161866e37072443b7")
        // Assert that firstRefRange contains the expected entry
        assertTrue(firstRefRange?.any { it == expectedEntryLineA } ?: false,
            "Expected entry LineA not found in firstRefRange")
        assertTrue(firstRefRange?.any { it == expectedEntryLineB } ?: false,
            "Expected entry LineB not found in firstRefRange")
        assertTrue(firstRefRange?.any { it == expectedEntryTestLine2 } ?: false,
            "Expected entry TestLine2 not found in firstRefRange")

        val lastRefRange = results.get("2:5501-6500")
        assertEquals(3, lastRefRange?.size) // entries for LineA, LineB, and TestLine2

        // NOTE _ these are haploid, so ID1 and ID2 are the same
        expectedEntryLineA = mapOf("SampleName" to "LineA", "ID1" to "50044914d5111c5b5ec58c9d720e3b2d","ID2" to "50044914d5111c5b5ec58c9d720e3b2d")
        expectedEntryLineB = mapOf("SampleName" to "LineB", "ID1" to "45b121547c7ae517a181fdd2621495c4", "ID2" to "45b121547c7ae517a181fdd2621495c4")
        // Many of the imputed valuesare the same as LineA
        expectedEntryTestLine2 = mapOf("SampleName" to "TestLine2", "ID1" to "50044914d5111c5b5ec58c9d720e3b2d", "ID2" to "50044914d5111c5b5ec58c9d720e3b2d")
        // Assert that firstRefRange contains the expected entry

        assertTrue(lastRefRange?.any { it == expectedEntryLineA } ?: false,
            "Expected entry LineA not found in lastRefRange")
        assertTrue(lastRefRange?.any { it == expectedEntryLineB } ?: false,
            "Expected entry LineB not found in lastRefRange")
        assertTrue(lastRefRange?.any { it == expectedEntryTestLine2 } ?: false,
            "Expected entry TestLine2 not found in lastRefRange")


        // convert to a dataFrame
        // NOTE: I"M not sure this is the format we want for the dataFrame -may need to
        // make adjustments in the convertQueryResultToDataFrame function to better split
        // this data into columns.  It is here to showcase the results and for other
        // developers to determine if this is useful and if so, what needs to be altered.
        val dataFrame = convertQueryResultToDataFrame(results)
        println("as dataFrame:")
        dataFrame.print()

        // delete the tiledbArray so next tests can recreate what they need
        File(dbPath).deleteRecursively()

    }

    @Test
    fun testDistinctValues() {

        // This test loads data from the imputed hvcf file and the LineA and LineB hvcf files
        // then queries the loaded db for distinct SampleNames and distinct refRanges
        println("Creating tiledb array")
        // ensure the top folder exists
        File(dbPath).mkdirs()
        TiledbCoreHvcfUtils.createTileDBCoreArrays(dbPath)

        // Load the imputed hvcf file

        var vcfReader = VCFFileReader(File(imputedHvcf), false)
        val altHeadersImputed = TiledbCoreHvcfUtils.parseTiledbAltHeaders(vcfReader)
        println("Finished parsing imputedHvcf alt headers ")
        val variantData = TiledbCoreHvcfUtils.parseTiledbVariantData(vcfReader)
        println("Finished parsing imputedHvcf variant data")

        // Now load the data
        val arrayName = "${dbPath}/alt_header_array"
        TiledbCoreHvcfUtils.writeAltDataToTileDB(arrayName, altHeadersImputed)
        val variantArrayName = "${dbPath}/hvcf_variants_array"
        TiledbCoreHvcfUtils.writeVariantsDataToTileDB(variantArrayName, variantData)

        // Do again for LineA and LineB
        vcfReader = VCFFileReader(File(lineAhvcf), false)
        val altHeadersLineA = TiledbCoreHvcfUtils.parseTiledbAltHeaders(vcfReader)
        TiledbCoreHvcfUtils.writeAltDataToTileDB(arrayName, altHeadersLineA)
        val variantDataLineA = TiledbCoreHvcfUtils.parseTiledbVariantData(vcfReader)
        TiledbCoreHvcfUtils.writeVariantsDataToTileDB(variantArrayName, variantDataLineA)
        println("Finished writing LineA data")

        vcfReader = VCFFileReader(File(lineBhvcf), false)
        val altHeadersLineB = TiledbCoreHvcfUtils.parseTiledbAltHeaders(vcfReader)
        TiledbCoreHvcfUtils.writeAltDataToTileDB(arrayName, altHeadersLineB)
        val variantDataLineB = TiledbCoreHvcfUtils.parseTiledbVariantData(vcfReader)
        TiledbCoreHvcfUtils.writeVariantsDataToTileDB(variantArrayName, variantDataLineB)

        // Query for distinct sampleNames from the hvcf_variants_array
        println("QUery distinct sample names from the hvcf_variants_array")
        val names = queryDistinctSampleNames(variantArrayName)
        println("Distinct sample names: $names")
        assertTrue(names.size == 3)
        assertTrue(names.contains("LineA"))
        assertTrue(names.contains("LineB"))
        assertTrue(names.contains("TestLine2"))

        //Query for distinct refRanges from the hvcf_variants_array
        println("Query distinct refRanges from the hvcf_variants_array")
        val refRangeList = queryDistinctRefRanges(variantArrayName)
        println("\nnumber of distinct ref ranges: ${refRangeList.size}")
        println("Distinct refRanges: $refRangeList")
        assertEquals(38, refRangeList.size)

        File(dbPath).deleteRecursively()
    }

}