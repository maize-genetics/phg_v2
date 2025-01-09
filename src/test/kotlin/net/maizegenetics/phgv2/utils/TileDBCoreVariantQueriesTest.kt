package net.maizegenetics.phgv2.utils

import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.cli.TestExtension
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File

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
        // This is currently not working correctly.

        // Make sure the top folder exists
        File(TiledbCoreHvcfUtilsTest.dbPath).mkdirs()

        // testing output from parseTiledbAltHeaders
        var vcfReader = VCFFileReader(File(TiledbCoreHvcfUtilsTest.lineAhvcf), false)

        val variantDataLineA = parseTiledbVariantData(vcfReader)
        println("Finished parsing lineAhvcf variant data")

        vcfReader = VCFFileReader(File(TiledbCoreHvcfUtilsTest.lineBhvcf), false)
        val variantDataLineB = parseTiledbVariantData(vcfReader)
        println("Finished parsing lineBhvcf variant data")

        // Create the tiledb array by calling the function createTileDBArray
        println("Creating tiledb array")
        //createTileDBArraySingleDimensionID(tiledbArrayName)
        createTileDBCoreArrays(TiledbCoreHvcfUtilsTest.dbPath)


        //Write to the variants array
        writeVariantsDataToTileDB(TiledbCoreHvcfUtilsTest.variantsArray, variantDataLineA)
        writeVariantsDataToTileDB(TiledbCoreHvcfUtilsTest.variantsArray, variantDataLineB)

        // Query the variants array
        val results = queryByRefRanges(TiledbCoreHvcfUtilsTest.variantsArray, listOf("1:1-1000","1:12001-16500","1:33001-34000","2:1-1000","2:5501-6500"))
        println("\nResults on variantsArray:")
        results.forEach { println(it) }


        // Manually verifying - need asserts here !!

        // delete the tiledbArray so next tests can recreate what they need
        File(TiledbCoreHvcfUtilsTest.dbPath).deleteRecursively()

    }

    @Test
    fun testVariantArray() {
        println("Creating tiledb array")
        // ensure the top folder exists
        File(TiledbCoreHvcfUtilsTest.dbPath).mkdirs()
        createTileDBCoreArrays(TiledbCoreHvcfUtilsTest.dbPath)

        // Load the imputed hvcf file
        // loadImputedHvcfToTileDB(imputedHvcf, dbPath)
        // testing output from parseTiledbAltHeaders
        var vcfReader = VCFFileReader(File(TiledbCoreHvcfUtilsTest.imputedHvcf), false)
        val altHeadersImputed = parseTiledbAltHeaders(vcfReader)
        println("Finished parsing imputedHvcf alt headers ")
        val variantData = parseTiledbVariantData(vcfReader)
        println("Finished parsing imputedHvcf variant data")

        // Now load the data
        val arrayName = "${TiledbCoreHvcfUtilsTest.dbPath}/alt_header_array"
        writeAltDataToTileDB(arrayName, altHeadersImputed)
        val variantArrayName = "${TiledbCoreHvcfUtilsTest.dbPath}/hvcf_variants_array"
        writeVariantsDataToTileDB(variantArrayName, variantData)

        // Do again for LineA and LineB
        vcfReader = VCFFileReader(File(TiledbCoreHvcfUtilsTest.lineAhvcf), false)
        val altHeadersLineA = parseTiledbAltHeaders(vcfReader)
        writeAltDataToTileDB(arrayName, altHeadersLineA)
        val variantDataLineA = parseTiledbVariantData(vcfReader)
        writeVariantsDataToTileDB(variantArrayName, variantDataLineA)
        println("Finished writing LineA data")

        vcfReader = VCFFileReader(File(TiledbCoreHvcfUtilsTest.lineBhvcf), false)
        val altHeadersLineB = parseTiledbAltHeaders(vcfReader)
        writeAltDataToTileDB(arrayName, altHeadersLineB)
        val variantDataLineB = parseTiledbVariantData(vcfReader)
        writeVariantsDataToTileDB(variantArrayName, variantDataLineB)

        // Query for distinct sampleNames from the hvcf_variants_array
        println("QUery distinct sample names from the hvcf_variants_array")
        val names = queryDistinctSampleNames(variantArrayName)
        println("Distinct sample names: $names")
        Assertions.assertTrue(names.size == 3)
        Assertions.assertTrue(names.contains("LineA"))
        Assertions.assertTrue(names.contains("LineB"))
        Assertions.assertTrue(names.contains("TestLine2"))

        //Query for distinct refRanges from the hvcf_variants_array
        println("Query distinct refRanges from the hvcf_variants_array")
        val refRangeList = queryDistinctRefRanges(variantArrayName)
        println("\nnumber of distinct ref ranges: ${refRangeList.size}")
        println("Distinct refRanges: $refRangeList")
        Assertions.assertEquals(38, refRangeList.size)


        // QUery for SampleName and ID based on a list of RefRanges
        // I can't get this to work.  It is fine when I try the query on the altHeaderArray
        //  but not on the variantsArray.

        //val refRanges = listOf("1:1-1000","1:12001-16500","1:33001-34000","2:1-1000","2:5501-6500")
        val refRanges = listOf("1:33001-34000","2:1-1000","2:5501-6500")
        //val results = queryByRefRanges(variantArrayName, refRanges)
        val results = queryByRefRanges(variantArrayName, listOf("1:1-1000","1:12001-16500","1:33001-34000","2:1-1000","2:5501-6500"))
        println("\nResults on variantsArray:")
        results.forEach { println(it) }

        // Query as we do the altHeaders
        val resultsByRefRange2 = queryIDsByRefRange(TiledbCoreHvcfUtilsTest.altHeaderArray, listOf("1:1-1000","1:12001-16500","1:33001-34000","2:1-1000","2:5501-6500"))

        println("\nResults with 5 refRanges:, altHeaderArray: ")
        results.forEach {println(it)}

        File(TiledbCoreHvcfUtilsTest.dbPath).deleteRecursively()
    }

    @Test
    fun testLoadImputedHvcf() {
        // TODO - this failes - no output!1 Just columns, no rows.
        // Load the imputed hvcf file into the tiledb core arrays
        // Create the tiledb array by calling the function createTileDBArray
        println("Creating tiledb array")
        // ensure the top folder exists
        File(TiledbCoreHvcfUtilsTest.dbPath).mkdirs()
        createTileDBCoreArrays(TiledbCoreHvcfUtilsTest.dbPath)

        // Load the imputed hvcf file
        // loadImputedHvcfToTileDB(imputedHvcf, dbPath)
        // testing output from parseTiledbAltHeaders
        var vcfReader = VCFFileReader(File(TiledbCoreHvcfUtilsTest.imputedHvcf), false)
        val altHeadersImputed = parseTiledbAltHeaders(vcfReader)
        println("Finished parsing imputedHvcf alt headers ")
        val variantData = parseTiledbVariantData(vcfReader)
        println("Finished parsing imputedHvcf variant data")

        // Now load the data
        val arrayName = "${TiledbCoreHvcfUtilsTest.dbPath}/alt_header_array"
        writeAltDataToTileDB(arrayName, altHeadersImputed)
        val variantArrayName = "${TiledbCoreHvcfUtilsTest.dbPath}/hvcf_variants_array"
        writeVariantsDataToTileDB(variantArrayName, variantData)

        println("QUery distinct sample names from the hvcf_variants_array")
        val names = queryDistinctSampleNames(variantArrayName)
        println("Distinct sample names: $names")

        // Query the data to verify it is loaded.
        // QUery only the variants data:
        val results = queryVariantsArray(variantArrayName, null, listOf("TestLine2"))

        // delete the tiledbArray so next tests can recreate what they need
        File(TiledbCoreHvcfUtilsTest.dbPath).deleteRecursively()
    }
}