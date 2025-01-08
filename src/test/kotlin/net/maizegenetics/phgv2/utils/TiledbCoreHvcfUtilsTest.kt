package net.maizegenetics.phgv2.utils

import biokotlin.util.bufferedReader
import htsjdk.variant.vcf.VCFFileReader
import io.tiledb.java.api.Array
import io.tiledb.java.api.QueryType
import io.tiledb.java.api.*
import net.maizegenetics.phgv2.cli.TestExtension
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File

@ExtendWith(TestExtension::class)
class TiledbCoreHvcfUtilsTest {

    private val myLogger = LogManager.getLogger(TiledbCoreHvcfUtilsTest::class.java)

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
    fun testTiledbArrayExists() {
        // This is mostly to  verify the commands that test for
        // tiledb array existence
        // make sure tiledb main folder exists
        File(dbPath).mkdirs()
        println("Creating tiledb array")
        createTileDBCoreArrays(dbPath)
        val array = Array(Context(), altHeaderArray, QueryType.TILEDB_READ)
        assertTrue(array.schema != null)
        array.close()

        // verify  the variants array exists
        val variantsArray = Array(Context(), variantsArray, QueryType.TILEDB_READ)
        assertTrue(variantsArray.schema != null)
        variantsArray.close()

        // Verify an exception is thrown if name exists but is not a tiledb array
        assertThrows<TileDBError>{
            Array(Context(), lineAhvcf, QueryType.TILEDB_READ)
        }

        // Verify exception is thrown if the array does not exist
        assertThrows<TileDBError>{
            Array(Context(), "badArray", QueryType.TILEDB_READ)
        }

        // delete the tiledbArray so next tests can recreate what they need
        File(dbPath).deleteRecursively()

    }

    @Test
    fun testParseAltHeadersTiledb() {
        // testing output from parseTiledbAltHeaders
        val vcfReader = VCFFileReader(File(lineAhvcf), false)
        val altHeaders = parseTiledbAltHeaders(vcfReader)
        println("Finished parsing smallSeqLineAHvcfFile")

        // Read the altheader data into a list to compare to what
        // is stored from the parseTiledbAltHeaders function
        // Use buffered reader to read the lines from lineAhvcf, but keep only the lines
        // that start with "##ALT"
        val altHeaderListFromHvcf = bufferedReader(lineAhvcf).readLines().filter { it.startsWith("##ALT") }
        assertTrue(altHeaderListFromHvcf.size == 38)
        assertEquals(altHeaders.size, altHeaderListFromHvcf.size)

        // Extract data
        // Note the ALT headers are stored in the h.vcf file ordered lexicographically by checksum (Hapid)
        val ids = altHeaders.map { it["ID"].orEmpty() }
        val sampleNames = altHeaders.map { it["SampleName"].orEmpty() }
        val regions = altHeaders.map { it["Regions"].orEmpty() }
        val refChecksums = altHeaders.map { it["RefChecksum"].orEmpty() }
        val refRanges = altHeaders.map { it["RefRange"].orEmpty() }

        // verify the number of entries from parseTiledbAltHeaders matches the number of ALT header lines in the VCF file
        assertEquals(altHeaders.size, altHeaderListFromHvcf.size)

        // THe first ALT header line has 1 region in it: 2:49501-50300, verify it is stored correctly
        assertEquals("2:49501-50300", regions[0])
        // The second ALT header line has 2 regions in it: Regions="1:1-950,1:952-1000"
        // Verify that the data is extracted correctly
        assertEquals("1:1-950,1:952-1000", regions[1])
        // Verify the sample name is LineA
        assertEquals("LineA", sampleNames[1])
        // Verify the hapID of the second entry is 12f0cec9102e84a161866e37072443b7
        assertEquals("12f0cec9102e84a161866e37072443b7", ids[1])
        // THe last line in the file is this: verify the values in the parsed header data
        //##ALT=<ID=f50fe6d6b3d9a9d305889db977969916,Description="haplotype data for line: LineA",
        // Source="/Users/lcj34/temp/phgv2Tests/tempDir/testTileDBURI//assemblies.agc",SampleName=LineA,
        // Regions=1:11001-12000,Checksum=f50fe6d6b3d9a9d305889db977969916,
        // RefChecksum=2b4590f722ef9229c15d29e0b4e51a0e,RefRange=1:11001-12000>

        assertEquals("f50fe6d6b3d9a9d305889db977969916", ids.last())
        assertEquals("LineA", sampleNames.last())
        assertEquals("1:11001-12000", regions.last())
        assertEquals("2b4590f722ef9229c15d29e0b4e51a0e", refChecksums.last())
        assertEquals("1:11001-12000", refRanges.last())

        println("Finished extracting data")
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

            val rangesToQuery = listOf("1:1-1000","1:6501-11000","2:1-1000","2:5501-6500")

        // Try with dynamic buffers
        val streamingResults = queryWithStreaming_sampleNameIdByRefRange(altHeaderArray, rangesToQuery)
        println("\nStreaming buffer results:")
        streamingResults.forEach { println(it) }

            // delete the tiledbArray so next tests can recreate what they need
            File(dbPath).deleteRecursively()
    }

    @Test
    fun testQueryByRefRange() {

        // testing output from parseTiledbAltHeaders
        var vcfReader = VCFFileReader(File(lineAhvcf), false)
        val altHeadersLineA = parseTiledbAltHeaders(vcfReader)
        println("Finished parsing LineA alt headers ")
        vcfReader = VCFFileReader(File(lineBhvcf), false)
        val altHeadersLineB = parseTiledbAltHeaders(vcfReader)
        println("Finished parsing lineBhvcf alt headers")

        // Create the tiledb array by calling the function createTileDBArray
        println("Creating tiledb array")
        //createTileDBArraySingleDimensionID(tiledbArrayName)
        createTileDBCoreArrays(dbPath)

        // Write the altHeaders to the tiledb array
        println("Writing altHeaders to tiledb array")
        writeAltDataToTileDB(altHeaderArray, altHeadersLineA)
        writeAltDataToTileDB(altHeaderArray, altHeadersLineB)
        // Try to query by refRange
        val resultsByRefRange = queryIDsByRefRange(altHeaderArray, listOf("1:1-1000"))
        println("Results with 1 refRange: $resultsByRefRange")
        val resultsByRefRange2 = queryIDsByRefRange(altHeaderArray, listOf("1:1-1000","2:12001-16500"))
        println("Results with 2 refRanges: $resultsByRefRange2")

        val resultsByRefRange5Ranges = queryIDsByRefRange(altHeaderArray, listOf("1:1-1000","1:12001-16500","1:33001-34000","2:1-1000","2:5501-6500"))
        println("Results with 5 refRanges: $resultsByRefRange5Ranges")

        // Need some asserts here !!

        // delete the tiledbArray so next tests can recreate what they need
        File(dbPath).deleteRecursively()

    }

    @Test
    fun testQueryByRefRange_2() {

        // Trying to see why I can filter on specific refRanges and get it correct for the alt header array,
        // but cannot get it correct for the variants array  Added the variants array processing here,
        // which is a copy of the test above.  Wanted to be sure I didn't mess up the test above which is
        // the reason for the duplicate test case.

        // testing output from parseTiledbAltHeaders
        var vcfReader = VCFFileReader(File(lineAhvcf), false)
        val altHeadersLineA = parseTiledbAltHeaders(vcfReader)
        println("Finished parsing LineA alt headers ")
        val variantDataLineA = parseTiledbVariantData(vcfReader)
        println("Finished parsing lineAhvcf variant data")

        vcfReader = VCFFileReader(File(lineBhvcf), false)
        val altHeadersLineB = parseTiledbAltHeaders(vcfReader)
        println("Finished parsing lineBhvcf alt headers")
        val variantDataLineB = parseTiledbVariantData(vcfReader)
        println("Finished parsing lineBhvcf variant data")

        // Create the tiledb array by calling the function createTileDBArray
        println("Creating tiledb array")
        //createTileDBArraySingleDimensionID(tiledbArrayName)
        createTileDBCoreArrays(dbPath)

        // Write the altHeaders to the tiledb array
        println("Writing altHeaders to tiledb array")
        writeAltDataToTileDB(altHeaderArray, altHeadersLineA)
        writeAltDataToTileDB(altHeaderArray, altHeadersLineB)

        //Write to the variants array
        writeVariantsDataToTileDB(variantsArray, variantDataLineA)
        writeVariantsDataToTileDB(variantsArray, variantDataLineB)

        // Try to query by refRange
        val resultsByRefRange = queryIDsByRefRange(altHeaderArray, listOf("1:1-1000"))
        println("Results with 1 refRange: $resultsByRefRange")
        val resultsByRefRange2 = queryIDsByRefRange(altHeaderArray, listOf("1:1-1000","2:12001-16500"))
        println("Results with 2 refRanges: $resultsByRefRange2")

        val resultsByRefRange5Ranges = queryIDsByRefRange(altHeaderArray, listOf("1:1-1000","1:12001-16500","1:33001-34000","2:1-1000","2:5501-6500"))
        println("Results with 5 refRanges: $resultsByRefRange5Ranges")

        // Need some asserts here !!

        // Query the variants array
        val results = queryByRefRanges(variantsArray, listOf("1:1-1000","1:12001-16500","1:33001-34000","2:1-1000","2:5501-6500"))
        println("\nResults on variantsArray:")
        results.forEach { println(it) }


        // delete the tiledbArray so next tests can recreate what they need
        File(dbPath).deleteRecursively()

    }

    @Test
    fun testVariantArray() {
        println("Creating tiledb array")
        // ensure the top folder exists
        File(dbPath).mkdirs()
        createTileDBCoreArrays(dbPath)

        // Load the imputed hvcf file
        // loadImputedHvcfToTileDB(imputedHvcf, dbPath)
        // testing output from parseTiledbAltHeaders
        var vcfReader = VCFFileReader(File(imputedHvcf), false)
        val altHeadersImputed = parseTiledbAltHeaders(vcfReader)
        println("Finished parsing imputedHvcf alt headers ")
        val variantData = parseTiledbVariantData(vcfReader)
        println("Finished parsing imputedHvcf variant data")

        // Now load the data
        val arrayName = "${dbPath}/alt_header_array"
        writeAltDataToTileDB(arrayName, altHeadersImputed)
        val variantArrayName = "${dbPath}/hvcf_variants_array"
        writeVariantsDataToTileDB(variantArrayName, variantData)

        // Do again for LineA and LineB
        vcfReader = VCFFileReader(File(lineAhvcf), false)
        val altHeadersLineA = parseTiledbAltHeaders(vcfReader)
        writeAltDataToTileDB(arrayName, altHeadersLineA)
        val variantDataLineA = parseTiledbVariantData(vcfReader)
        writeVariantsDataToTileDB(variantArrayName, variantDataLineA)
        println("Finished writing LineA data")

        vcfReader = VCFFileReader(File(lineBhvcf), false)
        val altHeadersLineB = parseTiledbAltHeaders(vcfReader)
        writeAltDataToTileDB(arrayName, altHeadersLineB)
        val variantDataLineB = parseTiledbVariantData(vcfReader)
        writeVariantsDataToTileDB(variantArrayName, variantDataLineB)

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
        val resultsByRefRange2 = queryIDsByRefRange(altHeaderArray, listOf("1:1-1000","1:12001-16500","1:33001-34000","2:1-1000","2:5501-6500"))

        println("\nResults with 5 refRanges:, altHeaderArray: ")
        results.forEach {println(it)}
    }

    @Test
    fun testLoadImputedHvcf() {
        // TODO - this failes - no output!1 Just columns, no rows.
        // Load the imputed hvcf file into the tiledb core arrays
        // Create the tiledb array by calling the function createTileDBArray
        println("Creating tiledb array")
        // ensure the top folder exists
        File(dbPath).mkdirs()
        createTileDBCoreArrays(dbPath)

        // Load the imputed hvcf file
       // loadImputedHvcfToTileDB(imputedHvcf, dbPath)
        // testing output from parseTiledbAltHeaders
        var vcfReader = VCFFileReader(File(imputedHvcf), false)
        val altHeadersImputed = parseTiledbAltHeaders(vcfReader)
        println("Finished parsing imputedHvcf alt headers ")
        val variantData = parseTiledbVariantData(vcfReader)
        println("Finished parsing imputedHvcf variant data")

        // Now load the data
        val arrayName = "${dbPath}/alt_header_array"
        writeAltDataToTileDB(arrayName, altHeadersImputed)
        val variantArrayName = "${dbPath}/hvcf_variants_array"
        writeVariantsDataToTileDB(variantArrayName, variantData)

        println("QUery distinct sample names from the hvcf_variants_array")
        val names = queryDistinctSampleNames(variantArrayName)
        println("Distinct sample names: $names")

        // Query the data to verify it is loaded.
        // QUery only the variants data:
        val results = queryVariantsArray(variantArrayName, null, listOf("TestLine2"))

        // delete the tiledbArray so next tests can recreate what they need
        File(dbPath).deleteRecursively()
    }

}