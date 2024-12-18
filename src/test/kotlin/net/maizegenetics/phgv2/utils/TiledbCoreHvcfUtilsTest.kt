package net.maizegenetics.phgv2.utils

import biokotlin.util.bufferedReader
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
        val lineAhvcf = "data/test/tiledbCoreHvcf/LineA.h.vcf"
        val lineBhvcf = "data/test/tiledbCoreHvcf/LineB.h.vcf"
        val dbPath = TestExtension.testTileDBURI
        val altHeaderArray = dbPath + "/alt_header_array"
        val variantsArray = dbPath + "/hvcf_variants_array"
        val edArray = dbPath + "/ed_hvcf_array"

        @BeforeAll
        @JvmStatic
        fun setup() {
            File(multiInputDir).mkdirs()
            File(outputHvcfDir).mkdirs()
            File(dbPath).mkdirs()

            File(lineAhvcf).copyTo(File(multiInputDir + File(lineAhvcf).name))
            File(lineBhvcf).copyTo(File(multiInputDir + File(lineBhvcf).name))
        }

        @AfterAll
        @JvmStatic
        fun teardown() {
            File(outputHvcfDir).deleteRecursively()
            File(multiInputDir).deleteRecursively()
            // delete the tempDir
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

//    @Test
//    fun testEdTiledbArray() {
//        createTileDBEdArray(dbPath)
//        val array = Array(Context(), edArray, QueryType.TILEDB_READ)
//        assertTrue(array.schema != null)
//        array.close()
//
//    }
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
    fun testParseALtHeadersTiledb() {
        // testing output from parseTiledbAltHeaders
        val altHeaders = parseTiledbAltHeaders(lineAhvcf)
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
        val altHeadersLineA = parseTiledbAltHeaders(lineAhvcf)
        println("Finished parsing lineAhvcf")

        val altHeadersLineB = parseTiledbAltHeaders(lineBhvcf)
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
        // These are the first 2 IDs in the LineA.h.vcf file followed by the first 2 IDs in the LineB.h.vcf file
        //val rangesToQuery = listOf("=1:1-1000","=2:1-1000")
        val rangesToQuery = listOf("1:1-1000","2:1-1000")
        //val idsToQuery = listOf("0eb9029f3896313aebc69c8489923141","12f0cec9102e84a161866e37072443b7","00f297caa4a0fa5a6f8e76d388393fa7","05efe15d97db33185b64821791b01b0f")

        // Read data from TileDB.  This gets the SampleName and Regions for the
        // given ID
        val results = querySampleNamesAndIDsByRefRange(altHeaderArray, rangesToQuery)

        // Print results
        println("Results:")
        results.forEach { println(it) }

        // delete the tiledbArray so next tests can recreate what they need
        File(dbPath).deleteRecursively()
    }
    @Test
    fun testQueryByRefRange() {

        // testing output from parseTiledbAltHeaders
        val altHeadersLineA = parseTiledbAltHeaders(lineAhvcf)
        println("Finished parsing LineA alt headers ")
        val altHeadersLineB = parseTiledbAltHeaders(lineBhvcf)
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

        // delete the tiledbArray so next tests can recreate what they need
        File(dbPath).deleteRecursively()

    }

}