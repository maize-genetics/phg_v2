package net.maizegenetics.phgv2.utils

import net.maizegenetics.phgv2.cli.TestExtension
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
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
        val tiledbArrayName = "${TestExtension.tempDir}/alt_header_array"

        @BeforeAll
        @JvmStatic
        fun setup() {
            File(multiInputDir).mkdirs()
            File(outputHvcfDir).mkdirs()

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

    @Test
    fun testWritingToTiledb() {
        // testing output from parseTiledbAltHeaders
        val altHeadersLineA = parseTiledbAltHeaders(lineAhvcf)
        println("Finished parsing smallSeqLineAHvcfFile")

        val altHeadersLineB = parseTiledbAltHeaders(lineBhvcf)
        println("Finished parsing smallseqLineBHvcfFile")

        // Create the tiledb array by calling the function createTileDBArray
        println("Creating tiledb array")
        createTileDBArray2Dimension(tiledbArrayName)

        // Write the altHeaders to the tiledb array
        println("Writing LineA altHeaders to tiledb array")
        writeAltDataToTileDB(tiledbArrayName, altHeadersLineA)
        println("Writing LineB altHeaders to tiledb array")
        writeAltDataToTileDB(tiledbArrayName, altHeadersLineB)
    }
    @Test
    fun testParseHeadersWriteFile() {
        // testing output from parseTiledbAltHeaders
        val altHeadersLineA = parseTiledbAltHeaders(lineAhvcf)
        println("Finished parsing lineAhvcf")

        val altHeadersLineB = parseTiledbAltHeaders(lineBhvcf)
        println("Finished parsing lineBhvcf")

        // Create the tiledb array by calling the function createTileDBArray
        println("Creating tiledb array")
        createTileDBArray2Dimension(tiledbArrayName)

        // Write the altHeaders to the tiledb array
        println("Writing LineA altHeaders to tiledb array")
        writeAltDataToTileDB(tiledbArrayName, altHeadersLineA)
        println("Writing LineB altHeaders to tiledb array")
        writeAltDataToTileDB(tiledbArrayName, altHeadersLineB)

        // extract data
        println("TestCase: Extracting data from tiledb array")
        val idsToQuery = listOf("0eb9029f3896313aebc69c8489923141")

        // Read data from TileDB
        val results = readSampleNameRegionsForID(tiledbArrayName, idsToQuery)

        // Print results
        println("Results:")
        results.forEach { println(it) }

        // delete the tiledbArray so next tests can recreate what they need
        File(tiledbArrayName).deleteRecursively()

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
        createTileDBArray2Dimension(tiledbArrayName)

        // Write the altHeaders to the tiledb array
        println("Writing altHeaders to tiledb array")
        writeAltDataToTileDB(tiledbArrayName, altHeadersLineA)
        writeAltDataToTileDB(tiledbArrayName, altHeadersLineB)
        // Try to query by refRange
        val resultsByRefRange = queryIDsByRefRange(tiledbArrayName, listOf("1:1-1000"))
        println("Results by refRange: $resultsByRefRange")

        // delete the tiledbArray so next tests can recreate what they need
        File(tiledbArrayName).deleteRecursively()

    }

}