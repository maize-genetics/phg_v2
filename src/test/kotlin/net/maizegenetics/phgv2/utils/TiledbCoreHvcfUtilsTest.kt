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
        val tiledbArrayName = "${TestExtension.tempDir}/alt_header_array"

        @BeforeAll
        @JvmStatic
        fun setup() {
            File(multiInputDir).mkdirs()
            File(outputHvcfDir).mkdirs()

            File(TestExtension.smallseqRefHvcfFile).copyTo(File(multiInputDir + File(TestExtension.smallseqRefHvcfFile).name))
            File(TestExtension.smallseqLineAHvcfFile).copyTo(File(multiInputDir + File(TestExtension.smallseqLineAHvcfFile).name))
            File(TestExtension.smallseqLineBHvcfFile).copyTo(File(multiInputDir + File(TestExtension.smallseqLineBHvcfFile).name))


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
    fun testParseHeadersWriteFile() {
        // testing output from parseTiledbAltHeaders
        val altHeaders = parseTiledbAltHeaders(TestExtension.smallseqLineAHvcfFile)
        println("Finished testing parseTiledbAltHeaders")

        // Create the tiledb array by calling the function createTileDBArray
        println("Creating tiledb array")
        createTileDBArray(tiledbArrayName)

        // Write the altHeaders to the tiledb array
        println("Writing altHeaders to tiledb array")
        writeAltDataToTileDB(tiledbArrayName, altHeaders)

        // extract data
        println("Extracting data from tiledb array")
        val idsToQuery = listOf("0eb9029f3896313aebc69c8489923141")

        // Read data from TileDB
        val results = readDataFromTileDB(tiledbArrayName, idsToQuery)

        // Print results
        println("Results:")
        results.forEach { println(it) }

    }

}