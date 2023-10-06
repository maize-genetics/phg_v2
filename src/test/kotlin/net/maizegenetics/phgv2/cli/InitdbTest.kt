package net.maizegenetics.phgv2.cli

import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.Ignore

class InitdbTest {
    companion object {
        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"

        @JvmStatic
        @BeforeAll
        fun setup() {
            File(tempDir).mkdirs()
        }

        // Comment these out if you need to look at the logs files created by ProcessBuilder
        // comda commands.
        @JvmStatic
        @AfterAll
        fun teardown() {
            File(tempDir).deleteRecursively()
        }
    }
    @Ignore
    @Test
    fun testInitdb() {

        // When running inidb, debug is printed to the screen.
        // IF an error occurs, the user is instructed to check the log files,
        // If the datasets, or folders with the dataset names already exist,
        // nothing is done.  Datasets are not overwritten.

        val dbPath = tempDir + "tiledb_datasets/"
        val initdb = Initdb()
        val result = initdb.createDataSets( dbPath)
        println("result = $result")
        // To verify, we  check that the outputDir contains the expected files
        val expectedFiles = listOf("gvcf_dataset", "hvcf_dataset")
        expectedFiles.forEach { assert(File(dbPath + it).exists()) }

        // Run this again to verify that the datasets are not overwritten.
        val result2 = initdb.createDataSets( dbPath)
        println("result2 = $result2")
        // To verify, we check the output log for a message that the files
        // already existed.  HHH ... this goes to the myLogger file, not
        // a process builder output file.  Can I delete any existing output
        // files before call initDB?  Then if the are not recreated,  but
        // the datasets exist, it means they were there previously and were not
        // overwritten.
        val gvcfOutput = dbPath + "/tiledb_gvcf_createURI_output.log"

    }
}