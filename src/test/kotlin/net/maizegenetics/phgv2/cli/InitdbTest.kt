package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.assertEquals

class InitdbTest {
    companion object {
        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"

        @JvmStatic
        @BeforeAll
        fun setup() {
            File(tempDir).mkdirs()
        }

        // Comment these out if you need to look at the logs files created by ProcessBuilder
        // commands.
        @JvmStatic
        @AfterAll
        fun teardown() {
            File(tempDir).deleteRecursively()
        }
    }
    //@Ignore
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
        val expectedFiles = listOf("${dbPath}gvcf_dataset", "${dbPath}hvcf_dataset", "${dbPath}temp/tiledb_gvcf_createURI_error.log", "${dbPath}temp/tiledb_gvcf_createURI_output.log", "${dbPath}temp/tiledb_hvcf_createURI_error.log", "${dbPath}temp/tiledb_hvcf_createURI_output.log")
        expectedFiles.forEach { assert(File( it).exists()) }

        // Delete the log files created by the ProcessBuilder commands.
        // Rerun the command .  The log files should NOT be created as
        // the code will realize the dataset already exists and does not try to recreate it.
        val expectedLogFiles = listOf("${dbPath}temp/tiledb_gvcf_createURI_output.log", "${dbPath}temp/tiledb_hvcf_createURI_error.log")
        expectedLogFiles.forEach { File(it).delete() }

        // Run this again to verify that the datasets are not overwritten.
        val result2 = initdb.createDataSets( dbPath)
        println("result2 = $result2")

        // The output and error log files are not recreated because the
        //  datasets exist,

        // verify the expectedLogFiles (ie the ProcessBuilder log files) do NOT exist
        expectedLogFiles.forEach { assert(!(File( it).exists())) }

        // verify the datasets still exist
        val expectedDatasetFiles = listOf("${dbPath}gvcf_dataset", "${dbPath}hvcf_dataset")
        expectedDatasetFiles.forEach { assert(File(it).exists()) }
    }

    @Test
    fun testAnchorGaps() {
        // testing that we can change the anchor gap size on both the gvcf and hvcf databases

        val dbPath = tempDir + "tiledb_datasets_anchor/"
        val initdb = Initdb()
        val gGap_truth = 10000
        val hGap_truth = 5000

        val result = initdb.test("--db-path $dbPath -g $gGap_truth -h $hGap_truth")

        assertEquals(result.statusCode, 0)

        // get stats to show that anchor-gap was set correctly for gvcf
        var builder = ProcessBuilder("conda","run","-n","phgv2-conda","tiledbvcf","stat","--uri", "$dbPath/gvcf_dataset/")
        var process = builder.start()
        var stats = process.inputStream.readAllBytes().decodeToString()
        process.waitFor()
        val gGap = stats.split("\n").filter{it.startsWith("- Anchor gap: ")}.map{it.replace("- Anchor gap: ", "").replace(",","")}[0].toInt()

        assertEquals(gGap_truth, gGap)

        // get stats to show that anchor-gap was set correctly for hvcf
        builder = ProcessBuilder("conda","run","-n","phgv2-conda","tiledbvcf","stat","--uri", "$dbPath/hvcf_dataset/")
        process = builder.start()
        stats = process.inputStream.readAllBytes().decodeToString()
        process.waitFor()
        val hGap = stats.split("\n").filter{it.startsWith("- Anchor gap: ")}.map{it.replace("- Anchor gap: ", "").replace(",","")}[0].toInt()

        assertEquals(hGap_truth, hGap)

    }
}
