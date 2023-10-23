package net.maizegenetics.phgv2.utils

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.AgcCompress
import net.maizegenetics.phgv2.cli.CreateRefVcf
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.cli.TestExtension.Companion.testOutputFastaDir
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import java.nio.file.Files
import java.nio.file.Paths
import kotlin.test.assertEquals
import kotlin.test.assertTrue
import org.junit.jupiter.api.assertThrows

class SeqUtilsTest {
    companion object {

        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"
        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.tempDir).mkdirs()

            // create the agc compressed file from which we'll pull data
            val fastaInputDir = "data/test/smallseq"
            val fastaOutputDir = TestExtension.testOutputFastaDir
            val dbPath = TestExtension.testTileDBURI

            // Why isn't this handled in TestExtension - I have added the command "File(testTileDBURI).mkdirs()"
            // to TestExtension
            Files.createDirectories(Paths.get(dbPath))

            // copy files with extension .fa from data/test/smallseq to fastaOutputDir
            val fastaFiles = File(fastaInputDir).listFiles { file -> file.extension == "fa" }
            fastaFiles.forEach { file -> file.copyTo(File(fastaOutputDir, file.name)) }


            // get the full path fasta file names  from the fastaInput, write to fileList
            val fileList = mutableListOf<String>()
            File(fastaOutputDir).walkTopDown().filter{it.name.endsWith(".fa")}.forEach {
                fileList.add(it.toString())
            }

            // The list of fasta files to load should not contain the reference.  The reference is listed separately.
            fileList.removeIf { it.contains("Ref") }

            // write the full path fasta file names  from the fastaOutputDir to a single file, one per line, in tempDir
            // This file will be used as input to the agc-compress command
            val fastaCreateFileNamesFile = File(dbPath, "fastaCreateFileNames.txt")
            fastaCreateFileNamesFile.writeText(fileList.joinToString("\n"))

            val refFasta = File(fastaOutputDir, "Ref.fa").toString()

            val agcCompress = AgcCompress()
            // Create the initial compressed file
            println("Calling agcCompress for CREATE")
            var result = agcCompress.test("--fasta-list ${fastaCreateFileNamesFile} --db-path ${dbPath} --ref-fasta ${refFasta}")

        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testBuildAgcCommandFromList() {
        // test command with single range on the list

        val dbPath = TestExtension.testTileDBURI
        var rangeList = mutableListOf<String>()
        val range1 = "1@LineA:20-40"
        rangeList.add(range1)
        var command = buildAgcCommandFromList(dbPath,rangeList)

        //The command should look like:
        val expectedCommand = arrayOf("conda","run","-n","phgv2-conda","agc","getctg","${dbPath}/assemblies.agc",range1)

        // Verify that command contains all elements from expected command and that the order is the same
        kotlin.test.assertEquals(expectedCommand.size, command.size)
        for (i in 0..expectedCommand.size-1) {
            kotlin.test.assertEquals(expectedCommand[i], command[i])
        }

        // Add multiple ranges to the list.  Verify the command contains all range requests
        val range2 = "1@LineA:60-80"
        val range3 = "1@LineA:100-120"

        rangeList.add(range2)
        rangeList.add(range3)
        command = buildAgcCommandFromList(dbPath,rangeList)
        val expectedCommand2 = arrayOf("conda","run","-n","phgv2-conda","agc","getctg","${dbPath}/assemblies.agc",range1,range2,range3)
        kotlin.test.assertEquals(expectedCommand2.size, command.size)
        for (i in 0..expectedCommand2.size-1) {
            kotlin.test.assertEquals(expectedCommand2[i], command[i])
        }

    }

    @Test
    fun testRetrieveAgcData() {
        // This tests the function retrieveAgcData.  retrieveAgcData() calls buildAgcCommandFromList() and
        // and then queryAgc() to get the data from the agc compressed file.  The data is returned as a
        // Map<String,NucSeq> where "String" is idline from the AGC created fasta, and NucSeq is the sequence.

        val dbPath = "${TestExtension.testTileDBURI}" // just the path, code appends "assemblies.agc"
        var rangeList = mutableListOf<String>()
        val range1 = "1@LineA:0-19" // AGC queries are 0-based !!
        rangeList.add(range1)
        var agcResult = retrieveAgcContigs(dbPath,rangeList)

        // this has only a single query, so the result should have 1 entry: 1 key and 1 sequence
        assertEquals(1, agcResult.size)
        assertEquals(1, agcResult.keys.size)
        // verify the key as "1:0-19".  The genome is not included in AGC's id line
        assertEquals("1:0-19", agcResult.keys.first())
        // verify the sequence - this is copied from LineA:1-20 in the fasta file
        assertEquals("GCGCGGGGACCGAGAAACCC", agcResult.values.first().toString())

        // add additional ranges to the list
        val range2 = "1@LineA:20-39"
        val range3 = "1@LineA:40-59"

        rangeList.add(range2)
        rangeList.add(range3)
        agcResult = retrieveAgcContigs(dbPath,rangeList)
        // this has 3 queries, so the result should have 3 entries: 3 keys and 3 sequences
        assertEquals(3, agcResult.size)
        assertEquals(3, agcResult.keys.size)

        // This is a map, the keys will not be in any specific order, so verify
        // that each key exists in the list of expected keys
        val keyList = mutableListOf<String>("1:0-19","1:20-39","1:40-59")
        // verify the map returned contains all keys from the keyList
        assertTrue(agcResult.keys.containsAll(keyList))

        // verify the sequences: these are copied from LineA:1-20, LineA:21-40, and LineA:41-60 in the fasta file
        // used in the setup() function
        assertEquals("GCGCGGGGACCGAGAAACCC", agcResult["1:0-19"].toString())
        assertEquals("GGCGGGGCAGGACGAACCGG", agcResult["1:20-39"].toString())
        assertEquals("GCGAGAACGGACAACCCCCC", agcResult["1:40-59"].toString())

    }

    @Test
    fun testBadAGCParameter() {
        // This function tests what happens when a bad parameter is passed to the agc command
        // The agc command should return an error message, and the function should throw an exception

        // This db path does not exist
        var dbPath = "/my/bad/dbPath"
        var rangeList = mutableListOf<String>()
        val range1 = "1@LineA:0-19" // AGC queries are 0-based !!
        rangeList.add(range1)
        assertThrows<IllegalStateException> {
            //Check that an error is thrown when the dbPath is not a directory that contains the assemblies.agc file
            retrieveAgcContigs(dbPath,rangeList)
        }

        // Test with a valid directory name, but the assemblies.agc file does not exist

        dbPath = testOutputFastaDir
        assertThrows<IllegalStateException> {
            //Check that an error is thrown when the dbPath is not a directory that contains the assemblies.agc file
            retrieveAgcContigs(dbPath,rangeList)
        }

    }
}