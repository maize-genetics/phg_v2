package net.maizegenetics.phgv2.utils

import biokotlin.genome.fastaToNucSeq
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.AgcCompress
import net.maizegenetics.phgv2.cli.CreateMafVcf
import net.maizegenetics.phgv2.cli.Initdb
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.cli.TestExtension.Companion.testOutputFastaDir
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import java.io.File
import java.nio.file.Files
import java.nio.file.Paths
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue
import org.junit.jupiter.api.assertThrows
import org.junit.jupiter.api.extension.ExtendWith

@ExtendWith(TestExtension::class)
class SeqUtilsTest {

    companion object {
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
            File(fastaOutputDir).walkTopDown().filter { it.name.endsWith(".fa") }.forEach {
                fileList.add(it.toString())
            }

            // The list of fasta files to load should not contain the reference.  The reference is listed separately.
            fileList.removeIf { it.contains("Ref") }

            // write the full path fasta file names  from the fastaOutputDir to a single file, one per line, in tempDir
            // This file will be used as input to the agc-compress command
            val fastaCreateFileNamesFile = File(dbPath, "fastaCreateFileNames.txt")
            fastaCreateFileNamesFile.writeText(fileList.joinToString("\n"))

            val refFasta = File(fastaOutputDir, "Ref.fa").toString()

            Initdb().createDataSets(TestExtension.testTileDBURI)
            // contain fastas with the sampleName, and once when it does contain the sampleName
            val agcCompress = AgcCompress()
            // Create the initial compressed file
            println("Calling agcCompress for CREATE")
            var result =
                agcCompress.test("--fasta-list ${fastaCreateFileNamesFile} --db-path ${dbPath} --reference-file ${refFasta}")
            println("result output: ${result.output}")
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(TestExtension.tempDir).deleteRecursively()
        }
    }


    @Test
    fun testBuildAgcCommandFromList() {
        // test getctg command with single range on the list

        val dbPath = TestExtension.testTileDBURI
        var rangeList = mutableListOf<String>()
        val range1 = "1@LineA:20-40"
        rangeList.add(range1)
        var command = buildAgcCommandFromList(dbPath, "getctg",rangeList)

        //The command should look like:
        val expectedCommand =
            arrayOf("conda", "run", "-n", "phgv2-conda", "agc", "getctg", "${dbPath}/assemblies.agc", range1)

        // Verify that command contains all elements from expected command and that the order is the same
        kotlin.test.assertEquals(expectedCommand.size, command.size)
        for (i in 0..expectedCommand.size - 1) {
            kotlin.test.assertEquals(expectedCommand[i], command[i])
        }

        // Add multiple ranges to the list.  Verify the getctg command contains all range requests
        val range2 = "1@LineA:60-80"
        val range3 = "1@LineA:100-120"

        rangeList.add(range2)
        rangeList.add(range3)
        command = buildAgcCommandFromList(dbPath, "getctg",rangeList)
        val expectedCommand2 = arrayOf(
            "conda",
            "run",
            "-n",
            "phgv2-conda",
            "agc",
            "getctg",
            "${dbPath}/assemblies.agc",
            range1,
            range2,
            range3
        )
        assertEquals(expectedCommand2.size, command.size)
        for (i in 0..expectedCommand2.size - 1) {
            kotlin.test.assertEquals(expectedCommand2[i], command[i])
        }

    }

    @Test
    fun testRetrieveAgcContigsNoSampleName() {
        //This test is to verify queryAgc() throws an exception when an idLine is missing "sampleName="
        // The first idLine in the file has "sampleName", so AgcCompress() will take the file.
        // But we will have problems processing queries related to the contig with the idLine that does not contain "sampleName="
        // I don't expect this to happen in practice, but it is a possibility so adding this test to show an exception
        // is thrown.

        val dbPath = TestExtension.tempDir
        val refFasta = "data/test/smallseq/Ref.fa"
        val fastaCreateFileNamesFile = File(dbPath, "fastaBadNames.txt")
        fastaCreateFileNamesFile.writeText("data/test/agcTestBad/LineD_someSNMissing.fa\n")

        Initdb().createDataSets(TestExtension.tempDir)
        val agcCompress = AgcCompress()
        // Create the initial compressed file
        val agcCompressResult = agcCompress.test("--fasta-list ${fastaCreateFileNamesFile} --db-path ${dbPath} --reference-file ${refFasta}")

        // verify agcCompressResult
        assertEquals(0, agcCompressResult.statusCode)

        val rangeList = mutableListOf<String>()
        // contig 2 idline does not contain "sampleName="
        val range1 = "2@LineD_someSNMissing:0-19" // AGC queries are 0-based !!
        rangeList.add(range1)

        assertThrows<IllegalStateException> {
            //Check that an exception is thrown when the idline does not contain "sampleName="
            var agcResult = retrieveAgcContigs(dbPath, rangeList)
        }

    }

    @Test
    fun testRetrieveAgcContigs() {
        // This tests the function retrieveAgcContigs.  retrieveAgcContigs() calls buildAgcCommandFromList() and
        // and then queryAgc() to get the data from the agc compressed file.  The data is returned as a
        // Map<String,NucSeq> where "String" is idline from the AGC created fasta, and NucSeq is the sequence.

        val dbPath = "${TestExtension.testTileDBURI}" // just the path, code appends "assemblies.agc"
        var rangeList = mutableListOf<String>()
        val range1 = "1@LineA:0-19" // AGC queries are 0-based !!
        rangeList.add(range1)
        var agcResult = retrieveAgcContigs(dbPath, rangeList)

        // this has only a single query, so the result should have 1 entry: 1 key and 1 sequence
        assertEquals(1, agcResult.size)
        assertEquals(1, agcResult.keys.size)
        // verify the key as "1:0-19".  The genome is not included in AGC's id line
        //assertEquals("1:0-19", agcResult.keys.first())
        assertEquals(Pair("LineA","1:0-19"), agcResult.keys.first())
        // verify the sequence - this is copied from LineA:1-20 in the fasta file
        assertEquals("GCGCGGGGACCGAGAAACCC", agcResult.values.first().toString())

        // add additional ranges to the list
        val range2 = "1@LineA:20-39"
        val range3 = "1@LineA:40-59"

        rangeList.add(range2)
        rangeList.add(range3)
        agcResult = retrieveAgcContigs(dbPath, rangeList)
        // this has 3 queries, so the result should have 3 entries: 3 keys and 3 sequences
        assertEquals(3, agcResult.size)
        assertEquals(3, agcResult.keys.size)

        // This is a map, the keys will not be in any specific order, so verify
        // that each key exists in the list of expected keys
        //val keyList = mutableListOf<String>("1:0-19", "1:20-39", "1:40-59")
        val keyList = mutableListOf<Pair<String,String>>(Pair("LineA","1:0-19"), Pair("LineA","1:20-39"), Pair("LineA","1:40-59"))
        // verify the map returned contains all keys from the keyList
        assertTrue(agcResult.keys.containsAll(keyList))

        // verify the sequences: these are copied from LineA:1-20, LineA:21-40, and LineA:41-60 in the fasta file
        // used in the setup() function
        assertEquals("GCGCGGGGACCGAGAAACCC", agcResult[Pair("LineA","1:0-19")].toString())
        assertEquals("GGCGGGGCAGGACGAACCGG", agcResult[Pair("LineA","1:20-39")].toString())
        assertEquals("GCGAGAACGGACAACCCCCC", agcResult[Pair("LineA","1:40-59")].toString())

    }

    @Test
    fun testPullFullContigs() {
        val dbPath = "${TestExtension.testTileDBURI}" // just the path, code appends "assemblies.agc"
        var rangeList = mutableListOf<String>()
        val chr1 = "1@LineA" // AGC queries are 0-based !!
        val chr2 = "2@LineA"
        rangeList.add(chr1)
        rangeList.add(chr2)
        var agcResult = retrieveAgcContigs(dbPath, rangeList)

        assertEquals(2, agcResult.size)
        assertEquals(2, agcResult.keys.size)

        // Need to read the sequences from the fasta file to compare to the agcResult
        // The sequence lives in fasta file in the testOutputFastaDir directory in LineA.fa
        // The fasta file is copied from the testInputFastaDir directory in the setup() function
        val fastaFile = File("${testOutputFastaDir}/LineA.fa").toString()
        val lineA = fastaToNucSeq(fastaFile)

        // Verify that agcResult contains the same sequences as the fasta file for the 2 chromosomes
        // Might have to deal with newlines in the sequences stored in LineA.fa
        assertEquals(lineA["1"]!!.toString(), agcResult[Pair("LineA","1")]!!.toString())
        assertEquals(lineA["2"]!!.toString(), agcResult[Pair("LineA","2")]!!.toString())
    }

    @Test
    fun testBadContigName() {
        val dbPath = "${TestExtension.testTileDBURI}" // just the path, code appends "assemblies.agc"
        var rangeList = mutableListOf<String>()
        val chr1 = "1@LineZ" // LineZ does not exist in the AGC compressed fasta
        rangeList.add(chr1)

        assertThrows<IllegalArgumentException> {
            //Check that an error is thrown when a bad genome is passed.
            retrieveAgcContigs(dbPath, rangeList)
        }
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
            retrieveAgcContigs(dbPath, rangeList)
        }

        // Test with a valid directory name, but the assemblies.agc file does not exist

        dbPath = testOutputFastaDir
        assertThrows<IllegalStateException> {
            //Check that an error is thrown when the dbPath is not a directory that contains the assemblies.agc file
            retrieveAgcContigs(dbPath, rangeList)
        }
    }

    @Test
    fun testNoCommandsToQueryAgc() {
        // this tests verifies an exception is thrown when no commands are passed to queryAgc()
        val command = mutableListOf<String>()
        assertThrows<IllegalStateException> {
            //Check that an error is thrown when the commands list is empty
            queryAgc(command.toTypedArray())
        }
    }

    @Test
    fun testRetrieveAgcGenomes() {
        val dbPath = "${TestExtension.testTileDBURI}" // just the path, code appends "assemblies.agc"
        var genomeList = mutableListOf<String>()
        val gn1 = "LineA"

        genomeList.add(gn1)
        var agcResult = retrieveAgcGenomes(dbPath, genomeList)

        assertEquals(2, agcResult.size)
        assertEquals(2, agcResult.keys.size)
        val fastaFile = File("${testOutputFastaDir}/LineA.fa").toString()
        val lineA = fastaToNucSeq(fastaFile)

        // Check the sequences match
        assertEquals(lineA["1"]!!.toString(), agcResult[Pair("LineA","1")]!!.toString())
        assertEquals(lineA["2"]!!.toString(), agcResult[Pair("LineA","2")]!!.toString())

        // try again with 2 genomes on the list.  This should now work
        // as we have updated the idline to contain the samplename
        val gn2 ="LineC"
        genomeList.add(gn2)
        agcResult = retrieveAgcGenomes(dbPath, genomeList)
        assertEquals(4, agcResult.size)
        assertEquals(4, agcResult.keys.size)

        //TODO LCJ - add test to verify the genome names are included in the results
        val keys = agcResult.keys
        assertTrue(keys.contains(Pair("LineA","1")))
        assertTrue(keys.contains(Pair("LineA","2")))
        assertTrue(keys.contains(Pair("LineC","1")))
        assertTrue(keys.contains(Pair("LineC","2")))
    }

    @Test
    fun testRetrieveAgcData_listset() {
        // This command tests the retrieveAgcData() function with a listset
        // the list is valid, and the listset is valid
        val dbPath = "${TestExtension.testTileDBURI}" // just the path, code appends "assemblies.agc"
        val commands = mutableListOf<String>("listset")

        var agcResult = retrieveAgcData(dbPath, commands)
        val expectedResult = listOf("LineA", "LineB","LineC","Ref")
        assertEquals(expectedResult, agcResult)

        // Verify an exception is thrown is an invalid command is sent to retrieveAgcData()
        assertThrows<IllegalStateException> {
            //Check that an error is thrown if the command is invalid
            retrieveAgcData(dbPath, listOf("happy"))
        }
    }

    @Test
    fun testRetrieveAgcData_listctg() {
        // This command tests the retrieveAgcData() function with a listctg command

        // Test with bad dbPath
        val badDbPath = "/my/bad/dbPath"
        assertThrows<IllegalStateException> {
            //Check that an error is thrown if the dbPath is invalid
            retrieveAgcData(badDbPath, listOf("listctg"))
        }

        // First, test without a geomome list
        val dbPath = "${TestExtension.testTileDBURI}" // just the path, code appends "assemblies.agc"
        val commands = mutableListOf<String>()

        // Verify an exception is thrown if no commands are on the list
        assertThrows<IllegalStateException> {
            //Check that an error is thrown if no cammand  is given
            retrieveAgcData(dbPath, commands)
        }

        commands.add("listctg")
        assertThrows<IllegalStateException> {
            //Check that an error is thrown if no genomes are included in the query.
            retrieveAgcData(dbPath, commands)
        }

        // Now test with a genome list
        commands.add("LineA")
        commands.add("LineB")
        val agcResult = retrieveAgcData(dbPath, commands)
        println(" agcResult = $agcResult")
        assertEquals(2, agcResult!!.size)

        println("Results from agcResult")
        for (entry in agcResult) {
            println("entry =${entry}")
        }
        // Verify agcResult contains the entry"LineA:1 sampleName=LineA,2 sampleName=LineA"
        // Remember - AGC returns the full idLine, including all comments
        assertTrue(agcResult.contains("LineA:1 sampleName=LineA,2 sampleName=LineA"))
        // Verify agcResult contains the entry"LineB:1 sampleName=LineB,2 sampleName=LineB"
        assertTrue(agcResult.contains("LineB:1 sampleName=LineB,2 sampleName=LineB"))

    }

    @Test
    fun testVerifySampleNameBad() {
        // first file has good annotations, second file is missing "sampleName="
        val fastaList = listOf("data/test/smallSeq/LineA.fa", "data/test/agcTestBad/LineA_noSN.fa")
        assertThrows<IllegalStateException> {
            //Check that an exception is thrown when the idline does not contain "sampleName="
            // the first file is fine, the second file is missing "sampleName="
            // Manually verified that both files were read, and the code went to the second file
            // after the first line of the first file was verified as good.
            val verifiedResults = AgcCompress().verifyFileAnnotation(fastaList)
        }
    }

    @Test
    fun testVerifySampleNameGood() {
        // 2 files, both are properly annotated.
        val badFastaList = listOf("data/test/smallSeq/LineC.fa", "data/test/smallSeq/Ref.fa")
        val goodResults = AgcCompress().verifyFileAnnotation(badFastaList)
        println("verifiedResults = $goodResults")
        assertEquals(true, goodResults)
    }
}