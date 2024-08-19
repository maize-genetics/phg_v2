package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import java.io.File
import java.nio.file.Files
import java.nio.file.Paths
import java.time.LocalDate
import kotlin.test.assertEquals

class AgcCompressTest {
    companion object {

        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"
        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.testTileDBURI).mkdirs()
            Initdb().createDataSets(TestExtension.testTileDBURI,"",false)
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testCliktParams() {
        val agcCompress = AgcCompress()
        val refFasta = "data/test/smallseq/Ref.fa"
        // Test missing fasta-dir parameter
        val resultMissingFastaList = agcCompress.test("--db-path ${TestExtension.testTileDBURI} --reference-file ${refFasta}")
        assertEquals(resultMissingFastaList.statusCode, 1)
        assertEquals("Usage: agc-compress [<options>]\n" +
                "\n" +
                "Error: missing option --fasta-list\n",resultMissingFastaList.output)


        // Test missing refFasta parameter
        val resultRefFasta = agcCompress.test("--fasta-list ${TestExtension.testInputFastaDir} --db-path ${TestExtension.testTileDBURI}")
        assertEquals(resultRefFasta.statusCode, 1)
        assertEquals("Usage: agc-compress [<options>]\n" +
                "\n" +
                "Error: missing option --reference-file\n",resultRefFasta.output)
    }

    @Test
    fun testAgcCommand() {
        // Test the AgcCompress class, both create and append to compressed file
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
        File(fastaOutputDir).walk(FileWalkDirection.TOP_DOWN).filter{it.name.endsWith(".fa")}.forEach {
            fileList.add(it.toString())
        }

        // create a list of files to append - just LineC for now
        val appendFileList = mutableListOf<String>()
        fileList.forEach{
            if (it.contains("LineC.fa")) appendFileList.add(it)
        }

        // print the fileList and appendFileList contents
        println("fileList contents: \n${fileList.joinToString("\n")}")
        println("appendFileList contents: \n${appendFileList.joinToString("\n")}")

        // The list of fasta files to load should not contain the reference.  The reference is listed separately.
        // LineC is also excluded as we will use that in the append command.
        fileList.removeIf { it.contains("Ref") || it.contains("LineC") }

        println("\nfileList contents after removing Ref and LineC: \n${fileList.joinToString("\n")}")

        // write the full path fasta file names  from the fastaOutputDir to a single file, one per line, in tempDir
        // This file will be used as input to the agc-compress command
        val fastaCreateFileNamesFile = File(dbPath, "fastaCreateFileNames.txt")
        fastaCreateFileNamesFile.writeText(fileList.joinToString("\n"))

        val fastaAppendFileNamesFile = File(dbPath, "fastaAppendFileNames.txt")
        fastaAppendFileNamesFile.writeText(appendFileList.joinToString("\n"))

        // print the contents of the fastaFileNamesFile
        println("fastaFileNamesFile contents: \n${fastaCreateFileNamesFile.readText()}")
        println("fastaAppendFileNamesFile contents: \n${fastaAppendFileNamesFile.readText()}")

        val refFasta = File(fastaOutputDir, "Ref.fa").toString()

        val agcCompress = AgcCompress()
        // Create the initial compressed file
        println("Calling agcCompress for CREATE")
        var result = agcCompress.test("--fasta-list ${fastaCreateFileNamesFile} --db-path ${dbPath} --reference-file ${refFasta}")
        assertEquals(result.statusCode, 0)

        // Verify that file dbPath/assemblies.agc exists
        val agcFile = File(dbPath, "assemblies.agc")
        assertEquals(true, agcFile.exists())

        // Verify the samples in the agc file are as expected,
        val agcInitialSampleList = agcCompress.getSampleListFromAGC(agcFile.toString(),"${dbPath}/temp","",false)
        val initialSampleSet = listOf("LineA","LineB","Ref")
        assertEquals(true, agcInitialSampleList.containsAll(initialSampleSet))

        // Append LineC to the compressed file
        println("Calling agcCompress for APPEND")
        result = agcCompress.test("--fasta-list ${fastaAppendFileNamesFile.toString()} --db-path ${dbPath} --reference-file ${refFasta}")
        assertEquals(result.statusCode, 0)

        // Verify that file dbPath/assemblies.agc exists
        // verify that file dbPath/assemblies_backup_${LocalDate.now()}.agc exists
        // assemblies.agc should be the updated file, assemblies_backup_${LocalDate.now()}.agc should be the original
        val agcBackupFile = File(dbPath, "assemblies_backup_${LocalDate.now()}.agc")
        assertEquals(true, agcFile.exists())
        assertEquals(true, agcBackupFile.exists())

        // Verify the samples in both agc files are as expected,
        // do this using the agcCompress getSampleListFromAGC method
        val agcSampleList = agcCompress.getSampleListFromAGC(agcFile.toString(),"${dbPath}/temp","",false)
        val sampleSetForAGC = listOf("LineA","LineB","LineC","Ref") // all files plus the file (LineC) that was appended
        assertEquals(true, agcSampleList.containsAll(sampleSetForAGC))

        val agcBackupSampleList = agcCompress.getSampleListFromAGC(agcBackupFile.toString(),"${dbPath}/temp","",false)
        // Verify the agcBackupSampleLIst contains only the samples in fileList and does not
        // cantain the samples in appendFileList
        val sampleSetForBackup = listOf("LineA","LineB","Ref") // all files minus the file (LineC) that was appended
        assertEquals(true, agcBackupSampleList.containsAll(sampleSetForBackup))

        // Call agcCompress once more with LineC.fa in the fastaAppendFileNamesFile
        // This should return without loading anything.  A message to that effect will be printed
        // to the console.

        println("Calling second agcCompress for APPEND")
        result = agcCompress.test("--fasta-list ${fastaAppendFileNamesFile.toString()} --db-path ${dbPath} --reference-file ${refFasta}")
        assertEquals(result.statusCode, 0)

    }

    @Test
    fun testBadAGCOption() {
        // agc load will only accept "create" or "append" as options
        val fastaList = "FastaListFile.txt"
        val refFasta = "Ref.fa"
        val dbPath = TestExtension.testTileDBURI

        val agcCompress = AgcCompress()
         // "badOption" is not valid option, should return false
        val success = agcCompress.loadAGCFiles(fastaList, "badOption",dbPath,refFasta, tempDir,"",false)
        assertEquals(false, success)
    }

    @Test
    fun testAgcCompressNoSampleName() {
        //This test is to verify queryAgc() throws an exception when there is no "sampleName=" in the idline

        val dbPath = TestExtension.tempDir
        val refFasta = "data/test/smallseq/Ref.fa"

        val fastaCreateFileNamesFile = File(dbPath, "fastaBadNames.txt")
        fastaCreateFileNamesFile.writeText("data/test/agcTestBad/LineA_noSN.fa\n")

        Initdb().createDataSets(TestExtension.tempDir,"",false)
        val agcCompress = AgcCompress()
        assertThrows<IllegalStateException> {
            // Create the initial compressed file
            val agcCompressResult = agcCompress.test("--fasta-list ${fastaCreateFileNamesFile} --db-path ${dbPath} --reference-file ${refFasta}")
        }

    }
}