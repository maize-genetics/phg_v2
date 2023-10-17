package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import java.lang.Exception
import java.nio.file.Files
import java.nio.file.Paths
import java.util.stream.Collectors
import kotlin.test.assertEquals

class AgcCompressTest {
    companion object {

        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"
        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.tempDir).mkdirs()
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
        val resultMissingFastaList = agcCompress.test("--db-path ${TestExtension.testTileDBURI} --ref-fasta ${refFasta}")
        assertEquals(resultMissingFastaList.statusCode, 1)
        assertEquals("Usage: agc-compress [<options>]\n" +
                "\n" +
                "Error: invalid value for --fasta-list: --fasta-list must not be blank\n",resultMissingFastaList.output)

        // Test missing db-path parameter
        val resultMissingDB = agcCompress.test("--fasta-list ${TestExtension.testInputFastaDir} --ref-fasta ${refFasta}")
        assertEquals(resultMissingDB.statusCode, 1)
        assertEquals("Usage: agc-compress [<options>]\n" +
                "\n" +
                "Error: invalid value for --db-path: --db-path must not be blank\n",resultMissingDB.output)

        // Test missing refFasta parameter
        val resultRefFasta = agcCompress.test("--fasta-list ${TestExtension.testInputFastaDir} --db-path ${TestExtension.testTileDBURI}")
        assertEquals(resultRefFasta.statusCode, 1)
        assertEquals("Usage: agc-compress [<options>]\n" +
                "\n" +
                "Error: invalid value for --ref-fasta: --ref-fasta must not be blank\n",resultRefFasta.output)
    }

    @Test
    fun testAgcCommand() {
        // this method is just to test that the command runs without error
        val fastaInputDir = "data/test/smallseq"
        val fastaOutputDir = TestExtension.testOutputFastaDir
        val dbPath = TestExtension.testTileDBURI

        // WHy isn't this handled in TestExtension - I have added the command for TestExtension.testTileDBURI
        Files.createDirectories(Paths.get(dbPath))

        // copy files with extension .fa from data/test/smallseq to fastaOutputDir
        val fastaFiles = File(fastaInputDir).listFiles { file -> file.extension == "fa" }
        fastaFiles.forEach { file -> file.copyTo(File(fastaOutputDir, file.name)) }

        val fileList = mutableListOf<String>()
        val appendFileList = mutableListOf<String>()
        // get the full path fasta file names  from the fastaInputDir
        File(fastaOutputDir).walk(FileWalkDirection.TOP_DOWN).filter{it.name.endsWith(".fa")}.forEach {
            fileList.add(it.toString())
        }

        // create a list of files to append - just LineC for now
        fileList.forEach{
            if (it.contains("LineC.fa")) appendFileList.add(it)
        }

        // The list of fasta files to load should not contain the reference.  The reference is listed separately.
        // LineC is also excluded as we will use that in the append command.
        fileList.removeIf { it.contains("Ref") || it.contains("LineC") }

        // write the full path fasta file names  from the fastaOutputDir to a single file, one per line, in tempDir
        val fastaCreateFileNamesFile = File(dbPath, "fastaCreateFileNames.txt")
        fastaCreateFileNamesFile.writeText(fileList.joinToString("\n"))

        val fastaAppendFileNamesFile = File(dbPath, "fastaAppendFileNames.txt")
        fastaAppendFileNamesFile.writeText(appendFileList.joinToString("\n"))

        // print the contents of the fastaFileNamesFile
        println("fastaFileNamesFile contents: \n${fastaCreateFileNamesFile.readText()}")
        println("fastaAppendFileNamesFile contents: \n${fastaAppendFileNamesFile.readText()}")

        val refFasta = File(fastaOutputDir, "Ref.fa").toString()

        val agcCompress = AgcCompress()
        println("Calling agcCompress for CREATE")
        var result = agcCompress.test("--fasta-list ${fastaCreateFileNamesFile.toString()} --db-path ${dbPath} --ref-fasta ${refFasta}")

        println("Calling agcCompress for APPEND")
        result = agcCompress.test("--fasta-list ${fastaAppendFileNamesFile.toString()} --db-path ${dbPath} --ref-fasta ${refFasta}")

    }
}