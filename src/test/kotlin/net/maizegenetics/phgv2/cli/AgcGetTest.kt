package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension.Companion.testTileDBURI
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import java.io.BufferedReader
import java.io.File
import java.nio.file.Files
import java.nio.file.Paths
import kotlin.test.assertEquals

class AgcGetTest {
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
            File(fastaOutputDir).walk(FileWalkDirection.TOP_DOWN).filter{it.name.endsWith(".fa")}.forEach {
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

//        @JvmStatic
//        @AfterAll
//        fun teardown() {
//            File(TestExtension.tempDir).deleteRecursively()
//        }
    }

    @Test
    fun testCliktParams() {
        val agcGet = AgcGet()

        val sampleNames = "B73,KY27"
        val contig = "1"
        val start = 20
        val end = 40

        // Test missing db-path parameter
        val resultMissingDB = agcGet.test("--contigs ${contig} --sample-names ${sampleNames} --start 1 --end 20")
        assertEquals(1, resultMissingDB.statusCode )
        assertEquals("Usage: agc-get [<options>]\n" +
                "\n" +
                "Error: invalid value for --db-path: --db-path must not be blank\n",resultMissingDB.output)

        // Test missing db-path parameter
        val resultMissingSampleNames = agcGet.test("--db-path ${TestExtension.testTileDBURI}  --contigs ${contig} --start 1 --end 20")
        assertEquals(resultMissingSampleNames.statusCode,1 )
        assertEquals("Usage: agc-get [<options>]\n" +
                "\n" +
                "Error: invalid value for --sample-names: --sample-names must not be blank\n",resultMissingSampleNames.output)

        // Test missing end parameter
        val resultMissingEnd = agcGet.test("--db-path ${TestExtension.testTileDBURI} --contigs ${contig} --sample-names ${sampleNames} --start 1")
        assertEquals(1,resultMissingDB.statusCode)
        assertEquals("Usage: agc-get [<options>]\n" +
                "\n" +
                "Error: invalid value for --end: if --start is present, end must be greater than or equal to --start\n",resultMissingEnd.output)


        // Test missing start parameter
        // This passes, even when "end" parameter is present.  That is because
        // without creating an OptionsGroup, these cannot be linked or I get a
        // "recursive problem" error.  But when an options group is created and you make them
        // cooccurring (ie both must be there), then one of the 2 has to be marked
        // as required, and in this case, neither is required.  They are only required
        // to be present as a pair.
        assertThrows<IllegalArgumentException> {
            //Check that an error is thrown when the end param is present but start is missing
            agcGet.test("--db-path ${TestExtension.testTileDBURI} --contigs ${contig} --sample-names ${sampleNames} --end 20")
        }

        assertThrows<IllegalArgumentException> {
            //Check that an error is thrown when the start/end params are present but contigs is missing
            agcGet.test("--db-path ${TestExtension.testTileDBURI}  --sample-names ${sampleNames} --start 2 --end 20")
        }

    }

    @Test
    fun testGoodCall() {
        val agcGet = AgcGet()

        val dbPath = "${testTileDBURI}/assemblies.agc"
        val sampleNames = "LineA,LineB"
        val contig = "1"
        val start = 20
        val end = 40
        // The assemblies agc file was created in the setup() above

        println("testGOodCall - calling agcGet.test")
        val command = agcGet.buildAgcCommandFromParams(dbPath,sampleNames, contig, start, end)
        assertEquals(command.size,9)

        val agcResult = agcGet.retrieveAgcData(command)
        val content = agcResult.bufferedReader().use(BufferedReader::readText)

        assertTrue(content.contains(">1:20-40"))

        //val content = agcGet.test("--db-path ${dbPath}  --sample-names ${sampleNames} --contig ${contig}--start 2 --end 20")
        println("testGoodCall - agcResult = \n${content}")
    }
}