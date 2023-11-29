package net.maizegenetics.phgv2.agc

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class AnnotateFastaTest {
    companion object {

        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"
        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.testOutputFastaDir).mkdirs()
        }

//        @JvmStatic
//        @AfterAll
//        fun teardown() {
//            File(TestExtension.tempDir).deleteRecursively()
//        }
    }

    @Test
    fun testCliktParams() {
        val annotateFasta = AnnotateFasta()

        // Test missing fasta-list parameter
        val resultMissingFastaList = annotateFasta.test(" --output-dir ${TestExtension.testOutputFastaDir}")
        assertEquals(resultMissingFastaList.statusCode, 1)
        assertEquals("Usage: annotate-fasta [<options>]\n" +
                "\n" +
                "Error: invalid value for --fasta-list: --fasta-list must not be blank\n",resultMissingFastaList.output)

        // Test missing output-dir parameter
        val resultMissingOutDir = annotateFasta.test("--fasta-list ${TestExtension.testInputFastaDir} ")
        assertEquals(resultMissingOutDir.statusCode, 1)
        assertEquals("Usage: annotate-fasta [<options>]\n" +
                "\n" +
                "Error: invalid value for --output-dir: --output-dir must not be blank\n",resultMissingOutDir.output)

    }

    @Test
    fun testAnnotateFastaCommand() {
        val fastaInputDir = "data/test/smallseq"
        val fastaOutputDir = TestExtension.testOutputFastaDir

        // Create a List<String> of fasta files in the fastaInputDir
        val fileList = File(fastaInputDir).listFiles().filter { it.extension == "fa" || it.extension == "fasta" }.map { it.absolutePath }

        // Create a file named ${fastaOutputDir}/inputFastaFiles containing the full path name for each fasta file in the fastaInputDir
        // this should be written to the fastaOutputDir
        val filesToUpdate = File(fastaOutputDir, "fastaCreateFileNames.txt")
        filesToUpdate.writeText(fileList.joinToString("\n"))

        // Test the AnnotateFasta class
        val annotateFasta = AnnotateFasta()
        val result = annotateFasta.test( "--fasta-list ${filesToUpdate} --threads 2 --output-dir ${TestExtension.testOutputFastaDir}")
        assertEquals(result.statusCode, 0)

        // get a list of fasta files created in the fastaOutputDir, as a List<String> in variable named updatedFiles
        val updatedFiles = File(fastaOutputDir).listFiles().filter { it.extension == "fa" || it.extension == "fasta" }.map { it.absolutePath }

        // verify the idlines of each fasta file were updated to include
        // "sampleName=${sampleName}" where sampleName is the fasta file name minus the extension
        updatedFiles.forEach { fastaFile ->
            val sampleName = File(fastaFile).nameWithoutExtension
            val newFilename = "${fastaOutputDir}/${File(fastaFile).name}"
            File(newFilename).forEachLine { line ->
                if (line.startsWith(">")) {
                    assertTrue(line.contains("sampleName=${sampleName}"))
                }
            }
        }

        // run this again on the newly updated files - verify the idlines are not updated if sampleName= is already present
        // Create a new fasta list file with the updated fasta files
        val filesToUpdate2 = File(fastaOutputDir, "fastaCreateFileNames2.txt")
        filesToUpdate2.writeText(updatedFiles.joinToString("\n"))
        // run the annotateFasta command again on the updated fasta files.  This will overwrite the existing files
        // in the outputDir
        val secondOutputDir = "${TestExtension.testOutputFastaDir}/secondOutputDir"
        // create the new outputDir
        File(secondOutputDir).mkdirs()
        val result2 = annotateFasta.test( "--fasta-list ${filesToUpdate2} --output-dir ${secondOutputDir}")

        // Get list of fastas files in the newOutputDir
        val updatedFiles2 = File(secondOutputDir).listFiles().filter { it.extension == "fa" || it.extension == "fasta" }.map { it.absolutePath }
        assertEquals(result2.statusCode, 0)
        // verify the idlines lines of each fasta files contain only a single "sampleName=${sampleName}" string
        updatedFiles2.forEach { fastaFile ->
            val sampleName = File(fastaFile).nameWithoutExtension
            val newFilename = "${fastaOutputDir}/${File(fastaFile).name}"
            File(newFilename).forEachLine { line ->
                if (line.startsWith(">")) {
                    assertTrue(line.contains("sampleName=${sampleName}"))
                    assertTrue(line.indexOf("sampleName=${sampleName}") == line.lastIndexOf("sampleName=${sampleName}"))
                }
            }
        }

        // verify that each file in the updatedFiles list is the same as the corresponding file in the updatedFiles2 list
        updatedFiles.forEachIndexed { index, fastaFile ->

            val newFilename = "${fastaOutputDir}/${File(fastaFile).name}"
            val newFilename2 = "${secondOutputDir}/${File(fastaFile).name}"

            // compare the contents of the two files
            val file1 = File(newFilename).readLines()
            val file2 = File(newFilename2).readLines()
            assertEquals(file1, file2)
        }
    }
}