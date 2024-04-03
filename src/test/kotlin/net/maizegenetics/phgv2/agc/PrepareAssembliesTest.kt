package net.maizegenetics.phgv2.agc

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class PrepareAssembliesTest {
    companion object {

        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"
        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.testInputFastaDir).mkdirs()
            File(TestExtension.testOutputFastaDir).mkdirs()

            val fastaOrigDir = "data/test/smallseq"
            val fastaInputDir = TestExtension.testInputFastaDir
            // copy files with extension .fa from data/test/smallseq to fastaInputDir for this test
            val fastaFiles = File(fastaOrigDir).listFiles { file -> file.extension == "fa" }
            fastaFiles.forEach { file -> file.copyTo(File(fastaInputDir, file.name)) }

        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testCliktParams() {
        val prepareAssemblies = PrepareAssemblies()

        // Test missing fasta-list parameter
        val resultMissingKeyfile = prepareAssemblies.test(" --output-dir ${TestExtension.testOutputFastaDir}")
        assertEquals(resultMissingKeyfile.statusCode, 1)
        assertEquals("Usage: prepare-assemblies [<options>]\n" +
                "\n" +
                "Error: invalid value for --keyfile: --keyfile must not be blank\n",resultMissingKeyfile.output)

        // Test missing output-dir parameter
        val resultMissingOutDir = prepareAssemblies.test("--keyfile ${TestExtension.testInputFastaDir} ")
        assertEquals(resultMissingOutDir.statusCode, 1)
        assertEquals("Usage: prepare-assemblies [<options>]\n" +
                "\n" +
                "Error: invalid value for --output-dir: --output-dir must not be blank\n",resultMissingOutDir.output)

    }

    @Test
    fun testPrepareAssembliesCommand() {
        val fastaInputDir = TestExtension.testInputFastaDir
        val fastaOutputDir = TestExtension.testOutputFastaDir

        // Create a List<String> of fasta files in the fastaInputDir
        val fileList = File(fastaInputDir).listFiles().filter { it.extension == "fa" || it.extension == "fasta" }.map { it.absolutePath }

        // Create a tab-delimited file of fasta file names and sample names
        // The fasta file names are the full path names for the fasta files in the fastaInputDir
        // the sample names are the fasta file names minus the extension
        // write the fasta file names in the first column and the sample names in the second column of a tab-delimited file
        // named ${fastaOutputDir}/fastaCreateFileNames.txt
        val filesToUpdate = File(fastaOutputDir, "fastaCreateFileNames.txt")
        filesToUpdate.writeText(fileList.joinToString("\n") { "${it}\t${File(it).nameWithoutExtension}" })

        // Test the AnnotateFasta class
        val prepareAssemblies = PrepareAssemblies()
        val result = prepareAssemblies.test( "--keyfile ${filesToUpdate} --threads 2 --output-dir ${TestExtension.testOutputFastaDir}")
        assertEquals(0, result.statusCode )

        // get a list of fasta files created in the fastaOutputDir, as a List<String>
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
        filesToUpdate2.writeText(fileList.joinToString("\n") { "${it}\t${File(it).nameWithoutExtension}" })

        // run the annotateFasta command again on the updated fasta files.  This will overwrite the existing files
        // in the outputDir
        val secondOutputDir = "${TestExtension.testOutputFastaDir}/secondOutputDir"
        File(secondOutputDir).mkdirs()
        val result2 = prepareAssemblies.test( "--keyfile ${filesToUpdate2} --output-dir ${secondOutputDir}")

        // Get list of fastas files in the newOutputDir
        val updatedFiles2 = File(secondOutputDir).listFiles().filter { it.extension == "fa" || it.extension == "fasta" }.map { it.absolutePath }
        assertEquals(0,result2.statusCode)
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

    @Test
    fun testChangedFastaName() {
        // This test verifies that the fasta file name is changed to the sample name
        // Use for testing the fasta file in folder data/test/annotateFastaTest
        // make a new directory as a subdirectory of TestExtension.tempdir and call it annotatedFastaTest
        // copy the fasta files from data/test/annotateFastaTest to the new directory

        val fastaOrigDir = "data/test/annotateFastaTest"
        val fastaInputDir = "${TestExtension.tempDir}/annotatedFastaTest"
        val fastaOutputDir = "${TestExtension.tempDir}/annotatedFastaTestOutput"

        File(fastaInputDir).mkdirs()
        File(fastaOutputDir).mkdirs()

        val fastaFiles = File(fastaOrigDir).listFiles { file -> file.extension == "fa" || file.extension == "fasta" || file.extension == "gz"}
        fastaFiles.forEach { file -> file.copyTo(File(fastaInputDir, file.name)) }

        // Create the key file for this fasta files
        val annotateKeyFile = File(fastaOutputDir, "fastaCreateFileNames.txt")
        // The file name is first, followed by the sampleName.
        // there is only 1 *.fa file and 1 *.fa.gz file, so
        // writing the keyfile with known values for the sample names.

        // Key file has 2 columns: the file name and the sample name
        BufferedWriter(FileWriter(annotateKeyFile)).use { writer ->
            for (file in fastaFiles) {
                if (file.extension == "gz") {
                    // THere is only 1 file compressed, that is LineB
                    writer.write("${file.absolutePath}\tLineB\n")
                } else  {
                    // THere is only 1 file not compressed, that is LineA
                    writer.write("${file.absolutePath}\tLineA\n")
                }

            }
        }

        // Run annotateFasta on the fasta files in the fastaInputDir
        val prepareAssemblies = PrepareAssemblies()
        val result = prepareAssemblies.test( "--keyfile ${annotateKeyFile} --threads 2 --output-dir ${fastaOutputDir}")
        assertEquals(0, result.statusCode)

        // Verify the id lines of each fasta file were updated to include
        // "sampleName=${sampleName}" where sampleName is the fasta file name minus the extension
        // And verify the output fasta is named LineA.fa
        //val fastaFiles = File(fastaOrigDir).listFiles { file -> file.extension == "fa" || file.extension == "fasta" || file.endsWith(".gz")}
        val updatedFiles = File(fastaOutputDir).listFiles { file -> file.extension == "fa" || file.endsWith(".gz")  }.map { it.absolutePath }

        updatedFiles.forEach { fastaFile ->
            val sampleName = File(fastaFile).nameWithoutExtension.substringBefore(".")
            val newFilename = "${fastaOutputDir}/${File(fastaFile).name}"
            File(newFilename).forEachLine { line ->
                if (line.startsWith(">")) {
                    assertTrue(line.contains("sampleName=${sampleName}"))
                }
            }
        }
        // verify files named LineA.fa and LineB.fa exist in the output directory
        // We no longer compress the fasta files, so the output files are named LineA.fa and LineB.fa
        // This test verifies the created annotated file for LineB is uncompressed though the original file was compressed.
        assertTrue(File(fastaOutputDir, "LineA.fa").exists())
        assertTrue(File(fastaOutputDir, "LineB.fa").exists())

    }
    @Test
    fun testAnnotateGzippedFastaCommand() {
        // This test is the same as above, but with gzipped files
        // It verifies we can read and write gzipped files and they are named properly
        val fastaInputDir = TestExtension.testInputFastaDir
        val fastaOutputDir = TestExtension.testOutputFastaDir

        // Create a List<String> of fasta files in the fastaInputDir
        var fileList = File(fastaInputDir).listFiles().filter { it.extension == "fa" || it.extension == "fasta" }.map { it.absolutePath }

        // run gzip on all the files in the fileList
        fileList.forEach { file ->
            val cmd = "gzip ${file}"
            val proc = Runtime.getRuntime().exec(cmd)
            proc.waitFor()
        }

        // Update the fileList with the gzipped names
        fileList = File(fastaInputDir).listFiles().filter { it.extension == "gz" }.map { it.absolutePath }

        // Create a tab-delimited file of fasta file names and sample names
        // The fasta file names are the full path names for the fasta files in the fastaInputDir
        // the sample names are the fasta file names just up to the first "." character
        // write the fasta file names in the first column and the sample names in the second column of  a tab-delimited
        // file  named ${fastaOutputDir}/fastaCreateFileNames.txt
        val filesToUpdate = File(fastaOutputDir, "fastaCreateFileNames.txt")
        // write the fasta file names in the first column and the sample names in the second column of  a tab-delimited
        filesToUpdate.writeText(fileList.joinToString("\n") { "${it}\t${File(it).nameWithoutExtension.substringBefore(".")}" })

        println("write keyfile to ${filesToUpdate.absolutePath}")

        // Test the AnnotateFasta class
        val prepareAssemblies = PrepareAssemblies()
        val result = prepareAssemblies.test( "--keyfile ${filesToUpdate} --threads 2 --output-dir ${TestExtension.testOutputFastaDir}")
        assertEquals(0,result.statusCode)

        // get a list of fasta files created in the fastaOutputDir, as a List<String> in variable named updatedFiles
        val updatedFiles = File(fastaOutputDir).listFiles().filter { it.extension == ".gz" }.map { it.absolutePath }

        // verify the idlines of each fasta file were updated to include
        // "sampleName=${sampleName}" where sampleName is the fasta file name minus the extension
        // These are gzipped files, so you must read the compressed file and decompress it to get the fasta file

        updatedFiles.forEach { fastaFile ->
            val sampleName = File(fastaFile).nameWithoutExtension.substringBefore(".")
            val newFilename = "${fastaOutputDir}/${File(fastaFile).name}"
            File(newFilename).forEachLine { line ->
                if (line.startsWith(">")) {
                    assertTrue(line.contains("sampleName=${sampleName}"))
                }
            }
        }

        println("Running the second time")
        // run this again on the newly updated files - verify the idlines are not updated if sampleName= is already present
        // Create a new fasta list file with the updated fasta files
        val filesToUpdate2 = File(fastaOutputDir, "fastaCreateFileNames2.txt")
        filesToUpdate2.writeText(fileList.joinToString("\n") { "${it}\t${File(it).nameWithoutExtension.substringBefore(".")}" })

        // run the annotateFasta command again on the updated fasta files.  This will overwrite the existing files
        // in the outputDir
        val secondOutputDir = "${TestExtension.testOutputFastaDir}/secondOutputDir"
        // create the new outputDir
        File(secondOutputDir).mkdirs()
        val result2 = prepareAssemblies.test( "--keyfile ${filesToUpdate2} --output-dir ${secondOutputDir}")

        // Get list of fastas files in the newOutputDir
        val updatedFiles2 = File(secondOutputDir).listFiles().filter { it.extension == ".gz" }.map { it.absolutePath }
        assertEquals(0,result2.statusCode )
        // verify the idlines lines of each fasta files contain only a single "sampleName=${sampleName}" string
        updatedFiles2.forEach { fastaFile ->
            val sampleName = File(fastaFile).nameWithoutExtension.substringBefore(".")
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