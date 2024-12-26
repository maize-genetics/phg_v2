package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.AgcCompress
import net.maizegenetics.phgv2.cli.CreateFastaFromHvcf
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter

class RopebwtIndexTest {


    companion object {
        val tempTestDir = "${TestExtension.tempDir}ropebwtTest/"
        val tempDBPathDir = "${TestExtension.testOutputFastaDir}dbPath/"


        //Setup/download  files
        //Resetting on both setup and teardown just to be safe.
        @JvmStatic
        @BeforeAll
        fun setup() {
            resetDirs()
            setupAgc()
            buildPangenomeFile()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            resetDirs()
        }

        private fun resetDirs() {

            File(TestExtension.tempDir).deleteRecursively()
            File(TestExtension.testOutputFastaDir).deleteRecursively()
            File(TestExtension.testOutputDir).deleteRecursively()
            File(tempTestDir).deleteRecursively()
            File(tempDBPathDir).deleteRecursively()

            File(TestExtension.tempDir).mkdirs()
            File(TestExtension.testOutputFastaDir).mkdirs()
            File(TestExtension.testOutputDir).mkdirs()
            File(tempTestDir).mkdirs()
            File(tempDBPathDir).mkdirs()


        }

        private fun setupAgc() {
            //create an AGC record with the Ref in it
            val altFileListFile = TestExtension.testOutputFastaDir+"/agc_altList.txt"
            BufferedWriter(FileWriter(altFileListFile)).use { writer ->
                writer.write("data/test/smallseq/LineA.fa\n")
                writer.write("data/test/smallseq/LineB.fa\n")
                writer.write("data/test/smallseq/Ref.fa\n")
            }

            val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"
            File(dbPath).mkdirs()

            //Call AGCCompress to create the AGC file
            val agcCompress = AgcCompress()
            agcCompress.processAGCFiles(dbPath,altFileListFile,"data/test/smallseq/Ref.fa","")
        }

        private fun buildPangenomeFile() {
            val createFasta = CreateFastaFromHvcf()
            createFasta.test("--db-path ${TestExtension.testOutputFastaDir}/dbPath --output-dir ${tempTestDir} --fasta-type pangenomeHaplotype --hvcf-dir data/test/smallseq/")
        }
    }

    @Test
    fun testCliktParams() {
        val ropebwtIndex = RopebwtIndex()
        ropebwtIndex.parse(arrayOf("--input-fasta","data/test/smallseq/Ref.fa","--index-file-prefix","${TestExtension.tempDir}ropebwtTest/testIndex","--num-threads","3","--delete-fmr-index","--conda-env-prefix",""))
        ropebwtIndex.run()
    }
}