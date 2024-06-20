package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.AgcCompress
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter

class SimulateReadsTest {

    companion object {
        //Setup/download  files
        //Resetting on both setup and teardown just to be safe.
        @JvmStatic
        @BeforeAll
        fun setup() {
            resetDirs()
            setupAgc()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            resetDirs()
        }

        fun resetDirs() {
            File(TestExtension.tempDir).deleteRecursively()

            File(TestExtension.testVCFDir).mkdirs()
            File(TestExtension.testInputFastaDir).mkdirs()
            File(TestExtension.testOutputFastaDir).mkdirs()
            File(TestExtension.testOutputDir).mkdirs()
        }

        private fun setupAgc() {
            //create an AGC record with the Ref in it
            val altFileListFile = TestExtension.testOutputFastaDir + "/agc_altList.txt"
            BufferedWriter(FileWriter(altFileListFile)).use { writer ->
                writer.write("data/test/smallseq/LineA.fa\n")
                writer.write("data/test/smallseq/LineB.fa\n")
                writer.write("data/test/smallseq/Ref.fa\n")
            }

            val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"
            File(dbPath).mkdirs()

            //Call AGCCompress to create the AGC file
            val agcCompress = AgcCompress()
            agcCompress.processAGCFiles(dbPath,altFileListFile,"data/test/smallseq/Ref.fa")

            //copy hvcf files to TestExtension.testVCFDir
            File(TestExtension.smallseqLineAHvcfFile).copyTo(File("${TestExtension.testVCFDir}LineA.h.vcf"))
            File(TestExtension.smallseqLineBHvcfFile).copyTo(File("${TestExtension.testVCFDir}LineB.h.vcf"))
            File(TestExtension.smallseqRefHvcfFile).copyTo(File("${TestExtension.testVCFDir}Ref.h.vcf"))
        }
    }

    @Test
    fun testSimulateHaploid() {
        println("Simulating haploid reads")
        val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"
        val haploidArgs = "--db-path $dbPath --hvcf-dir ${TestExtension.testVCFDir} --sample-names LineA --fastq-out-dir ${TestExtension.testOutputDir}"
        SimulateReads().test(haploidArgs)

    }

    @Test
    fun testSimulateDiploid() {

    }
}