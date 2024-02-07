package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.utils.bgzipAndIndexGVCFfile
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.assertEquals

class HvcfFromPHGGvcfTest {
    companion object {

        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.testVCFDir).mkdirs()
            File(TestExtension.testTileDBURI).mkdirs()
        }

        @JvmStatic
        @AfterAll
        fun tearDown() {
            File(TestExtension.testVCFDir).deleteRecursively()
            File(TestExtension.testTileDBURI).deleteRecursively()
        }
    }

    @Test
    fun testCliktParams() {
        val hvcfFromPhgGvcf = HvcfFromPhgGvcf()

        // There are only 3 required parameters - test for missing each one
        val resultMissingBed =
            hvcfFromPhgGvcf.test("--db-path ${TestExtension.testTileDBURI} --gvcf-dir ${TestExtension.testVCFDir} --reference-file ${TestExtension.testRefFasta} ")
        assertEquals(resultMissingBed.statusCode, 1)
        assertEquals(
            "Usage: hvcf-from-phg-gvcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --bed: --bed must not be blank\n", resultMissingBed.output
        )
        val resultMissingRef =
            hvcfFromPhgGvcf.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} --gvcf-dir ${TestExtension.testMafDir}")
        assertEquals(resultMissingRef.statusCode, 1)
        assertEquals(
            "Usage: hvcf-from-phg-gvcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --reference-file: --reference-file must not be blank\n", resultMissingRef.output
        )


        val resultMissingGvcfDir =
            hvcfFromPhgGvcf.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} --reference-file ${TestExtension.testRefFasta}")
        assertEquals(resultMissingGvcfDir.statusCode, 1)
        assertEquals(
            "Usage: hvcf-from-phg-gvcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --gvcf-dir: --gvcf-dir must not be blank\n", resultMissingGvcfDir.output
        )

    }

    @Test
    fun testSimpleHvcfFromGvcf() {
        // Copy the gvcf files from data/test/smallseq to the testVCFDir
        val gvcfDir = TestExtension.testVCFDir
        val gvcfFiles = File("data/test/smallseq").listFiles { _, name -> name.endsWith(".g.vcf") }
        gvcfFiles.forEach { file -> file.copyTo(File(gvcfDir, file.name)) }

        // HvcfFromPhgGvcf can use compressed or uncompressed file.
        // Maybe have some of each in here.
        // Run the bgzipAndIndexGVCFfile function on each file
        val gvcfFilesList = File(gvcfDir).listFiles { _, name -> name.endsWith(".g.vcf") }
        gvcfFilesList.forEach { file -> bgzipAndIndexGVCFfile(file.toString()) }

        //Need to create the agc record before we run this:
        // We will not be loading to tiledb, but we will be pulling sequence from AGC.
        // AGC lives in the same directory as the tiledb datasets.  Need all the fasta files in the agc record
        // for each sample in the gvcf files

        // copy the fasta files
        val fastaInputDir = "data/test/smallseq"
        val fastaOutputDir = TestExtension.testOutputFastaDir
        val fastaFiles = File(fastaInputDir).listFiles { file -> file.extension == "fa" }
        fastaFiles.forEach { file -> file.copyTo(File(fastaOutputDir, file.name) ) }

        val dbPath = TestExtension.testTileDBURI
        val refFasta = TestExtension.smallseqRefFile
        println("refFasta: $refFasta")

        // get the full path fasta file names  from the fastaInput, write to fileList
        val fileList = mutableListOf<String>()
        File(fastaOutputDir).walk(FileWalkDirection.TOP_DOWN).filter{it.name.endsWith(".fa")}.forEach {
            fileList.add(it.toString())
        }
        // write the full path fasta file names  from the fastaOutputDir to a single file, one per line, in tempDir
        // This file will be used as input to the agc-compress command
        val fastaCreateFileNamesFile = File(dbPath, "agcFileList.txt")
        fastaCreateFileNamesFile.writeText(fileList.joinToString("\n"))

        val agcCompress = AgcCompress()
        // Create the compressed file
        val agcResult =
            agcCompress.test("--fasta-list ${fastaCreateFileNamesFile} --db-path ${dbPath} --reference-file ${refFasta}")
        println(agcResult.output)

        val bedFile = TestExtension.smallseqAnchorsBedFile
        val hvcfFromGvcf = HvcfFromPhgGvcf()
        val result =
            hvcfFromGvcf.test("--db-path ${dbPath} --bed ${bedFile} --reference-file ${refFasta} --gvcf-dir ${gvcfDir} ")
        println(result.output)

        // TODO verify the h.vcf is correct! Should I run CreatMafVcf and compare the outputs?
        // FIX the name of the file !! (LineA.g.vcf.h.vcf.gz is not a good name)

    }

    @Test
    fun compareToCreateMafVcfOutput() {
        // This test will compare the output of HvcfFromPhgGvcf to the hvcf output of CreateMafVcf
        // the gvcf used to create the hvcf via the new HvcfFromPhgGvcf code is the gvcf
        // created by CreateMafVcf.

        // Setup and run CreateMafVcf to get the gvcf file for testing
        val fastaCreateFileNamesFile = "data/test/buildMAFVCF/fastaCreateFileNames.txt"
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/buildMAFVCF/B73_Test.fa"


        val agcCompress = AgcCompress()
        // Create the initial compressed file
        val agcResult = agcCompress.test("--fasta-list ${fastaCreateFileNamesFile} --db-path ${dbPath} --reference-file ${refFasta}")
        println(agcResult.output)

        val createMAFVCF = CreateMafVcf()
        val result = createMAFVCF.test("--db-path ${dbPath} --bed data/test/buildMAFVCF/B73_Test.bed --reference-file ${refFasta} --maf-dir data/test/buildMAFVCF/mafs/ -o ${TestExtension.testVCFDir}")
        println(result.output)


    }
}