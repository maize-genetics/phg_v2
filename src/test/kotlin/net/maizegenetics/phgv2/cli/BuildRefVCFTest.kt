package net.maizegenetics.phgv2.cli

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class BuildRefVCFTest {
    companion object {
        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"

        @JvmStatic
        @BeforeAll
        fun setup() {
            File(tempDir).mkdirs()
        }

        // Comment these out if you need to look at the logs files created by ProcessBuilder
        // commands.
        @JvmStatic
        @AfterAll
        fun teardown() {
            File(tempDir).deleteRecursively()
        }
    }

    @Test
    fun testCliktParams() {
        val buildRefVCF = BuildRefVcf()
        //val resultGood = buildRefVCF.test("--bed ${TestExtension.testBEDFile} --ref-fasta ${TestExtension.testRefFasta} --ref-name ${TestExtension.refLineName} -o ${TestExtension.testVCFDir}")

        // Send request without bed file parameter
        val resultMissingBed = buildRefVCF.test("--ref-name ${TestExtension.refLineName} --ref-fasta ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingBed.statusCode, 1)
        assertEquals("Usage: build-ref-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --bed: --bed must not be blank\n",resultMissingBed.output)

        // Send request without reference file parameter
        val resultMissingRef = buildRefVCF.test("--ref-name ${TestExtension.refLineName} --bed ${TestExtension.testBEDFile} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingRef.statusCode, 1)
        assertEquals("Usage: build-ref-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --ref-fasta: --ref-fasta must not be blank\n",resultMissingRef.output)

        // Send request without reference name parameter
        val resultMissingRefName = buildRefVCF.test("--ref-fasta ${TestExtension.testRefFasta} --bed ${TestExtension.testBEDFile} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingRef.statusCode, 1)
        assertEquals("Usage: build-ref-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --ref-name: --ref-name must not be blank\n",resultMissingRef.output)


        // Send request without output directory parameter
        val resultMissingOutput = buildRefVCF.test("--bed ${TestExtension.testBEDFile} --ref-name ${TestExtension.refLineName} --ref-fasta ${TestExtension.testRefFasta} ")
        assertEquals(resultMissingOutput.statusCode, 1)
        assertEquals("Usage: build-ref-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --output-dir: --output-dir/-o must not be blank\n",resultMissingOutput.output)
    }

    @Test
    fun testBuildRefVCF() {
        val vcfDir = tempDir
        val refName = "Ref"
        // TODO: Replace ranges/genome with the smallseq files when the test case branch is merged.

        // val ranges = "data/test/smallseq/anchors.bed"
        // val genome = "data/test/smallseq/Ref.fa"
        val ranges = "/Users/lcj34/notes_files/phg_v2/smallSeq_final/anchors.bed"
        val genome = "/Users/lcj34/notes_files/phg_v2/smallSeq_final/Ref.fa"

        val buildRefVCF = BuildRefVcf()
        val result = BuildRefVcf().createRefHvcf(ranges,genome,refName,vcfDir)

        val outFileCompressed = "${tempDir}Ref.hvcf.gz"
        val outFileIndexed = "${tempDir}Ref.hvcf.gz.csi"
        // Verify the outputFile exists
        assertEquals(true, File(outFileCompressed).exists())
        assertEquals(true, File(outFileIndexed).exists())


        // test cases to verify the ALT headers lines in the file
        // This is the reference HVCF, so each vcf data line is a single haplotype,
        // and these haplotypes match 1-1 with the entries in the bed file.
        // IE, there is a single hvcf entry for each bed entry.

        // Verify the number of hvcf "data" lines match the number of bed file lines
        val numberOfLinesInBed = bufferedReader(ranges).readLines().size
        val numberOfDataLinesInHvcf = bufferedReader(outFileCompressed).readLines().filter { !it.startsWith("#") }.size
        assertEquals(numberOfLinesInBed, numberOfDataLinesInHvcf)

        // Verify the number of hvcf "##ALT" lines match the number of bed file lines
        val numberOfAltLinesInHvcf = bufferedReader(outFileCompressed).readLines().filter { it.startsWith("##ALT") }.size
        assertEquals(numberOfLinesInBed, numberOfAltLinesInHvcf)

    }
}