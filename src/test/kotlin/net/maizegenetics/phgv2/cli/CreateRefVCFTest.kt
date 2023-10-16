package net.maizegenetics.phgv2.cli

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class CreateRefVCFTest {
    companion object {
        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"

        @JvmStatic
        @BeforeAll
        fun setup() {
            File(tempDir).mkdirs()
        }

        // Comment out the tearDown()if you need to look at the logs files created by ProcessBuilder
        // commands.
        @JvmStatic
        @AfterAll
        fun teardown() {
            File(tempDir).deleteRecursively()
        }
    }

    @Test
    fun testCliktParams() {
        val createRefVCF = CreateRefVCF()

        // Testing the "good" case happens in the actual test case below

        // Test missing bed file parameter, also missing refurl
        // it is the missing bed file parameter that will be flagged
        val resultMissingBed = createRefVCF.test("--refname ${TestExtension.refLineName} --referencefile ${TestExtension.testRefFasta}  -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingBed.statusCode, 1)
        assertEquals("Usage: create-ref-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --bed: --bed must not be blank\n",resultMissingBed.output)

        // Test missing reference file
        val resultMissingRef = createRefVCF.test("--bed ${TestExtension.testBEDFile} --refurl ${TestExtension.refURL} --refname ${TestExtension.refLineName} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingRef.statusCode, 1)
        assertEquals("Usage: create-ref-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --referencefile: --referencefile must not be blank\n",resultMissingRef.output)


        // Test missing output directory
        val resultMissingOutput = createRefVCF.test("--bed ${TestExtension.testBEDFile} --refurl ${TestExtension.refURL} --refname ${TestExtension.refLineName} --referencefile ${TestExtension.testRefFasta}")
        assertEquals(resultMissingOutput.statusCode, 1)
        assertEquals("Usage: create-ref-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --output-dir: --output-dir/-o must not be blank\n",resultMissingOutput.output)

        // Test missing ref name parameter
        val resultMissingRefName = createRefVCF.test("--referencefile ${TestExtension.testRefFasta} --refurl ${TestExtension.refURL} --bed ${TestExtension.testBEDFile} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingRefName.statusCode, 1)
        println("resultMissingRefName.output = \n${resultMissingRefName.output}")
        assertEquals("Usage: create-ref-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --refname: --refname must not be blank\n",resultMissingRefName.output)

    }

    @Test
    fun testBuildRefVCF_badIntervals() {

        val anchorFile = "${tempDir}/testAnchorFile.txt"
        File(anchorFile).bufferedWriter().use {
            // Lines 4 and 5 overlap line 3
            // Line 8 overlaps line 7 (First chr2 anchor)
            // Line 13 overlaps line 12 (but different chrom, so should not count)
            val anchorContents = """
                chr1	10575	13198
                chr1	20460	29234
                chr1	141765	145627
                chr1	143661	146302
                chr1	144363	148863
                chr1	175219	177603
                chr2	81176	84375
                chr2	82776	87024
                chr2	108577	113286
                chr3	116671	122941
                chr3	140393	145805
                chr3	159053	164568
                chr5	158053	164568
            """.trimIndent()
            it.write(anchorContents)
        }

        val vcfDir = tempDir
        val refName = "Ref"
        val refUrl = TestExtension.refURL

        val genome = "data/test/smallseq/Ref.fa"

        assertThrows<IllegalArgumentException> {
            //Check that an error is thrown when the bed file has overlapping intervals
            CreateRefVCF().createRefHvcf(anchorFile,genome,refName,refUrl,vcfDir)
        }

    }

    // Test is ignored until we can figure out how to get the test to run on the github actions server.
    // ie - conda setup is working correctly.
    //@Ignore
    @Test
    fun testBuildRefVCF() {
        val vcfDir = tempDir
        val refName = "Ref"
        val refUrl = TestExtension.refURL

        val ranges = "data/test/smallseq/anchors.bed"
        val genome = "data/test/smallseq/Ref.fa"

        val createRefVcf = CreateRefVCF()

        // This could also be called via:
        //createRefVcf.createRefHvcf(ranges,genome,refName,refUrl,vcfDir)
        val result = CreateRefVCF().test("--bed $ranges --refname $refName --referencefile $genome --refurl ${refUrl} -o $vcfDir")

        val outFileCompressed = "${tempDir}Ref.h.vcf.gz"
        val outFileIndexed = "${tempDir}Ref.h.vcf.gz.csi"
        // Verify the outputFiles exist
        assertEquals(true, File(outFileCompressed).exists())
        assertEquals(true, File(outFileIndexed).exists())

        // verify the ALT headers lines in the file
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

        // Verify the reference header was added
        val numberOfReferenceLines = bufferedReader(outFileCompressed).readLines().filter { it.startsWith("##reference") }.size
        assertEquals(1, numberOfReferenceLines)

        // Verify the RefAllele for the first haplotype for each chromosome is correct
        val genomeLines = bufferedReader(genome).readLines()
        val firstChrom = genomeLines.filter { it.startsWith(">") }[0].removePrefix(">")
        val secondChrom = genomeLines.filter { it.startsWith(">") }[1].removePrefix(">")

        // get the vcf data lines
        val vcfDataLines = bufferedReader(outFileCompressed).readLines().filter { !it.startsWith("#") }

        // verify the first dataline is the first chrom
        assertEquals(firstChrom,vcfDataLines[0].split("\t")[0])

        // verify the 21th dataline is the second chrom
        assertEquals(secondChrom,vcfDataLines[20].split("\t")[0])

        // verify there are 20 data lines for the first chrom (10 genes, 10 intergenic regions)
        assertEquals(20,vcfDataLines.filter { it.split("\t")[0] == firstChrom }.size)

        //verify there are 20 data lines for the second chrom(10 genes, 10 intergenic regions)
        assertEquals(20,vcfDataLines.filter { it.split("\t")[0] == secondChrom }.size)

        // verify the ref allele for the first data line matches the first allele in the reference fasta, e.g. the genomeLines files
        assertEquals(genomeLines[1].get(0).toString(),vcfDataLines[0].split("\t")[3])
    }
}