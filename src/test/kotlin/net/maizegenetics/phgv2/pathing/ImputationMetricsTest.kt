package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File

/**
 * THe test files for this class were created from both the FullPipelineIT test and the
 * FindPaths Test..  From FullPipelineTest, the LineB.h.vcf file created from the
 * create-maf-vcf call was used as the parentSampleHvcf file.  The imputationHvcf file was
 * created from the FIndPathsTest:testHaploidPathFinding() test.  In this test, I temporaritly
 * altered the input for the
 *   val switchResult = FindPaths().test(switchTestArgs)
 *
 * The above used the hvcf created from FullPipelineIT instead of the hvcf stored in the repository.
 * This was to ensure the h.vcf had ALT headers with the RefRange values set to contig:start-end
 * vs a haplotype check sum.  This format is needed to get valid metrics.
 */
@ExtendWith(TestExtension::class)
class ImputationMetricsTest {
    companion object {
        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.testOutputDir).deleteRecursively()
            File(TestExtension.testOutputDir).mkdirs()

        }

        @JvmStatic
        @AfterAll
        fun tearDown() {
            //comment out the following line to inspect the test results after the tests have been run
            File(TestExtension.testOutputDir).deleteRecursively()
        }
    }


    @Test
    fun testImputationMetrics() {


        val imputationHvcf = "data/test/resequenceHaplotypeVCF/Imputation.h.vcf"

        val parentSampleHvcf = "data/test/resequenceHaplotypeVCF/LineB.h.vcf"
        val readMappingFiles = "data/test/resequenceHaplotypeVCF/readMappingFileList.txt"

        var outputFile = "${TestExtension.testOutputDir}/smallSeq_verificationScriptResults_chr1.txt"
        var chrom = "1"
        val bedFile = "data/test/smallseq/anchors.bed"

        // Test chr1
        val imputationMetricsResultsChr1 = ImputationMetrics().test("--imputation-hvcf $imputationHvcf --sample-hvcf $parentSampleHvcf --read-mapping-files $readMappingFiles " +
                "--output-file $outputFile --chrom ${chrom} --bed-file $bedFile")
        Assertions.assertEquals(
            0,
            imputationMetricsResultsChr1.statusCode,
            "verifyImputationResults status code was ${imputationMetricsResultsChr1.statusCode}"
        )

        var resultsFileLines = bufferedReader(outputFile).readLines()
        var headerLine = resultsFileLines[0]
        var headerColumns = headerLine.split("\t")
        // Verify the number of lines is 21
        Assertions.assertEquals(21, resultsFileLines.size)

        // for lines in the resultsFileLines, skip the header line, then verify all other lines start with "1"
        resultsFileLines = resultsFileLines.drop(1)

        for (line in resultsFileLines) {
            val cols = line.split("\t")
            val chr = cols[0].split(":")[0]
            Assertions.assertEquals("1", chr)
            val start = cols[0].split(":")[1].split("-")[0].toInt()
            val imputationSampleIdx = headerColumns.indexOf("ImputationSampleName")

            // This was setup so positions < 25000 are LineA and positions >= 25000 are LineB on either chrom 1 or 2
            // There were no reads mapping to the last range, so the sample name there is NA
            if (start <= 25000) assert(cols[imputationSampleIdx] == "LineA") {"Imputation sample is not LineA"}
            else if (start < 50501) assert(cols[imputationSampleIdx] == "LineB") {"Imputation sample is not LineB"}
        }

        // run again with only chrom2
        chrom = "2"
        outputFile = "${TestExtension.testOutputDir}/smallSeq_verificationScriptResults_chr2.txt"
        val imputationMetricsResultsChr2 = ImputationMetrics().test("--imputation-hvcf $imputationHvcf --sample-hvcf $parentSampleHvcf --read-mapping-files $readMappingFiles " +
                "--output-file $outputFile --chrom ${chrom} --bed-file $bedFile")
        Assertions.assertEquals(
            0,
            imputationMetricsResultsChr2.statusCode,
            "verifyImputationResults status code was ${imputationMetricsResultsChr2.statusCode}"
        )
        resultsFileLines = bufferedReader(outputFile).readLines()
        headerLine = resultsFileLines[0]
        headerColumns = headerLine.split("\t")
        // Verify the number of lines is 21
        Assertions.assertEquals(21, resultsFileLines.size)

        // for lines in the resultsFileLines, skip the header line, then verify all other lines start with "2"
        resultsFileLines = resultsFileLines.drop(1)

        for (line in resultsFileLines) {
            val cols = line.split("\t")
            val chr = cols[0].split(":")[0]
            Assertions.assertEquals("2", chr)
            val imputationSampleIdx = headerColumns.indexOf("ImputationSampleName")
            val start = cols[0].split(":")[1].split("-")[0].toInt()
            // This was setup so positions < 25000 are LineA and positions >= 25000 are LineB on both chrom 1 and 2
            // There were no reads mapping to the last range, so the sample name there is NA
            if (start <= 25000) assert(cols[imputationSampleIdx] == "LineA") {"Imputation sample is not LineA"}
            else if (start < 50501) assert(cols[imputationSampleIdx] == "LineB") {"Imputation sample is not LineB"}
        }

        // run again with all chroms
        outputFile = "${TestExtension.testOutputDir}/smallSeq_verificationScriptResults_allChroms.txt"
        val imputationMetricsResultsAllChroms = ImputationMetrics().test("--imputation-hvcf $imputationHvcf --sample-hvcf $parentSampleHvcf --read-mapping-files $readMappingFiles " +
                "--output-file $outputFile  --bed-file $bedFile")
        Assertions.assertEquals(
            0,
            imputationMetricsResultsAllChroms.statusCode,
            "verifyImputationResults status code was ${imputationMetricsResultsAllChroms.statusCode}"
        )
        // verify both chr1 and chr2 show up in the file, total lines shouldbe 41
        resultsFileLines = bufferedReader(outputFile).readLines()
        headerLine = resultsFileLines[0]
        headerColumns = headerLine.split("\t")
        Assertions.assertEquals(41, resultsFileLines.size)

        // for lines in the resultsFileLines, skip the header line, process the rest
        resultsFileLines = resultsFileLines.drop(1)
        var chr1Count = 0
        var chr2Count = 0
        for (line in resultsFileLines) {
            if (line.startsWith("1")) chr1Count++
            if (line.startsWith("2")) chr2Count++

            val cols = line.split("\t")
            val imputationSampleIdx = headerColumns.indexOf("ImputationSampleName")
            val start = cols[0].split(":")[1].split("-")[0].toInt()
            // This was setup so positions < 25000 are LineA and positions >= 25000 are LineB for both chr1 and chr2
            // There were no reads mapping to the last range, so the sample name there is NA
            if (start <= 25000) assert(cols[imputationSampleIdx] == "LineA") {"Imputation sample is not LineA"}
            else if (start < 50501) assert(cols[imputationSampleIdx] == "LineB") {"Imputation sample is not LineB"}
        }
        // verify chr1Count = 20 and chr2Count = 20
        Assertions.assertEquals(20, chr1Count)
        Assertions.assertEquals(20, chr2Count)

    }
}