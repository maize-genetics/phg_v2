package net.maizegenetics.phgv2.cli

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.testing.test
import com.google.common.collect.Range
import com.google.common.collect.TreeRangeSet
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.utils.Position
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter
import kotlin.test.assertEquals
import kotlin.test.fail

@ExtendWith(TestExtension::class)
class BuildMAFVCFTest {
    companion object {

    }

    @Test
    fun testCliktParams() {
        val buildMAFVCF = BuildMafVcf()

        val resultMissingBed =
            buildMAFVCF.test("--maf-dir ${TestExtension.testMafDir} --reference ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingBed.statusCode, 1)
        assertEquals(
            "Usage: build-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --bed: --bed must not be blank\n", resultMissingBed.output
        )
        val resultMissingRef =
            buildMAFVCF.test("--bed ${TestExtension.testBEDFile} --maf-dir ${TestExtension.testMafDir} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingRef.statusCode, 1)
        assertEquals(
            "Usage: build-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --reference: --reference must not be blank\n", resultMissingRef.output
        )

        val resultMissingOutput =
            buildMAFVCF.test("--bed ${TestExtension.testBEDFile} --maf-dir ${TestExtension.testMafDir} --reference ${TestExtension.testRefFasta}")
        assertEquals(resultMissingOutput.statusCode, 1)
        assertEquals(
            "Usage: build-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --output-dir: --output-dir/-o must not be blank\n", resultMissingOutput.output
        )

        val resultMissingMafDir =
            buildMAFVCF.test("--bed ${TestExtension.testBEDFile} --reference ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingMafDir.statusCode, 1)
        assertEquals(
            "Usage: build-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --maf-dir: --maf-dir must not be blank\n", resultMissingMafDir.output
        )
    }

    @Test
    fun testSimpleBuildMAFVCF() {
        val buildMAFVCF = BuildMafVcf()
        buildMAFVCF.test("--bed data/test/buildMAFVCF/B73_Test.bed --reference data/test/buildMAFVCF/B73_Test.fa --maf-dir data/test/buildMAFVCF/mafs/ -o ${TestExtension.testVCFDir}")

        //compare the contents of the output gVCF files to the expected output
        compareTwoGVCFFiles("data/test/buildMAFVCF/truthGVCFs/B97_truth.g.vcf", "${TestExtension.testVCFDir}/B97.g.vcf")
    }

    fun compareTwoGVCFFiles(truthGVCFFile : String, generatedFile: String) {
        //Load in the output GVCF  and the truth GVCF and verify that the output is correct
        val truthVariantIterator = VCFFileReader(File(truthGVCFFile),false).iterator()
        val truthVariants = mutableListOf<VariantContext>()
        while(truthVariantIterator.hasNext()) {
            truthVariants.add(truthVariantIterator.next())
        }
        val truthMap = truthVariants.associateBy { Position(it.contig, it.start) }

        val outputVariantIterator = VCFFileReader(File(generatedFile), false).iterator()
        val outputVariants = mutableListOf<VariantContext>()
        while(outputVariantIterator.hasNext()) {
            outputVariants.add(outputVariantIterator.next())
        }

        assertEquals(truthVariants.size, outputVariants.size,"Number of Variants does not match:")

        for(variant in outputVariants) {
            if(!truthMap.containsKey(Position(variant.contig, variant.start))) {
                fail("No matching variant found: ${variant.contig}:${variant.start}")
            }
            val matchingTruth = truthMap[Position(variant.contig, variant.start)]!!

            //Check END
            assertEquals(matchingTruth.end, variant.end,
                "End position does not match: outputFile: ${variant.contig}:${variant.start}-${variant.end}, " +
                        "Truth: ${matchingTruth.contig}:${matchingTruth.start}-${matchingTruth.end}")

            //Check alleles
            for(element in matchingTruth.alleles.zip(variant.alleles)) {
                assertEquals(element.first, element.second, "Alleles do not match: ${variant.contig}:${variant.start}-${variant.end}")
            }
            //Check GT
            assertEquals( matchingTruth.getGenotype(0).genotypeString,variant.getGenotype(0).genotypeString,
                "GT Fields do not match: ${variant.contig}:${variant.start}-${variant.end}")
            //Check AD
            for(element in matchingTruth.getGenotype(0).ad.zip(variant.getGenotype(0).ad)) {
                assertEquals(element.first, element.second, "AD fields do not match: ${variant.contig}:${variant.start}-${variant.end}")
            }
            //Check ASM Contig
            assertEquals(matchingTruth.getAttribute("ASM_Chr"), variant.getAttribute("ASM_Chr"),
                "ASM_Contig does not match: ${variant.contig}:${variant.start}-${variant.end}")
            //Check ASM Start
            assertEquals(matchingTruth.getAttribute("ASM_Start"), variant.getAttribute("ASM_Start"),
                "ASM_Start does not match: ${variant.contig}:${variant.start}-${variant.end}")
            //Check ASM END
            assertEquals(matchingTruth.getAttribute("ASM_End"), variant.getAttribute("ASM_End"),
                "ASM_End does not match: ${variant.contig}:${variant.start}-${variant.end}")
            //Check ASM Strand
            assertEquals(matchingTruth.getAttribute("ASM_Strand"), variant.getAttribute("ASM_Strand"),
                "ASM_Strand does not match: ${variant.contig}:${variant.start}-${variant.end}")

        }
    }

}