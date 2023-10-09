package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class BuildRefVCFTest {
    companion object {

    }

    @Test
    fun testCliktParams() {
        val buildRefVCF = BuildRefVcf()
        val resultGood = buildRefVCF.test("--bed ${TestExtension.testBEDFile} --refFasta ${TestExtension.testRefFasta} --refName ${TestExtension.refLineName} -o ${TestExtension.testVCFDir}")

        // Send request without bed file parameter
        val resultMissingBed = buildRefVCF.test("--reference ${TestExtension.testRefFasta} --refFasta ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingBed.statusCode, 1)
        assertEquals("Usage: build-ref-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --bed: --bed must not be blank\n",resultMissingBed.output)

        // Send request without reference file parameter
        val resultMissingRef = buildRefVCF.test("--refFasta ${TestExtension.testRefFasta} --bed ${TestExtension.testBEDFile} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingRef.statusCode, 1)
        assertEquals("Usage: build-ref-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --reference: --reference must not be blank\n",resultMissingRef.output)


        // Send request without output directory parameter
        val resultMissingOutput = buildRefVCF.test("--bed ${TestExtension.testBEDFile} --refFasta ${TestExtension.testRefFasta} --refFasta ${TestExtension.testRefFasta}")
        assertEquals(resultMissingOutput.statusCode, 1)
        assertEquals("Usage: build-ref-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --output-dir: --output-dir/-o must not be blank\n",resultMissingOutput.output)
    }

    @Test
    fun testBuildRefVCF() {
        val vcfDir = "${TestExtension.testVCFDir}"
        val buildRefVCF = BuildRefVcf()
    }
}