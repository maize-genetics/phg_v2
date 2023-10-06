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
        val resultGood = buildRefVCF.test("--bed ${TestExtension.testBEDFile} --reference ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")

        val resultMissingBed = buildRefVCF.test("--reference ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingBed.statusCode, 1)
        assertEquals("Usage: build-ref-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --bed: --bed must not be blank\n",resultMissingBed.output)
        val resultMissingRef = buildRefVCF.test("--bed ${TestExtension.testBEDFile} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingRef.statusCode, 1)
        assertEquals("Usage: build-ref-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --reference: --reference must not be blank\n",resultMissingRef.output)


        val resultMissingOutput = buildRefVCF.test("--bed ${TestExtension.testBEDFile} --reference ${TestExtension.testRefFasta}")
        assertEquals(resultMissingOutput.statusCode, 1)
        assertEquals("Usage: build-ref-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --output-dir: --output-dir/-o must not be blank\n",resultMissingOutput.output)
    }
}