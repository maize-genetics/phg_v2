package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class BuildMAFVCFTest {
    companion object {

    }

    @Test
    fun testCliktParams() {
        val buildMAFVCF = BuildMafVcf()
        val resultGood =
            buildMAFVCF.test("--bed ${TestExtension.testBEDFile} --reference ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")

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
}