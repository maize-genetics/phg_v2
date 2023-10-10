package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.extension.ExtendWith
import kotlin.test.Test
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class LoadVCFTest {
    companion object {

    }

    @Test
    fun testCliktParams() {
        val loadVCF = LoadVcf()
        val resultGood = loadVCF.test("--vcf-dir ${TestExtension.testVCFDir} --db-path ${TestExtension.testTileDBURI}")

        // Test missing vcf-dir parameter
        val resultMissingVCFDir = loadVCF.test("--db-path ${TestExtension.testTileDBURI} --temp-dir ${TestExtension.tempDir}")
        assertEquals(resultMissingVCFDir.statusCode, 1)
        assertEquals("Usage: load-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --vcf-dir: --vcf-dir must not be blank\n",resultMissingVCFDir.output)

        // Test missing db-path parameter
        val resultMissingDB = loadVCF.test("--vcf-dir ${TestExtension.testVCFDir} --temp-dir ${TestExtension.tempDir}")
        assertEquals(resultMissingDB.statusCode, 1)
        assertEquals("Usage: load-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --db-path: --db-path must not be blank\n",resultMissingDB.output)

        // Test missing temp-dir parameter
        val resultMissingTempDir = loadVCF.test("--vcf-dir ${TestExtension.testVCFDir} --db-path ${TestExtension.testTileDBURI}")
        assertEquals(resultMissingTempDir.statusCode, 1)
        assertEquals("Usage: load-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --temp-dir: --temp-dir must not be blank\n",resultMissingTempDir.output)
    }
}