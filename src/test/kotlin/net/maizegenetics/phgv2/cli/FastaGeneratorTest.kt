package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class FastaGeneratorTest {
    companion object {

    }

    @Test
    fun testCliktParams() {
        val fastaGenerator = FastaGenerator()
        val resultGood = fastaGenerator.test("--db-path ${TestExtension.testTileDBURI} --agc-file ${TestExtension.testAGCFile} --sample-name Ref -o ${TestExtension.testOutputRefFasta}")

        val resultMissingDB = fastaGenerator.test("--agc-file ${TestExtension.testAGCFile} --sample-name Ref -o ${TestExtension.testOutputRefFasta}")
        assertEquals(resultMissingDB.statusCode, 1)
        assertEquals("Usage: fasta-generator [<options>]\n" +
                "\n" +
                "Error: invalid value for --db-path: --db-path must not be blank\n",resultMissingDB.output)
        val resultMissingAGC = fastaGenerator.test("--db-path ${TestExtension.testTileDBURI} --sample-name Ref -o ${TestExtension.testOutputRefFasta}")
        assertEquals(resultMissingAGC.statusCode, 1)
        assertEquals("Usage: fasta-generator [<options>]\n" +
                "\n" +
                "Error: invalid value for --agc-file: --agc-file must not be blank\n",resultMissingAGC.output)

        val resultMissingSampleName = fastaGenerator.test("--db-path ${TestExtension.testTileDBURI} --agc-file ${TestExtension.testAGCFile} -o ${TestExtension.testOutputRefFasta}")
        assertEquals(resultMissingSampleName.statusCode, 1)
        assertEquals("Usage: fasta-generator [<options>]\n" +
                "\n" +
                "Error: invalid value for --sample-name: --sample-name must not be blank\n",resultMissingSampleName.output)

        val resultMissingOutput = fastaGenerator.test("--db-path ${TestExtension.testTileDBURI} --agc-file ${TestExtension.testAGCFile} --sample-name Ref")
        assertEquals(resultMissingOutput.statusCode, 1)
        assertEquals("Usage: fasta-generator [<options>]\n" +
                "\n" +
                "Error: invalid value for --output: --output/-o must not be blank\n",resultMissingOutput.output)
    }
}