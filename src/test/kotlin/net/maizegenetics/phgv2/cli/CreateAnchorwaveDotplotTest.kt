package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class CreateAnchorwaveDotplotTest {
    @Test
    fun testCreateAnchorwaveDotplot() {
        // This functionality is well tested in the AlignAssembliesTest
        // Here we just test that the class can be called and that it runs without error
        // when invoked from CreateAnchorwaveDotplot

        // TODO _ lcj, finishe this!


    }

    @Test
    fun testCliktParams() {
        val createAWdp = CreateAnchorwaveDotplot()
        val inputFile = "data/test/smallseq/dummy_anchors_small.anchorspro"
        val outputFile = "${TestExtension.tempDir}/dummy_anchors_small.png"

        // Test missing output file
        val resultMissingOutput =
            createAWdp.test("--input-file ${inputFile}")
        assertEquals(resultMissingOutput.statusCode, 1)
        assertEquals(
            "Usage: create-anchorwave-dotplot [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --output-file: --output-file must not be blank\n", resultMissingOutput.output
        )

        // test missing input file
        val resultMissingInput =
            createAWdp.test("--output-file ${outputFile}")
        assertEquals(resultMissingInput.statusCode, 1)
        assertEquals(
            "Usage: create-anchorwave-dotplot [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --input-file: --input-file must not be blank\n", resultMissingInput.output
        )

    }
}