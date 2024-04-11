package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
class CreateAnchorwaveDotplotTest {
    @Test
    fun testCreateAnchorwaveDotplot() {
        // This functionality is well tested in the AlignAssembliesTest
        // Here we just test that the class can be called and that it runs without error
        // when invoked from CreateAnchorwaveDotplot

        val createAWdp = CreateAnchorwaveDotplot()
        val inputFile = "data/test/smallseq/dummy_anchors_small.anchorspro"
        val outputFile = "${TestExtension.tempDir}/dummy_anchors_small.svg"
        val result = createAWdp.test("--input-file ${inputFile} --output-file ${outputFile}")
        assertEquals(result.statusCode, 0)
        assertTrue(File(outputFile).exists(),"Output file not created")
    }

    @Test
    fun testDotPlotWithRelativePath() {
        val createAWdp = CreateAnchorwaveDotplot()
        val inputFile = "data/test/smallseq/dummy_anchors_small.anchorspro"
        val outputFile = "dummy_anchors_small.svg"

        val outputFilePath = "${System.getProperty("user.dir")}/$outputFile"
        val result = createAWdp.test("--input-file ${inputFile} --output-file ${outputFile}")
        assertEquals(result.statusCode, 0)
        assertTrue(File(outputFilePath).exists(),"Output file not created")

        // Delete the outputFilePath
        File(outputFilePath).delete()
    }
    @Test
    fun testCliktParams() {
        val createAWdp = CreateAnchorwaveDotplot()
        val inputFile = "data/test/smallseq/dummy_anchors_small.anchorspro"
        val outputFile = "${TestExtension.tempDir}/dummy_anchors_small.svg"

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