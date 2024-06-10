package net.maizegenetics.phgv2.api

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.MergeGVCFs
import net.maizegenetics.phgv2.cli.TestExtension
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith

@ExtendWith(TestExtension::class)
class MergeGVCFsTest {

    private val myLogger = LogManager.getLogger(MergeGVCFsTest::class.java)

    @Test
    fun testMergeGVCFs() {

        try {
            MergeGVCFs().test("--input-dir inputDir --output-file outputFile.vcf")
        } catch (e: Exception) {
            require(e.message == "Input GVCF directory does not exist: inputDir")
        }

        var result = MergeGVCFs().test("--output-file outputFile.vcf")
        require(result.output.contains("Error: missing option --input-dir"))

        result = MergeGVCFs().test("--input-dir inputDir")
        require(result.output.contains("Error: missing option --output-file"))

    }

}