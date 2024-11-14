package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.utils.getChecksum
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class SampleHapidByRangeTest {

    private val myLogger = LogManager.getLogger(SampleHapidByRangeTest::class.java)

    companion object {

        private val exportDir = "${TestExtension.tempDir}/sample-hapid-by-range/"
        private val multiInputDir = "${TestExtension.tempDir}/sample-hapid-by-range-input/"
        private val expectOutputFile = "data/test/sample-hapid-by-range-output.txt"

        @BeforeAll
        @JvmStatic
        fun setup() {
            File(exportDir).mkdirs()
            File(multiInputDir).mkdirs()

            File(TestExtension.smallseqRefHvcfFile).copyTo(File(multiInputDir + File(TestExtension.smallseqRefHvcfFile).name))
            File(TestExtension.smallseqLineAHvcfFile).copyTo(File(multiInputDir + File(TestExtension.smallseqLineAHvcfFile).name))
            File(TestExtension.smallseqLineBHvcfFile).copyTo(File(multiInputDir + File(TestExtension.smallseqLineBHvcfFile).name))
        }

    }

    @Test
    fun testSampleHapidByRangeTest() {

        val outputFile = "${exportDir}testSampleHapidByRangeTest.txt"

        // phg sample-hapid-by-range --input-dir hvcf-files --output-file sample-hapid-by-range.txt
        val result = SampleHapidByRange().test(
            "--input-dir $multiInputDir --output-file $outputFile"
        )

        myLogger.info("testSampleHapidByRangeTest: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        var checksum1 = getChecksum(expectOutputFile)
        var checksum2 = getChecksum(outputFile)

        myLogger.info("testSampleHapidByRangeTest expected checksum1: $checksum1")
        myLogger.info("testSampleHapidByRangeTest actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testSampleHapidByRangeTest.vcf checksums do not match")

    }

}