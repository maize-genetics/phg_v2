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
class MergeHvcfsTest {

    private val myLogger = LogManager.getLogger(MergeHvcfsTest::class.java)

    companion object {

        private val exportHvcfDir = "${TestExtension.tempDir}/merge-hvcfs/"
        private val singleInputDir = "${TestExtension.tempDir}/single-input/"
        private val multiInputDir = "${TestExtension.tempDir}/multi-input/"
        private val outputHvcfDir = "${exportHvcfDir}/output/"

        @BeforeAll
        @JvmStatic
        fun setup() {
            File(exportHvcfDir).mkdirs()
            File(outputHvcfDir).mkdirs()
            File(singleInputDir).mkdirs()

            File(TestExtension.smallseqRefHvcfFile).copyTo(File(singleInputDir + File(TestExtension.smallseqRefHvcfFile).name))

            File(TestExtension.smallseqRefHvcfFile).copyTo(File(multiInputDir + File(TestExtension.smallseqRefHvcfFile).name))
            File(TestExtension.smallseqLineAHvcfFile).copyTo(File(multiInputDir + File(TestExtension.smallseqLineAHvcfFile).name))
            File(TestExtension.smallseqLineBHvcfFile).copyTo(File(multiInputDir + File(TestExtension.smallseqLineBHvcfFile).name))
        }

    }

    @Test
    fun testSingleSampleMergeHvcfs() {

        val outputFile = "${outputHvcfDir}testSingleSampleHaplotypeGraph.vcf"

        // phg merge-hvcfs --input-dir <singleInputDir> --id-format CHECKSUM --reference-file TestExtension.smallseqRefFile --output-file exported-vcfs
        val result = MergeHvcfs().test(
            "--input-dir $singleInputDir --id-format CHECKSUM --reference-file ${TestExtension.smallseqRefFile} --output-file $outputFile"
        )

        myLogger.info("testSingleSampleMergeHvcfs: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        var checksum1 = getChecksum(TestExtension.exportGraphSingleSample)
        var checksum2 = getChecksum(outputFile)

        myLogger.info("testSingleSampleHaplotypeGraph.vcf expected checksum1: $checksum1")
        myLogger.info("testSingleSampleHaplotypeGraph.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testSingleSampleHaplotypeGraph.vcf checksums do not match")

    }

    @Test
    fun testMultipleFilesMergeHvcfs() {

        val outputFile = "${outputHvcfDir}testMultipleFilesHaplotypeGraph.vcf"

        val result = MergeHvcfs().test(
            "--input-dir $multiInputDir --id-format CHECKSUM --reference-file ${TestExtension.smallseqRefFile} --output-file $outputFile"
        )

        myLogger.info("testMultipleFilesMergeHvcfs: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        var checksum1 = getChecksum(TestExtension.exportGraphMultiSample)
        var checksum2 = getChecksum(outputFile)

        myLogger.info("testMultipleFilesHaplotypeGraph.vcf expected checksum1: $checksum1")
        myLogger.info("testMultipleFilesHaplotypeGraph.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testMultipleFilesHaplotypeGraph.vcf checksums do not match")

    }

    @Test
    fun testSingleSampleRoundTripMergeHvcfs() {

        val outputFile = "${outputHvcfDir}testSingleSampleRoundTripMergeHvcfs.vcf"

        val result = MergeHvcfs().test(
            "--input-dir $singleInputDir --id-format CHECKSUM --reference-file ${TestExtension.smallseqRefFile} --output-file $outputFile"
        )

        myLogger.info("testSingleSampleRoundTripMergeHvcfs: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        var checksum1 = getChecksum(TestExtension.exportGraphSingleSample)
        var checksum2 = getChecksum(outputFile)

        myLogger.info("testSingleSampleRoundTripMergeHvcfs.vcf expected checksum1: $checksum1")
        myLogger.info("testSingleSampleRoundTripMergeHvcfs.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testSingleSampleRoundTripMergeHvcfs.vcf checksums do not match")

    }

    @Test
    fun testMultipleSamplesRoundTripMergeHvcfs() {

        val outputFile = "${outputHvcfDir}testMultipleSamplesRoundTripMergeHvcfs.vcf"

        val result = MergeHvcfs().test(
            "--input-dir $multiInputDir --id-format CHECKSUM --reference-file ${TestExtension.smallseqRefFile} --output-file $outputFile"
        )

        myLogger.info("testMultipleSamplesRoundTripMergeHvcfs: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        var checksum1 = getChecksum(TestExtension.exportGraphMultiSample)
        var checksum2 = getChecksum(outputFile)

        myLogger.info("testMultipleSamplesRoundTripMergeHvcfs.vcf expected checksum1: $checksum1")
        myLogger.info("testMultipleSamplesRoundTripMergeHvcfs.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testMultipleSamplesRoundTripMergeHvcfs.vcf checksums do not match")

    }

    @Test
    fun testSingleSampleMergeHvcfsRangeSampleGamete() {

        val outputFile = "${outputHvcfDir}testSingleSampleMergeHvcfsRangeSampleGamete.vcf"

        val result = MergeHvcfs().test(
            "--input-dir $singleInputDir --id-format RANGE_SAMPLE_GAMETE --reference-file ${TestExtension.smallseqRefFile} --output-file $outputFile"
        )

        myLogger.info("testSingleSampleMergeHvcfsRangeSampleGamete: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        var checksum1 = getChecksum(TestExtension.exportGraphSingleSampleRangeSampleGamete)
        var checksum2 = getChecksum(outputFile)

        myLogger.info("testSingleSampleMergeHvcfsRangeSampleGamete.vcf expected checksum1: $checksum1")
        myLogger.info("testSingleSampleMergeHvcfsRangeSampleGamete.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testSingleSampleMergeHvcfsRangeSampleGamete.vcf checksums do not match")

    }

}