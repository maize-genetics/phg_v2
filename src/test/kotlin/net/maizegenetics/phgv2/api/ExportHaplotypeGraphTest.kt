package net.maizegenetics.phgv2.api

import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.getChecksum
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class ExportHaplotypeGraphTest {

    private val myLogger = LogManager.getLogger(ExportHaplotypeGraphTest::class.java)

    companion object {

        private val exportHvcfDir = "${TestExtension.tempDir}/export-vcfs/"
        private val outputHvcfDir = "${exportHvcfDir}/output/"

        @BeforeAll
        @JvmStatic
        fun setup() {
            File(exportHvcfDir).mkdirs()
            File(outputHvcfDir).mkdirs()
        }

    }

    @Test
    fun testSingleSampleHaplotypeGraph() {

        val graph = HaplotypeGraph(listOf(TestExtension.smallseqRefHvcfFile))

        val outputFile = "${outputHvcfDir}testSingleSampleHaplotypeGraph.vcf"

        exportMultiSampleHVCF(graph, outputFile, TestExtension.smallseqRefFile)

        var checksum1 = getChecksum(TestExtension.exportGraphSingleSample)
        var checksum2 = getChecksum(outputFile)

        myLogger.info("testSingleSampleHaplotypeGraph.vcf expected checksum1: $checksum1")
        myLogger.info("testSingleSampleHaplotypeGraph.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testSingleSampleHaplotypeGraph.vcf checksums do not match")

    }

    @Test
    fun testMultipleFilesHaplotypeGraph() {

        val graph = HaplotypeGraph(
            listOf(
                TestExtension.smallseqLineAHvcfFile,
                TestExtension.smallseqLineBHvcfFile,
                TestExtension.smallseqRefHvcfFile
            )
        )

        val outputFile = "${outputHvcfDir}testMultipleFilesHaplotypeGraph.vcf"

        exportMultiSampleHVCF(graph, outputFile, TestExtension.smallseqRefFile)

        var checksum1 = getChecksum(TestExtension.exportGraphMultiSample)
        var checksum2 = getChecksum(outputFile)

        myLogger.info("testMultipleFilesHaplotypeGraph.vcf expected checksum1: $checksum1")
        myLogger.info("testMultipleFilesHaplotypeGraph.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testMultipleFilesHaplotypeGraph.vcf checksums do not match")

    }

    @Test
    fun testMultipleFilesHaplotypeGraphRangeSampleGamete() {

        val graph = HaplotypeGraph(
            listOf(
                TestExtension.smallseqLineAHvcfFile,
                TestExtension.smallseqLineBHvcfFile,
                TestExtension.smallseqRefHvcfFile
            )
        )

        val outputFile = "${outputHvcfDir}testMultipleFilesHaplotypeGraphRangeSampleGamete.vcf"

        exportMultiSampleHVCF(graph, outputFile, TestExtension.smallseqRefFile, SymbolicAllele.RANGE_SAMPLE_GAMETE)

        var checksum1 = getChecksum(TestExtension.exportGraphMultiSampleRangeSampleGamete)
        var checksum2 = getChecksum(outputFile)

        myLogger.info("testMultipleFilesHaplotypeGraphRangeSampleGamete.vcf expected checksum1: $checksum1")
        myLogger.info("testMultipleFilesHaplotypeGraphRangeSampleGamete.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testMultipleFilesHaplotypeGraphRangeSampleGamete.vcf checksums do not match")

    }

    @Test
    fun testRoundTripHaplotypeGraphRangeSampleGamete() {

        val graph = HaplotypeGraph(
            listOf(
                TestExtension.exportGraphMultiSampleRangeSampleGamete
            )
        )

        val outputFile = "${outputHvcfDir}testRoundTripHaplotypeGraphRangeSampleGamete.vcf"

        exportMultiSampleHVCF(graph, outputFile, TestExtension.smallseqRefFile, SymbolicAllele.CHECKSUM)

        var checksum1 = getChecksum(TestExtension.exportGraphMultiSample)
        var checksum2 = getChecksum(outputFile)

        myLogger.info("testRoundTripHaplotypeGraphRangeSampleGamete.vcf expected checksum1: $checksum1")
        myLogger.info("testRoundTripHaplotypeGraphRangeSampleGamete.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testRoundTripHaplotypeGraphRangeSampleGamete.vcf checksums do not match")

    }

    @Test
    fun testSingleSampleRoundTripHaplotypeGraph() {

        val graph = HaplotypeGraph(listOf(TestExtension.smallseqRefHvcfFile))

        val outputFile = "${outputHvcfDir}testSingleSampleRoundTripHaplotypeGraph.vcf"

        exportMultiSampleHVCF(graph, outputFile, TestExtension.smallseqRefFile)

        var checksum1 = getChecksum(TestExtension.exportGraphSingleSample)
        var checksum2 = getChecksum(outputFile)

        myLogger.info("testSingleSampleRoundTripHaplotypeGraph.vcf expected checksum1: $checksum1")
        myLogger.info("testSingleSampleRoundTripHaplotypeGraph.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testSingleSampleRoundTripHaplotypeGraph.vcf checksums do not match")

    }

    @Test
    fun testMultipleSamplesRoundTripHaplotypeGraph() {

        val graph = HaplotypeGraph(listOf(TestExtension.exportGraphMultiSample))

        val outputFile = "${outputHvcfDir}testMultipleSamplesRoundTripHaplotypeGraph.vcf"

        exportMultiSampleHVCF(graph, outputFile, TestExtension.smallseqRefFile)

        var checksum1 = getChecksum(TestExtension.exportGraphMultiSample)
        var checksum2 = getChecksum(outputFile)

        myLogger.info("testMultipleSamplesRoundTripHaplotypeGraph.vcf expected checksum1: $checksum1")
        myLogger.info("testMultipleSamplesRoundTripHaplotypeGraph.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testMultipleSamplesRoundTripHaplotypeGraph.vcf checksums do not match")

    }

    @Test
    fun testSingleSampleHaplotypeGraphRangeSampleGamete() {

        val graph = HaplotypeGraph(listOf(TestExtension.smallseqRefHvcfFile))

        val outputFile = "${outputHvcfDir}testSingleSampleHaplotypeGraphRangeSampleGamete.vcf"

        exportMultiSampleHVCF(graph, outputFile, TestExtension.smallseqRefFile, SymbolicAllele.RANGE_SAMPLE_GAMETE)

        var checksum1 = getChecksum(TestExtension.exportGraphSingleSampleRangeSampleGamete)
        var checksum2 = getChecksum(outputFile)

        myLogger.info("testSingleSampleHaplotypeGraphRangeSampleGamete.vcf expected checksum1: $checksum1")
        myLogger.info("testSingleSampleHaplotypeGraphRangeSampleGamete.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testSingleSampleHaplotypeGraphRangeSampleGamete.vcf checksums do not match")

    }

    @Test
    fun testSingleSampleRoundTripHaplotypeGraphWithBedfile() {

        val graph = HaplotypeGraph(listOf(TestExtension.smallseqRefHvcfFile))

        val outputFile = "${outputHvcfDir}testSingleSampleRoundTripHaplotypeGraphWithBedfile.vcf"

        exportMultiSampleHVCF(
            graph,
            outputFile,
            TestExtension.smallseqRefFile,
            SymbolicAllele.RANGE_SAMPLE_GAMETE,
            rangeBedfile = TestExtension.smallseqAnchorsBedFile
        )

        var checksum1 = getChecksum(TestExtension.exportGraphSingleSampleRangeSampleGamete)
        var checksum2 = getChecksum(outputFile)

        myLogger.info("testSingleSampleRoundTripHaplotypeGraphWithBedfile.vcf expected checksum1: $checksum1")
        myLogger.info("testSingleSampleRoundTripHaplotypeGraphWithBedfile.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testSingleSampleRoundTripHaplotypeGraph.vcf checksums do not match")

    }

}