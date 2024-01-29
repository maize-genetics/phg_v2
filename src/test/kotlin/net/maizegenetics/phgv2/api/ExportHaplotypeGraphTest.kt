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

        exportMultiSampleHVCF(
            graph,
            "/Users/tmc46/git/phg_v2/testMultipleFilesHaplotypeGraph.vcf",
            TestExtension.smallseqRefFile
        )

    }

}