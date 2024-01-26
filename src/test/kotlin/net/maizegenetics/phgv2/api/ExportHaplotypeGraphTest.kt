package net.maizegenetics.phgv2.api

import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import kotlin.test.assertEquals
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
class ExportHaplotypeGraphTest {

    @Test
    fun testSingleSampleHaplotypeGraph() {

        val graph = HaplotypeGraph(listOf(TestExtension.smallseqRefHvcfFile))

        exportMultiSampleHVCF(graph, "/Users/tmc46/git/phg_v2/testHaplotypeGraph.vcf", TestExtension.smallseqRefFile)

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
            "/Users/tmc46/git/phg_v2/testHaplotypeGraphMultiSample.vcf",
            TestExtension.smallseqRefFile
        )

    }

}