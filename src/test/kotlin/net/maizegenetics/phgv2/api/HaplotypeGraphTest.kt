package net.maizegenetics.phgv2.api

import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import kotlin.test.assertEquals
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
class HaplotypeGraphTest {

    @Test
    fun testHaplotypeGraph() {

        val graph = HaplotypeGraph(listOf(TestExtension.smallseqRefHvcfFile))

        assertEquals(40, graph.numberOfRanges(), "numOfRanges not 40: ${graph.numberOfRanges()}")

        val ranges = graph.ranges()

        assertTrue(ranges.isSorted(), "ranges not sorted")

    }

    @Test
    fun testMultipleFilesHaplotypeGraph() {

        val graph = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFile, TestExtension.smallseqLineBHvcfFile))

        assertEquals(38, graph.numberOfRanges(), "numOfRanges not 38: ${graph.numberOfRanges()}")

        assertEquals(2, graph.numberOfSamples(), "numOfSamples not 2: ${graph.numberOfSamples()}")

        val ranges = graph.ranges()

        assertTrue(ranges.isSorted(), "ranges not sorted")

        assertEquals(graph.numberOfRanges(), ranges.size, "ranges size not equal to numberOfRanges")

        ranges.forEach { range ->
            println("$range: ${graph.hapIdToSamples(range)}")
            graph.hapIdToSamples(range).forEach { (hapId, samples) ->
                println("$hapId: ${samples.joinToString(", ")}")
            }
        }

    }

    @Test
    fun testBadInputHaplotypeGraph() {
        val graph = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFileBadAltTag))
        assert(graph == null || graph.numberOfRanges() == 0) { "graph not null or empty" }
    }

    /**
     * Returns true if the list of ReferenceRange is sorted.
     */
    private fun List<ReferenceRange>.isSorted(): Boolean {
        for (i in 1 until this.size) {
            if (this[i - 1] > this[i]) {
                return false
            }
        }
        return true
    }

}