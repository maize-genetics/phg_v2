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

        // tests hapIdToSamples() method

        var hapIdToSamples = graph.hapIdToSamples(ranges[0])

        assertEquals(2, hapIdToSamples.size, "hapIdToSamples size not 2: ${hapIdToSamples.size}")

        var hapid = "12f0cec9102e84a161866e37072443b7"
        var samples = hapIdToSamples[hapid]
        assertEquals("LineA", samples?.get(0), "sample not LineA: ${samples?.get(0)}")

        hapid = "4fc7b8af32ddd74e07cb49d147ef1938"
        samples = hapIdToSamples[hapid]
        assertEquals("LineB", samples?.get(0), "sample not LineB: ${samples?.get(0)}")

        hapIdToSamples = graph.hapIdToSamples(ranges[ranges.size - 1])
        hapid = "0eb9029f3896313aebc69c8489923141"
        samples = hapIdToSamples[hapid]
        assertEquals("LineA", samples?.get(0), "sample not LineA: ${samples?.get(0)}")

        hapid = "5031218d4ac709dd51a946acd0550356"
        samples = hapIdToSamples[hapid]
        assertEquals("LineB", samples?.get(0), "sample not LineB: ${samples?.get(0)}")

        // tests for sampleToHapId() method

        var checksum = graph.sampleToHapId(ranges[0], "LineA")
        assertEquals(
            "12f0cec9102e84a161866e37072443b7",
            checksum,
            "sampleToHapId: checksum not 12f0cec9102e84a161866e37072443b7: $checksum"
        )

        checksum = graph.sampleToHapId(ranges[0], "LineB")
        assertEquals(
            "4fc7b8af32ddd74e07cb49d147ef1938",
            checksum,
            "sampleToHapId: checksum not 4fc7b8af32ddd74e07cb49d147ef1938: $checksum"
        )

        checksum = graph.sampleToHapId(ranges[ranges.size - 1], "LineA")
        assertEquals(
            "0eb9029f3896313aebc69c8489923141",
            checksum,
            "sampleToHapId: checksum not 0eb9029f3896313aebc69c8489923141: $checksum"
        )

        checksum = graph.sampleToHapId(ranges[ranges.size - 1], "LineB")
        assertEquals(
            "5031218d4ac709dd51a946acd0550356",
            checksum,
            "sampleToHapId: checksum not 5031218d4ac709dd51a946acd0550356: $checksum"
        )

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