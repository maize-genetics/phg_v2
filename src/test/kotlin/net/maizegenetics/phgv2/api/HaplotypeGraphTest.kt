package net.maizegenetics.phgv2.api

import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import kotlin.test.assertEquals
import kotlin.test.assertTrue
import kotlin.test.fail

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
        //include RefHvcf to test what happens when a sample has no haplotype in at least one ReferenceRange
        val graph = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFile, TestExtension.smallseqLineBHvcfFile, TestExtension.smallseqRefHvcfFile))

        assertEquals(40, graph.numberOfRanges(), "numOfRanges not 40: ${graph.numberOfRanges()}")

        assertEquals(3, graph.numberOfSamples(), "numOfSamples not 3: ${graph.numberOfSamples()}")

        val ranges = graph.ranges()

        assertTrue(ranges.isSorted(), "ranges not sorted")

        assertEquals(graph.numberOfRanges(), ranges.size, "ranges size not equal to numberOfRanges")

        // tests hapIdToSamples() method

        var hapIdToSamples = graph.hapIdToSampleGametes(ranges[0])

        assertEquals(3, hapIdToSamples.size, "hapIdToSamples size not 3: ${hapIdToSamples.size}")

        var hapid = "12f0cec9102e84a161866e37072443b7"
        var samples = hapIdToSamples[hapid]
        assertEquals("LineA", samples?.get(0)?.name, "sample not LineA: ${samples?.get(0)?.name}")

        hapid = "4fc7b8af32ddd74e07cb49d147ef1938"
        samples = hapIdToSamples[hapid]
        assertEquals("LineB", samples?.get(0)?.name, "sample not LineB: ${samples?.get(0)?.name}")

        hapIdToSamples = graph.hapIdToSampleGametes(ranges[ranges.size - 2])
        hapid = "0eb9029f3896313aebc69c8489923141"
        samples = hapIdToSamples[hapid]
        assertEquals("LineA", samples?.get(0)?.name, "sample not LineA: ${samples?.get(0)?.name}")

        hapid = "5031218d4ac709dd51a946acd0550356"
        samples = hapIdToSamples[hapid]
        assertEquals("LineB", samples?.get(0)?.name, "sample not LineB: ${samples?.get(0)?.name}")

        // tests for sampleToHapId() method

        var checksum = graph.sampleToHapId(ranges[0], SampleGamete("LineA"))
        assertEquals(
            "12f0cec9102e84a161866e37072443b7",
            checksum,
            "sampleToHapId: checksum not 12f0cec9102e84a161866e37072443b7: $checksum"
        )

        checksum = graph.sampleToHapId(ranges[0], SampleGamete("LineB"))
        assertEquals(
            "4fc7b8af32ddd74e07cb49d147ef1938",
            checksum,
            "sampleToHapId: checksum not 4fc7b8af32ddd74e07cb49d147ef1938: $checksum"
        )

        checksum = graph.sampleToHapId(ranges[ranges.size - 2], SampleGamete("LineA"))
        assertEquals(
            "0eb9029f3896313aebc69c8489923141",
            checksum,
            "sampleToHapId: checksum not 0eb9029f3896313aebc69c8489923141: $checksum"
        )

        checksum = graph.sampleToHapId(ranges[ranges.size - 2], SampleGamete("LineB"))
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

    @Test
    fun testHapIdToRefRangeMap() {
        val graph = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFile, TestExtension.smallseqLineBHvcfFile))
        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()
        assertEquals(76, hapIdToRefRangeMap.size, "hapIdToRefRangeMap size not 38: ${hapIdToRefRangeMap.size}")
        for (range in graph.ranges()) {
            val hapIdLineA = graph.sampleToHapId(range, SampleGamete("LineA"))
            val refRangeLineA = hapIdToRefRangeMap[hapIdLineA]
            assertEquals(range, refRangeLineA, "hapIdToRefRangeMap does not contain correct range $refRangeLineA for LineA")
            val hapIdLineB = graph.sampleToHapId(range, SampleGamete("LineB"))
            val refRangeLineB = hapIdToRefRangeMap[hapIdLineB]
            assertEquals(range, refRangeLineB, "hapIdToRefRangeMap does not contain correct range $refRangeLineB for LineB")
        }
    }

    @Test
    fun testRefRangeToHapIdMap() {
        val graph = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFile, TestExtension.smallseqLineBHvcfFile))
        val refRangeToHapidMap = graph.refRangeToHapIdMap()

        assertEquals(38, refRangeToHapidMap.size, "refRangeToHapidMap size not 38: ${refRangeToHapidMap.size}")
        for(range in graph.ranges()) {
            val hapIds = refRangeToHapidMap[range]!!.keys
            val hapIdLineA = graph.sampleToHapId(range, SampleGamete("LineA"))
            assertTrue(hapIdLineA in hapIds, "refRangeToHapidMap does not contain correct hapId $hapIdLineA for LineA")

            val hapIdLineB = graph.sampleToHapId(range, SampleGamete("LineB"))
            assertTrue(hapIdLineB in hapIds, "refRangeToHapidMap does not contain correct hapId $hapIdLineB for LineB")
        }

    }

    @Test
    fun testRefRangeToIndexMap() {
        val graph = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFile, TestExtension.smallseqLineBHvcfFile))
        val refRangeToIndexMap = graph.refRangeToIndexMap()

        assertEquals(38, refRangeToIndexMap.size, "refRangeToIndexMap size not 38: ${refRangeToIndexMap.size}")
        for((index,range) in graph.ranges().withIndex()) {
            val mappedIndex = refRangeToIndexMap[range]
            assertEquals(index, mappedIndex, "refRangeToIndexMap does not contain correct index $index for range $range")
        }
    }

    @Test
    fun testRefRangeIdToHapIdMap() {
        val graph = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFile, TestExtension.smallseqLineBHvcfFile))
        val refRangeIdToHapIdMap = graph.refRangeIdToHapIdMap()

        assertEquals(38, refRangeIdToHapIdMap.size, "refRangeIdToHapIdMap size not 38: ${refRangeIdToHapIdMap.size}")
        for((index, range) in graph.ranges().withIndex()) {
            val hapIdLineA = graph.sampleToHapId(range, SampleGamete("LineA"))
            val hapIdLineB = graph.sampleToHapId(range, SampleGamete("LineB"))
            val hapIds = refRangeIdToHapIdMap[index]!!
            assertTrue(hapIdLineA in hapIds, "refRangeIdToHapIdMap does not contain correct hapId $hapIdLineA for LineA")
            assertTrue(hapIdLineB in hapIds, "refRangeIdToHapIdMap does not contain correct hapId $hapIdLineB for LineB")
        }
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