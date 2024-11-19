package net.maizegenetics.phgv2.api

import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
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

        assertTrue(graph.numberOfSamples() == 1, "numOfSamples not 1: ${graph.numberOfSamples()}")

        assertTrue(
            graph.checksum == "1566df5fad42a5752c5720e1aef039ab",
            "checksum not 1566df5fad42a5752c5720e1aef039ab: ${graph.checksum}"
        )

    }

    @Test
    fun testMultipleFilesHaplotypeGraph() {
        //include RefHvcf to test what happens when a sample has no haplotype in at least one ReferenceRange
        val graph = HaplotypeGraph(
            listOf(
                TestExtension.smallseqLineAHvcfFile,
                TestExtension.smallseqLineBHvcfFile,
                TestExtension.smallseqRefHvcfFile
            )
        )

        assertTrue(
            graph.checksum == "1e3d5c6ed7d43c1725542d826a515b16",
            "checksum not 1e3d5c6ed7d43c1725542d826a515b16: ${graph.checksum}"
        )

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
        // ALT tag does not contain sample name
        assertThrows<IllegalStateException> {
            HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFileBadAltTag))
        }
        // This is the assert needed if we switch back to co-routines in HaplotypGraph
        // Coroutines suppress the exception, hence we need to check the graph for null
        // val graph = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFileBadAltTag))
        // assert(graph == null || graph.numberOfRanges() == 0) { "graph not null or empty" }
    }

    @Test
    fun testHapIdToRefRangeMap() {
        val graph = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFile, TestExtension.smallseqLineBHvcfFile))
        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()
        assertEquals(76, hapIdToRefRangeMap.size, "hapIdToRefRangeMap size not 38: ${hapIdToRefRangeMap.size}")
        for (range in graph.ranges()) {
            val hapIdLineA = graph.sampleToHapId(range, SampleGamete("LineA"))
            val refRangeLineA = hapIdToRefRangeMap[hapIdLineA]
            assertTrue(refRangeLineA != null, "$hapIdLineA does not map to a reference range")
            assertEquals(1, refRangeLineA.size)
            assertEquals(
                range,
                refRangeLineA[0],
                "hapIdToRefRangeMap does not contain correct range $refRangeLineA for LineA"
            )
            val hapIdLineB = graph.sampleToHapId(range, SampleGamete("LineB"))
            val refRangeLineB = hapIdToRefRangeMap[hapIdLineB]
            assertTrue(refRangeLineB != null, "$hapIdLineB does not map to a reference range")
            assertEquals(1, refRangeLineB.size)
            assertEquals(
                range,
                refRangeLineB[0],
                "hapIdToRefRangeMap does not contain correct range $refRangeLineB for LineB"
            )
        }
    }

    @Test
    fun testRefRangeToHapIdMap() {
        val graph = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFile, TestExtension.smallseqLineBHvcfFile))
        val refRangeToHapidMap = graph.refRangeToHapIdMap()

        assertEquals(38, refRangeToHapidMap.size, "refRangeToHapidMap size not 38: ${refRangeToHapidMap.size}")
        for (range in graph.ranges()) {
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
        for ((index, range) in graph.ranges().withIndex()) {
            val mappedIndex = refRangeToIndexMap[range]
            assertEquals(
                index,
                mappedIndex,
                "refRangeToIndexMap does not contain correct index $index for range $range"
            )
        }
    }

    @Test
    fun testRefRangeIdToHapIdMap() {
        val graph = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFile, TestExtension.smallseqLineBHvcfFile))
        val refRangeIdToHapIdMap = graph.refRangeIdToHapIdMap()

        assertEquals(38, refRangeIdToHapIdMap.size, "refRangeIdToHapIdMap size not 38: ${refRangeIdToHapIdMap.size}")
        for ((index, range) in graph.ranges().withIndex()) {
            val hapIdLineA = graph.sampleToHapId(range, SampleGamete("LineA"))
            val hapIdLineB = graph.sampleToHapId(range, SampleGamete("LineB"))
            val hapIds = refRangeIdToHapIdMap[index]!!
            assertTrue(
                hapIdLineA in hapIds,
                "refRangeIdToHapIdMap does not contain correct hapId $hapIdLineA for LineA"
            )
            assertTrue(
                hapIdLineB in hapIds,
                "refRangeIdToHapIdMap does not contain correct hapId $hapIdLineB for LineB"
            )
        }
    }

    @Test
    fun testSampleGameteToHaplotypeId() {
        val graph = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFile, TestExtension.smallseqLineBHvcfFile))
        val sampleGametes = graph.sampleGametesInGraph()

        val expectedHapIds = listOf(
            "12f0cec9102e84a161866e37072443b7",
            "4fc7b8af32ddd74e07cb49d147ef1938"
        )

        sampleGametes.forEachIndexed { index, sampleGamete ->
            val testHapIds = graph.sampleGameteToHaplotypeId(sampleGamete)
            assertEquals(38, testHapIds.size)
            assertTrue { testHapIds[0] is String }
            assertEquals(expectedHapIds[index], testHapIds[0])
        }
    }

    @Test
    fun testHapIdsToSampleGametes() {
        val graph = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFile, TestExtension.smallseqLineBHvcfFile))
        val hapIdsToSampleGametes = graph.hapIdsToSampleGametes()

        val ranges = graph.ranges()

        for (range in ranges) {
            val hapIdToSampleGametes = graph.hapIdToSampleGametes(range)
            for ((hapId, sampleGametes) in hapIdToSampleGametes) {
                val testSampleGametes = hapIdsToSampleGametes[hapId]
                assertEquals(sampleGametes, testSampleGametes)
            }
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

    @Test
    fun testRefChecksum() {

        val hvcfFile = "${TestExtension.smallSeqInputDir}LineB_kmer_index_test.h.vcf"

        val graph = HaplotypeGraph(listOf(hvcfFile))

        var range = ReferenceRange("1", 1001, 5500)
        assertTrue(
            graph.refChecksum(range) == "57705b1e2541c7634ea59a48fc52026f",
            "checksum not 57705b1e2541c7634ea59a48fc52026f: ${graph.checksum}"
        )

        range = ReferenceRange("2", 6501, 11000)
        assertTrue(
            graph.refChecksum(range) == "b8843efbd6adaa261a01518dc2a39aa2",
            "checksum not b8843efbd6adaa261a01518dc2a39aa2: ${graph.checksum}"
        )

        range = ReferenceRange("2", 49501, 50500)
        assertTrue(
            graph.refChecksum(range) == "39f96726321b329964435865b3694fd2",
            "checksum not 39f96726321b329964435865b3694fd2: ${graph.checksum}"
        )

    }

}