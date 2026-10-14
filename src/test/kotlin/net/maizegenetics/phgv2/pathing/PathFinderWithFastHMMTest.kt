package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals

/**
 * Tests [PathFinderWithFastHMM] on the small-sequence graph of Ref, LineA and LineB: two contigs of
 * twenty ranges each, every founder with a distinct haplotype in every range, except the last range of
 * each contig, where only Ref has one.
 */
class PathFinderWithFastHMMTest {

    private val graph = HaplotypeGraph(listOf(TestExtension.smallseqRefHvcfFile,
        TestExtension.smallseqLineAHvcfFile, TestExtension.smallseqLineBHvcfFile))
    private val ranges = graph.rangesByContig()
    private val lineA = SampleGamete("LineA")
    private val lineB = SampleGamete("LineB")
    private val ref = SampleGamete("Ref")

    /** Read mappings for one range: [count] reads hitting only [gamete]'s haplotype, for each pair. */
    private fun reads(range: ReferenceRange, vararg counts: Pair<SampleGamete, Int>): Map<List<String>, Int> =
        counts.associate { (gamete, count) -> listOf(graph.sampleToHapId(range, gamete)!!) to count }

    private fun pathFinder(isHaploid: Boolean, probSame: Double) =
        PathFinderWithFastHMM(graph, isHaploid, probCorrect = 0.99, sameGameteProbability = probSame)

    @Test
    fun haploidPathFollowsTheReadsAndFillsRangesWithoutReads() {
        val readMap = mutableMapOf<ReferenceRange, Map<List<String>, Int>>()
        ranges["1"]!!.forEachIndexed { index, range ->
            if (index < 10) readMap[range] = reads(range, lineA to 5)
            else if (index < 19) readMap[range] = reads(range, lineB to 5)
        }
        ranges["2"]!!.take(19).forEach { readMap[it] = reads(it, lineA to 5) }

        val path = pathFinder(true, 0.9).findBestPath(readMap).first
        assertEquals(graph.ranges(), path.map { it.refRange }, "one node per range, in order")
        for (node in path) {
            val index = ranges[node.refRange.contig]!!.indexOf(node.refRange)
            val expected = if (node.refRange.contig == "1" && index >= 10) lineB else lineA
            assertEquals(listOf(expected), node.sampleGametes, "${node.refRange}")
        }
    }

    /**
     * Ranges between two ranges with reads take the earlier state up to the midpoint and the later
     * one after it, and a contig without reads gets no path.
     */
    @Test
    fun diploidPathSplitsUnreadRangesAtTheMidpoint() {
        val chr1 = ranges["1"]!!
        val readMap = mutableMapOf<ReferenceRange, Map<List<String>, Int>>()
        for (index in 0..2) readMap[chr1[index]] = reads(chr1[index], lineA to 10, lineB to 10)
        for (index in 7..9) readMap[chr1[index]] = reads(chr1[index], ref to 20)

        val path = pathFinder(false, 0.9).findBestPath(readMap).first
        assertEquals(chr1, path.map { it.refRange }, "chr1 only: chr2 has no reads")
        path.forEachIndexed { index, node ->
            // ranges 3 and 4 are at or before the midpoint of 2 and 7
            val expected = if (index <= 4) setOf(lineA, lineB) else setOf(ref)
            assertEquals(expected, node.sampleGametes.toSet(), "range index $index")
            assertEquals(2, node.sampleGametes.size)
        }
    }

    /**
     * In the last range only Ref has a haplotype, so every read matches Ref alone. LineA is absent
     * there, so (LineA, Ref) must score the same as (Ref, Ref) and the heterozygous path continues
     * rather than switching to homozygous Ref for one range.
     */
    @Test
    fun aFounderWithoutAHaplotypeDoesNotForceHomozygosity() {
        val chr2 = ranges["2"]!!
        val readMap = chr2.associateWith { range ->
            if (graph.sampleToHapId(range, lineA) == null) reads(range, ref to 20)
            else reads(range, lineA to 10, ref to 10)
        }
        val path = pathFinder(false, 0.9).findBestPath(readMap).first
        for (node in path) assertEquals(setOf(lineA, ref), node.sampleGametes.toSet(), "${node.refRange}")
    }

    /**
     * Ranges without reads count in the transition: one LineB read right after three LineA reads is
     * not enough to switch, but the same read eighteen ranges later is, because the probability of
     * staying is 0.99 per range.
     */
    @Test
    fun rangesWithoutReadsCountTowardTheSwitchProbability() {
        val chr1 = ranges["1"]!!
        val chr2 = ranges["2"]!!
        val readMap = mapOf(
            chr1[0] to reads(chr1[0], lineA to 3), chr1[1] to reads(chr1[1], lineB to 1),
            chr2[0] to reads(chr2[0], lineA to 3), chr2[18] to reads(chr2[18], lineB to 1)
        )
        val path = pathFinder(true, 0.99).findBestPath(readMap).first
        val byRange = path.associate { it.refRange to it.sampleGametes.single() }
        assertEquals(lineA, byRange[chr1[1]], "adjacent read does not switch")
        assertEquals(lineA, byRange[chr2[9]], "at the midpoint, still LineA")
        assertEquals(lineB, byRange[chr2[10]], "after the midpoint, LineB")
        assertEquals(lineB, byRange[chr2[18]], "distant read switches")
    }
}
