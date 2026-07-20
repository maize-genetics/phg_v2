package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.utils.Position
import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test

class PathIntervalsTest {

    /** Builds a path on chr1 from (bin, call) pairs. */
    private fun <T> path(vararg bins: Pair<Int, T>) = bins.map { (bin, call) -> Pair(Position("chr1", bin), call) }

    @Test
    fun testEmptyPath() {
        assertTrue(pathToIntervals(emptyList<Pair<Position, String>>(), 100).isEmpty())
        assertTrue(pathToIntervals(emptyList<Pair<Position, String>>(), 100, mergeAdjacent = false).isEmpty())
    }

    @Test
    fun testSingleBin() {
        // A contig with one bin has no midpoint to cut at: it spans [0, pos * binSize).
        val intervals = pathToIntervals(path(50 to "lineA:0"), 100)
        assertEquals(listOf(PathInterval("chr1", 0, 5000, "lineA:0")), intervals)

        // Expanding cannot split a single bin either.
        assertEquals(intervals, pathToIntervals(path(50 to "lineA:0"), 100, mergeAdjacent = false))
    }

    @Test
    fun testNoMergeWhenAllCallsDiffer() {
        // Bins 10,20,30,40,50 at binSize=100. Cuts fall at the bin midpoints:
        //   (10+20)/2*100=1500, (20+30)/2*100=2500, (30+40)/2*100=3500, (40+50)/2*100=4500
        // and the last interval ends at 50*100=5000.
        val intervals = pathToIntervals(
            path(10 to "a", 20 to "b", 30 to "a", 40 to "b", 50 to "a"), 100
        )
        assertEquals(
            listOf(
                PathInterval("chr1", 0, 1500, "a"),
                PathInterval("chr1", 1500, 2500, "b"),
                PathInterval("chr1", 2500, 3500, "a"),
                PathInterval("chr1", 3500, 4500, "b"),
                PathInterval("chr1", 4500, 5000, "a")
            ),
            intervals
        )
    }

    @Test
    fun testExpandBinsMatchesOneRecordPerBin() {
        // Same bins, but every call identical: with mergeAdjacent = false the cuts are still made,
        // so the output is one record per bin with the same coordinates as the test above.
        val intervals = pathToIntervals(
            path(10 to "a", 20 to "a", 30 to "a", 40 to "a", 50 to "a"), 100, mergeAdjacent = false
        )
        assertEquals(5, intervals.size)
        assertEquals(listOf(0, 1500, 2500, 3500, 4500), intervals.map { it.start })
        assertEquals(listOf(1500, 2500, 3500, 4500, 5000), intervals.map { it.end })
        assertTrue(intervals.all { it.call == "a" })
    }

    @Test
    fun testMergeCollapsesIdenticalCalls() {
        val intervals = pathToIntervals(
            path(10 to "a", 20 to "a", 30 to "a", 40 to "a", 50 to "a"), 100
        )
        assertEquals(listOf(PathInterval("chr1", 0, 5000, "a")), intervals)
    }

    @Test
    fun testMergeWithRunsAtBothEnds() {
        // Runs a,a | b | c,c -- the cut between runs is the midpoint of the flanking bins.
        val intervals = pathToIntervals(
            path(10 to "a", 20 to "a", 30 to "b", 40 to "c", 50 to "c"), 100
        )
        assertEquals(
            listOf(
                PathInterval("chr1", 0, 2500, "a"),
                PathInterval("chr1", 2500, 3500, "b"),
                PathInterval("chr1", 3500, 5000, "c")
            ),
            intervals
        )
    }

    @Test
    fun testIntervalsAreGaplessAndNonOverlapping() {
        val bins = listOf(3, 11, 12, 40, 41, 90)
        val calls = listOf("a", "a", "b", "b", "a", "a")
        val input = bins.zip(calls).map { (bin, call) -> Pair(Position("chr1", bin), call) }

        for (merge in listOf(true, false)) {
            val intervals = pathToIntervals(input, 256, mergeAdjacent = merge)
            assertEquals(0, intervals.first().start, "first interval must start at 0 (merge=$merge)")
            assertEquals(90 * 256, intervals.last().end, "last interval must end at pos * binSize (merge=$merge)")
            intervals.zipWithNext { left, right ->
                assertEquals(left.end, right.start, "intervals must be gapless (merge=$merge)")
                assertTrue(left.start < left.end, "intervals must be non-empty (merge=$merge)")
            }
        }
    }

    @Test
    fun testDiploidPhaseSwitchIsNotMerged() {
        // (lineA:0, lineB:0) and (lineB:0, lineA:0) are the same genotype but opposite phase.
        // They must not merge: parent1 and parent2 both have to match.
        val ab = Pair("lineA:0", "lineB:0")
        val ba = Pair("lineB:0", "lineA:0")
        val intervals = pathToIntervals(
            path(10 to ab, 20 to ab, 30 to ba, 40 to ba), 100
        )
        assertEquals(
            listOf(
                PathInterval("chr1", 0, 2500, ab),
                PathInterval("chr1", 2500, 4000, ba)
            ),
            intervals
        )
    }

    @Test
    fun testDiploidIdenticalPairsAreMerged() {
        val ab = Pair("lineA:0", "lineB:0")
        val intervals = pathToIntervals(path(10 to ab, 20 to ab, 30 to ab), 100)
        assertEquals(listOf(PathInterval("chr1", 0, 3000, ab)), intervals)
    }
}
