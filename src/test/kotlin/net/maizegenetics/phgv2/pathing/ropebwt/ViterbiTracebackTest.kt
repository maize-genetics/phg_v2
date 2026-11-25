package net.maizegenetics.phgv2.pathing.ropebwt

import org.junit.jupiter.api.Assertions.assertArrayEquals
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.Test
import kotlin.math.ln
import kotlin.random.Random

/**
 * The block-wise traceback must give exactly the path the whole table gives, ties included, for every
 * recursion that uses it and every block length.
 */
class ViterbiTracebackTest {

    private val unlimited = Long.MAX_VALUE

    /**
     * Emissions drawn from a few integer values, so that many states tie exactly and the tie rules are
     * exercised, or from a continuous range when [ties] is false.
     */
    private fun emissions(positions: Int, states: Int, seed: Int, ties: Boolean): Array<DoubleArray> {
        val random = Random(seed)
        return Array(positions) {
            DoubleArray(states) { if (ties) -random.nextInt(3).toDouble() else -4.0 * random.nextDouble() }
        }
    }

    /** Budgets from "one position's worth" up to just short of the whole table. */
    private fun budgets(positions: Int, states: Int): List<Long> =
        listOf(1L, states.toLong(), 3L * states, (positions - 2L) * states)

    @Test
    fun testDiploidFastPathIsUnchangedByBlocks() {
        val nParents = 4
        val states = nParents * nParents
        for (seed in 0 until 20) for (ties in listOf(true, false)) {
            val positions = 30 + seed
            val emission = emissions(positions, states, seed, ties)
            val lnNoSwitch = DoubleArray(positions) { ln(0.9) }
            val lnSwitch = DoubleArray(positions) { if (it % 3 == 0) ln(0.05) else ln(0.02) }
            val initial = DoubleArray(states) { -ln(states.toDouble()) }
            fun run(budget: Long) = ViterbiHMM(0.0, 0.99, 0.98, maxStoredBackPointers = budget)
                .viterbiOptimizedForDiploid(nParents, positions, initial, { emission[it] }, lnNoSwitch, lnSwitch)
            val whole = run(unlimited)
            for (budget in budgets(positions, states)) {
                val blocked = run(budget)
                assertArrayEquals(whole.first, blocked.first, "seed $seed, ties $ties, budget $budget")
                assertEquals(whole.second, blocked.second, 0.0)
            }
        }
    }

    @Test
    fun testGeneralScanIsUnchangedByBlocks() {
        val states = 5
        for (seed in 0 until 20) for (ties in listOf(true, false)) {
            val positions = 25 + seed
            val emission = emissions(positions, states, seed, ties)
            val random = Random(1000 + seed)
            val transitions = DoubleArray(states * states) { if (ties) -random.nextInt(2).toDouble() else -random.nextDouble() }
            val initial = DoubleArray(states) { -ln(states.toDouble()) }
            fun run(budget: Long) = ViterbiHMM(0.5, 0.99, 0.98, maxStoredBackPointers = budget)
                .viterbiOptimized(states, positions, initial, transitions, { emission[it] })
            val whole = run(unlimited)
            for (budget in budgets(positions, states)) {
                val blocked = run(budget)
                assertArrayEquals(whole.first, blocked.first, "seed $seed, ties $ties, budget $budget")
                assertEquals(whole.second, blocked.second, 0.0)
            }
        }
    }

    @Test
    fun testHaploidPathIsUnchangedByBlocks() {
        val states = 6
        for (seed in 0 until 20) for (ties in listOf(true, false)) {
            val positions = 40 + seed
            val emission = emissions(positions, states, seed, ties)
            val initial = DoubleArray(states) { -ln(states.toDouble()) }
            fun run(budget: Long) = ViterbiHMM(1.0, 0.95, 0.98, maxStoredBackPointers = budget)
                .viterbiOptimizedForHaploid(states, positions, initial, { emission[it] })
            val whole = run(unlimited)
            for (budget in budgets(positions, states)) {
                val blocked = run(budget)
                assertArrayEquals(whole.first, blocked.first, "seed $seed, ties $ties, budget $budget")
                assertEquals(whole.second, blocked.second, 0.0)
            }
        }
    }

    @Test
    fun testBlockedTracebackHoldsOneBlockAndStepsEachPositionTwice() {
        // A path of 10,000 positions x 50 states would need 500,000 backpointers whole. With a budget
        // of 5,000 it should hold 100 rows at a time (the budget allows 100, and sqrt(9,999) is 100),
        // and call step twice for each of the 9,999 steps.
        val states = 50
        val positions = 10_000
        val random = Random(7)
        val emission = Array(positions) { DoubleArray(states) { -random.nextDouble() } }
        var steps = 0
        var widestRow = 0
        fun run(budget: Long): Pair<IntArray, Double> {
            steps = 0
            return traceBestPath(states, positions, budget,
                initialize = { scores -> for (s in 0 until states) scores[s] = emission[0][s] }
            ) { position, previous, current, backPointers ->
                steps++
                widestRow = maxOf(widestRow, backPointers.size)
                var best = 0
                for (s in 1 until states) if (previous[s] > previous[best]) best = s
                for (s in 0 until states) {
                    val stay = previous[s] + ln(0.99)
                    val move = previous[best] + ln(0.01 / (states - 1))
                    if (stay >= move) { current[s] = stay + emission[position][s]; backPointers[s] = s }
                    else { current[s] = move + emission[position][s]; backPointers[s] = best }
                }
            }
        }
        val whole = run(unlimited)
        assertEquals(positions - 1, steps)
        val blocked = run(5_000)
        assertEquals(2 * (positions - 1), steps)
        assertEquals(states, widestRow)
        assertEquals(100, backPointerBlockLength(states, positions - 1, 5_000))
        assertArrayEquals(whole.first, blocked.first)
        assertEquals(whole.second, blocked.second, 0.0)
    }

    @Test
    fun testBlockLength() {
        // everything fits: one block
        assertEquals(1000, backPointerBlockLength(10, 1000, 10_000))
        // the budget allows 50 positions, more than sqrt(1000) = 32
        assertEquals(50, backPointerBlockLength(10, 1000, 500))
        // the budget allows 1 position, so sqrt(1000), rounded up, is used
        assertEquals(32, backPointerBlockLength(10, 1000, 10))
        // the 82-founder diploid case past the old limit: 400,000 steps at 6,724 states
        val length = backPointerBlockLength(82 * 82, 400_000, 1L shl 24)
        assertEquals((1L shl 24).toInt() / (82 * 82), length)
        assertTrue(length.toLong() * 82 * 82 <= 1L shl 24)
        // a single position needs no steps at all
        val single = traceBestPath(3, 1, 1, initialize = { it[0] = -1.0; it[1] = 0.0; it[2] = -2.0 }) { _, _, _, _ ->
            error("no steps for one position")
        }
        assertArrayEquals(intArrayOf(1), single.first)
    }
}
