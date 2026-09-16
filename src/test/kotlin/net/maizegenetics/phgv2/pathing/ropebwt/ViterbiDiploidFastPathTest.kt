package net.maizegenetics.phgv2.pathing.ropebwt

import org.junit.jupiter.api.Test
import kotlin.math.ln
import kotlin.random.Random
import kotlin.test.assertEquals

/**
 * The O(nParents^2) diploid recursion must agree with the general O(nStates^2) scan it replaces.
 * At an inbreeding coefficient of zero the transition is a Kronecker product of two haploid
 * transitions, so the maximisation separates; these tests check that the separation is exact.
 */
class ViterbiDiploidFastPathTest {

    private fun transitionMatrix(nParents: Int, pNoSwitch: Double): DoubleArray {
        val nStates = nParents * nParents
        val matrix = DoubleArray(nStates * nStates)
        val calculator = DiploidTransitionProbability(pNoSwitch, 0.0, nParents)
        var pointer = 0
        for (a in 0 until nParents) for (b in 0 until nParents)
            for (c in 0 until nParents) for (d in 0 until nParents)
                matrix[pointer++] = calculator.calculateLn(Pair(a, b), Pair(c, d))
        return matrix
    }

    /** The factorisation the fast path relies on: pnn, psn and pss are a product of two haploids. */
    @Test
    fun theTransitionIsAKroneckerProductOfTwoHaploidTransitions() {
        val nParents = 5
        val pNoSwitch = 0.9999
        val matrix = transitionMatrix(nParents, pNoSwitch)
        val lnNo = ln(pNoSwitch)
        val lnSw = ln((1.0 - pNoSwitch) / (nParents - 1))
        for (a in 0 until nParents) for (b in 0 until nParents)
            for (c in 0 until nParents) for (d in 0 until nParents) {
                val expected = (if (a == c) lnNo else lnSw) + (if (b == d) lnNo else lnSw)
                val from = a * nParents + b
                val to = c * nParents + d
                assertEquals(expected, matrix[from * nParents * nParents + to], 1e-12)
            }
    }

    private fun compare(nParents: Int, positions: Int, seed: Int, discrete: Boolean) {
        val random = Random(seed)
        val nStates = nParents * nParents
        val pNoSwitch = 0.999999999
        val hmm = ViterbiHMM(0.0, pNoSwitch, 0.98)
        // discrete emissions are multiples of one constant, so exact ties are common -- which is
        // where a separated maximisation is most likely to diverge from the full scan
        val emissions = Array(positions) { DoubleArray(nStates) {
            if (discrete) random.nextInt(0, 4) * ln(0.98) else random.nextDouble(-6.0, 0.0)
        } }
        val initial = DoubleArray(nStates) { if (discrete) 0.0 else random.nextDouble(-1.0, 0.0) }
        val emissionFn = { p: Int -> emissions[p] }

        val general = hmm.viterbiOptimized(nStates, positions, initial,
            transitionMatrix(nParents, pNoSwitch), emissionFn)
        val fast = hmm.viterbiOptimizedForDiploid(nParents, positions, initial, emissionFn)

        assertEquals(general.second, fast.second, 1e-9,
            "path log probability, nParents=$nParents positions=$positions discrete=$discrete")
        assertEquals(general.first.toList(), fast.first.toList(),
            "state path, nParents=$nParents positions=$positions discrete=$discrete")
    }

    @Test
    fun matchesTheGeneralScanOnContinuousEmissions() {
        for (n in listOf(2, 3, 6, 9)) for (p in listOf(5, 60)) compare(n, p, n * 100 + p, false)
    }

    @Test
    fun matchesTheGeneralScanWhenTiesAreEverywhere() {
        for (n in listOf(2, 3, 6, 9)) for (p in listOf(30, 120)) compare(n, p, n * 7 + p, true)
    }

    /**
     * At F = 1 the general scan and the homozygous-only fast path must choose the same states.
     * The emissions are arbitrary over all nParents^2 states, so a fast path that read the wrong
     * cells of the emission array -- the diagonal is the only correct choice -- would diverge.
     */
    private fun compareInbred(nParents: Int, positions: Int, seed: Int, discrete: Boolean) {
        val random = Random(seed)
        val nStates = nParents * nParents
        val pNoSwitch = 0.999999
        val hmm = ViterbiHMM(1.0, pNoSwitch, 0.98)
        val emissions = Array(positions) { DoubleArray(nStates) {
            if (discrete) random.nextInt(0, 4) * ln(0.98) else random.nextDouble(-6.0, 0.0)
        } }
        val emissionFn = { p: Int -> emissions[p] }

        // Initial probabilities exactly as findDiploidPath builds them at F = 1: homozygous
        // states carry 1/nParents, heterozygous states are impossible.
        val initial = DoubleArray(nStates) { -1.0e6 }
        for (i in 0 until nParents) initial[i * nParents + i] = ln(1.0 / nParents)

        val calculator = DiploidTransitionProbability(pNoSwitch, 1.0, nParents)
        val matrix = DoubleArray(nStates * nStates)
        var ptr = 0
        for (a in 0 until nParents) for (b in 0 until nParents)
            for (c in 0 until nParents) for (d in 0 until nParents)
                matrix[ptr++] = calculator.calculateLn(Pair(a, b), Pair(c, d))

        val general = hmm.viterbiOptimized(nStates, positions, initial, matrix, emissionFn)

        val diagonalInit = DoubleArray(nParents) { initial[it * nParents + it] }
        val diagonalEmission = { p: Int ->
            val full = emissionFn(p); DoubleArray(nParents) { full[it * nParents + it] }
        }
        val fastHap = hmm.viterbiOptimizedForHaploid(nParents, positions, diagonalInit, diagonalEmission)
        val fast = IntArray(fastHap.first.size) { fastHap.first[it] * nParents + fastHap.first[it] }

        assertEquals(general.first.toList(), fast.toList(),
            "state path, nParents=$nParents positions=$positions discrete=$discrete")
        assertEquals(general.second, fastHap.second, 1e-9,
            "path log probability, nParents=$nParents positions=$positions discrete=$discrete")
    }

    @Test
    fun theInbredFastPathMatchesTheGeneralScanAtFOne() {
        for (n in listOf(2, 3, 6, 9)) for (p in listOf(5, 60)) compareInbred(n, p, n * 31 + p, false)
    }

    @Test
    fun theInbredFastPathMatchesTheGeneralScanAtFOneWithTies() {
        for (n in listOf(2, 3, 6, 9)) for (p in listOf(30, 120)) compareInbred(n, p, n * 13 + p, true)
    }

    @Test
    fun atFOneEveryTransitionIntoAHeterozygousStateIsImpossible() {
        // The premise the fast path rests on: no heterozygous state is reachable, from anywhere.
        val calculator = DiploidTransitionProbability(0.999, 1.0, 4)
        for (a in 0 until 4) for (b in 0 until 4) for (c in 0 until 4) for (d in 0 until 4) {
            val p = calculator.calculate(Pair(a, b), Pair(c, d))
            if (c != d) assertEquals(0.0, p, "(($a,$b) -> ($c,$d)) must be impossible at F = 1")
        }
    }

    @Test
    fun handlesASingleCandidateParent() {
        val hmm = ViterbiHMM(0.0, 0.9999, 0.98)
        val emissions = Array(10) { DoubleArray(1) { -1.0 } }
        val result = hmm.viterbiOptimizedForDiploid(1, 10, DoubleArray(1), { p -> emissions[p] })
        assertEquals(List(10) { 0 }, result.first.toList())
    }
}
