package net.maizegenetics.phgv2.pathing.ropebwt

import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test
import kotlin.math.exp
import kotlin.math.pow
import net.maizegenetics.phgv2.pathing.ropebwt.SelfingTransitions.Companion.AA
import net.maizegenetics.phgv2.pathing.ropebwt.SelfingTransitions.Companion.AB
import net.maizegenetics.phgv2.pathing.ropebwt.SelfingTransitions.Companion.BA
import net.maizegenetics.phgv2.pathing.ropebwt.SelfingTransitions.Companion.BB

class SelfingTransitionsTest {

    /** First-order rate per unit r, read off the exact transition at a small r. */
    private fun rate(model: SelfingTransitions, from: Int, to: Int, r: Double = 1e-9) =
        model.biparentalTransitions(r)[from * 4 + to] / r

    @Test
    fun theFirstOrderRatesMatchTheClosedForms() {
        for (n in 1..8) {
            val model = SelfingTransitions(n, 2)
            val tag = "N=$n"
            val twoN = 2.0.pow(n)
            assertEquals((2 * twoN - 2 * n - 2) / (twoN - 1), rate(model, AA, BB), 1e-6, tag)
            assertEquals(n / (twoN - 1), rate(model, AA, AB), 1e-6, tag)
            assertEquals(n / (twoN - 1), rate(model, AA, BA), 1e-6, tag)
            assertEquals(n.toDouble(), rate(model, AB, AA), 1e-6, tag)
            assertEquals(n.toDouble(), rate(model, AB, BB), 1e-6, tag)
            assertEquals(0.0, rate(model, AB, BA), 1e-6, "$tag: a phase swap needs two recombinations")
            assertEquals(model.biparentalTransitions(1e-9)[BB * 4 + AA], model.biparentalTransitions(1e-9)[AA * 4 + BB],
                1e-20, "$tag: the founders are symmetric")
        }
    }

    @Test
    fun breakpointsOverHalfAMorganMatchTheSelfingDocument() {
        // Expected breakpoints over a length L at small r per unit: the stationary flow out of each
        // class. "Crossovers Under Selfing" gives these per 50 cM arm.
        val expected = mapOf(3 to Pair(0.5, 0.75), 4 to Pair(11.0 / 16, 0.5), 5 to Pair(0.8125, 0.3125))
        for ((n, homHet) in expected) {
            val model = SelfingTransitions(n, 2)
            val h = model.heterozygosity
            val homozygous = (1 - h) * rate(model, AA, BB) * 0.5
            val heterozygous = ((1 - h) * 2 * rate(model, AA, AB) + h * 2 * rate(model, AB, AA)) * 0.5
            assertEquals(homHet.first, homozygous, 1e-4, "homozygous breakpoints, N=$n")
            assertEquals(homHet.second, heterozygous, 1e-4, "heterozygous breakpoints, N=$n")
        }
    }

    @Test
    fun everyRowIsADistribution() {
        for (n in listOf(1, 3, 5, 12)) for (r in listOf(1e-22, 1e-6, 0.01, 0.3, 0.5)) {
            val p = SelfingTransitions(n, 2).biparentalTransitions(r)
            for (from in 0 until 4) {
                assertEquals(1.0, (0 until 4).sumOf { p[from * 4 + it] }, 1e-12, "N=$n r=$r from=$from")
                assertTrue((0 until 4).all { p[from * 4 + it] >= 0.0 }, "N=$n r=$r from=$from")
            }
        }
        for (k in 1..4) {
            val model = SelfingTransitions(4, k)
            val matrix = model.logMatrix(0.02)
            val nStates = k * k
            for (from in 0 until nStates) {
                assertEquals(1.0, (0 until nStates).sumOf { exp(matrix[from * nStates + it]) }, 1e-12,
                    "K=$k from=$from")
            }
        }
    }

    @Test
    fun theInitialProbabilitiesAreTheStationaryDistribution() {
        for (k in 2..4) for (n in listOf(1, 4, 5)) {
            val model = SelfingTransitions(n, k)
            val nStates = k * k
            val initial = model.initialLogProbabilities().map { exp(it) }
            assertEquals(1.0, initial.sum(), 1e-12, "K=$k N=$n")
            val matrix = model.logMatrix(0.05)
            for (to in 0 until nStates) {
                val flow = (0 until nStates).sumOf { initial[it] * exp(matrix[it * nStates + to]) }
                assertEquals(initial[to], flow, 1e-12, "K=$k N=$n state=$to")
            }
        }
    }

    @Test
    fun manyGenerationsApproachTheRecombinantInbredLine() {
        // Haldane and Waddington: a RIL by selfing has R = 2r / (1 + 2r) between two loci.
        val model = SelfingTransitions(40, 2)
        for (r in listOf(0.01, 0.1, 0.3)) {
            assertEquals(2 * r / (1 + 2 * r), model.biparentalTransitions(r)[AA * 4 + BB], 1e-9, "r=$r")
        }
        assertEquals(0.0, model.biparentalTransitions(0.1)[AA * 4 + AB], 1e-9)
    }

    @Test
    fun unlinkedLociAreIndependent() {
        val model = SelfingTransitions(5, 2)
        val p = model.biparentalTransitions(0.5)
        val h = model.heterozygosity
        for (from in 0 until 4) {
            assertEquals((1 - h) / 2, p[from * 4 + AA], 1e-12)
            assertEquals(h / 2, p[from * 4 + AB], 1e-12)
        }
    }

    @Test
    fun theTinyRatesATunedProbSwitchGivesKeepTheirPrecision() {
        // At the default --prob-switch, 1e-16 per Mb, a 256 bp step is r of about 2.6e-20.
        val model = SelfingTransitions(5, 2)
        val r = SelfingTransitions.recombinationFraction(256.0, 1e-16, 1e6)
        assertEquals(2.56e-20, r, 1e-30)
        val logs = model.biparentalLogTransitions(r)
        assertEquals(kotlin.math.ln(rate(model, AA, BB) * r), logs[AA * 4 + BB], 1e-6)
        val leaving = rate(model, AA, BB) + rate(model, AA, AB) + rate(model, AA, BA)
        assertEquals(-leaving * r, logs[AA * 4 + AA], 1e-6 * leaving * r,
            "staying must be ln1p of the leaving total, not ln of a sum that rounds to 1")
        assertTrue(logs[AA * 4 + AA] < 0.0)
    }

    @Test
    fun theRecombinationFractionIsTheSwitchProbabilityWhenSmallAndNeverAboveOneHalf() {
        assertEquals(1e-4, SelfingTransitions.recombinationFraction(1e6, 1e-4, 1e6), 1e-8)
        assertEquals(0.5, SelfingTransitions.recombinationFraction(1e9, 0.9, 1e6), 1e-12)
    }

    @Test
    fun aThirdFounderTakesAShareOfEachChange() {
        val k = 3
        val model = SelfingTransitions(4, k)
        val bip = model.biparentalLogTransitions(0.01)
        val matrix = model.logMatrix(0.01)
        fun state(a: Int, b: Int) = a * k + b
        fun cell(from: Int, to: Int) = matrix[from * k * k + to]
        val half = kotlin.math.ln(2.0)
        assertEquals(bip[AA * 4 + BB] - half, cell(state(0, 0), state(1, 1)), 1e-12)
        assertEquals(bip[AA * 4 + BB] - half, cell(state(0, 0), state(2, 2)), 1e-12)
        assertEquals(bip[AA * 4 + AB] - half, cell(state(0, 0), state(0, 2)), 1e-12)
        assertEquals(bip[AB * 4 + BB], cell(state(0, 1), state(1, 1)), 1e-12)
        assertEquals(bip[AB * 4 + AA], cell(state(0, 1), state(0, 0)), 1e-12)
        assertEquals(bip[AB * 4 + BA], cell(state(0, 1), state(1, 0)), 1e-12)
        assertTrue(cell(state(0, 0), state(1, 2)) < -1e5, "a homozygote cannot gain two new founders")
        assertTrue(cell(state(0, 1), state(2, 1)) < -1e5, "a heterozygote stays within its own pair")
        assertTrue(cell(state(0, 1), state(2, 2)) < -1e5, "a heterozygote stays within its own pair")
        assertEquals(bip[AA * 4 + AA], cell(state(2, 2), state(2, 2)), 1e-15)
    }

    @Test
    fun badParametersAreRejected() {
        assertThrows(IllegalArgumentException::class.java) { SelfingTransitions(0, 2) }
        assertThrows(IllegalArgumentException::class.java) { SelfingTransitions(5, SelfingTransitions.MAX_PARENTS + 1) }
        assertThrows(IllegalArgumentException::class.java) { SelfingTransitions(5, 2).biparentalTransitions(0.6) }
    }

    @Test
    fun aHomozygousBreakpointIsCalledDirectlyNotThroughAHeterozygousSegment() {
        // Weak evidence -- one read's worth per position -- for (A, A) then (B, B). The F = 0 model
        // pays two switches for (A, A) -> (B, B) whichever way it goes, so a heterozygous run between
        // the two costs it nothing more in transitions. The selfing model makes the direct step about
        // as cheap as one switch and a heterozygous segment expensive.
        val nParents = 2
        val positions = 60
        val distances = IntArray(positions) { if (it == 0) 0 else 256 }
        val emissions = Array(positions) { position ->
            val towardA = position < 30
            DoubleArray(4) { state ->
                when (state) {
                    AA -> if (towardA) -0.1 else -2.5
                    BB -> if (towardA) -2.5 else -0.1
                    else -> -0.8
                }
            }
        }
        val model = SelfingTransitions(5, nParents)
        val hmm = ViterbiHMM(0.0, 1.0 - 1e-6, 0.98)
        val selfed = hmm.findDiploidStatePath(nParents, positions, { emissions[it] },
            transitionLogsByPosition = model.logMatricesByDistance(distances, 1e-6, 1e6),
            initialLogProbabilities = model.initialLogProbabilities()).first
        assertEquals(List(30) { AA } + List(30) { BB }, selfed.toList())
    }

    @Test
    fun thePerPositionScanMatchesAFixedMatrixWhenEveryStepIsTheSame() {
        val nParents = 3
        val positions = 25
        val random = kotlin.random.Random(7)
        val emissions = Array(positions) { DoubleArray(nParents * nParents) { random.nextDouble(-4.0, 0.0) } }
        val model = SelfingTransitions(4, nParents)
        val distances = IntArray(positions) { 10_000 }
        val r = SelfingTransitions.recombinationFraction(10_000.0, 1e-2, 1e6)
        val hmm = ViterbiHMM(0.0, 0.999, 0.98)
        val byPosition = hmm.findDiploidStatePath(nParents, positions, { emissions[it] },
            transitionLogsByPosition = model.logMatricesByDistance(distances, 1e-2, 1e6),
            initialLogProbabilities = model.initialLogProbabilities())
        val fixed = hmm.viterbiOptimized(nParents * nParents, positions, model.initialLogProbabilities(),
            model.logMatrix(r), { emissions[it] })
        assertEquals(fixed.first.toList(), byPosition.first.toList())
        assertEquals(fixed.second, byPosition.second, 1e-9)
    }

    @Test
    fun aFullTransitionCannotBeCombinedWithPerHaplotypeSwitches() {
        val hmm = ViterbiHMM(0.0, 0.999, 0.98)
        val model = SelfingTransitions(3, 2)
        assertThrows(IllegalArgumentException::class.java) {
            hmm.findDiploidStatePath(2, 3, { DoubleArray(4) }, null, DoubleArray(3), DoubleArray(3),
                model.logMatricesByDistance(intArrayOf(0, 5, 5), 1e-3, 1e6))
        }
    }

    /** The dense K-founder probability matrix, from the log matrix the Viterbi path uses. */
    private fun denseMatrix(model: SelfingTransitions, r: Double) = model.logMatrix(r).map { exp(it) }.toDoubleArray()

    @Test
    fun theProbabilitySpaceStepMatchesTheLogMatrix() {
        for (k in 1..3) for (n in listOf(1, 5)) for (r in listOf(1e-12, 1e-3, 0.2)) {
            val model = SelfingTransitions(n, k)
            val states = k * k
            val matrix = denseMatrix(model, r)
            val biparental = model.biparentalTransitions(r)
            val vector = DoubleArray(states) { 0.1 + it * 0.37 }
            val forward = DoubleArray(states)
            val backward = DoubleArray(states)
            model.multiply(biparental, vector, forward, forward = true)
            model.multiply(biparental, vector, backward, forward = false)
            for (j in 0 until states) {
                val expectedForward = (0 until states).sumOf { i -> vector[i] * matrix[i * states + j] }
                val expectedBackward = (0 until states).sumOf { i -> matrix[j * states + i] * vector[i] }
                val tag = "K=$k N=$n r=$r state $j"
                assertEquals(expectedForward, forward[j], 1e-12 * expectedForward + 1e-300, tag)
                assertEquals(expectedBackward, backward[j], 1e-12 * expectedBackward + 1e-300, tag)
            }
        }
    }

    @Test
    fun theProbabilityInitialDistributionIsStationaryAndMatchesTheLogs() {
        for (k in 1..4) for (n in listOf(1, 3, 6)) {
            val model = SelfingTransitions(n, k)
            val initial = model.initialProbabilities()
            assertEquals(1.0, initial.sum(), 1e-12, "K=$k N=$n")
            if (k > 1) {
                val logs = model.initialLogProbabilities()
                for (s in initial.indices) assertEquals(exp(logs[s]), initial[s], 1e-15)
            }
            val next = DoubleArray(initial.size)
            model.multiply(model.biparentalTransitions(0.05), initial, next, forward = true)
            for (s in initial.indices) assertEquals(initial[s], next[s], 1e-12, "K=$k N=$n state $s")
        }
    }

    @Test
    fun forwardBackwardWithTheSelfingStepMatchesEnumeratingEveryPath() {
        // Brute force over all 4^5 state paths of a two-founder line, against the forward-backward
        // posteriors through SelfingTransitionStep: checks both the forward and the backward step.
        val model = SelfingTransitions(3, 2)
        val distances = intArrayOf(0, 200_000, 50_000, 400_000, 100_000)
        val positions = distances.size
        val probSwitch = 0.3
        val referenceDistance = 1_000_000.0
        val lnEmissions = Array(positions) { p -> DoubleArray(4) { s -> -((p * 7 + s * 3) % 5).toDouble() } }
        val initial = model.initialProbabilities()
        val step = SelfingTransitionStep(model, distances, probSwitch, referenceDistance)
        val result = PositionalForwardBackward(4, positions, initial, step, { lnEmissions[it] }).run()

        val matrices = Array(positions) { p ->
            if (p == 0) DoubleArray(0)
            else denseMatrix(model, SelfingTransitions.recombinationFraction(distances[p].toDouble(), probSwitch, referenceDistance))
        }
        val marginal = Array(positions) { DoubleArray(4) }
        var total = 0.0
        val path = IntArray(positions)
        for (code in 0 until (1 shl (2 * positions))) {
            for (p in 0 until positions) path[p] = (code shr (2 * p)) and 3
            var weight = initial[path[0]] * exp(lnEmissions[0][path[0]])
            for (p in 1 until positions) weight *= matrices[p][path[p - 1] * 4 + path[p]] * exp(lnEmissions[p][path[p]])
            total += weight
            for (p in 0 until positions) marginal[p][path[p]] += weight
        }
        assertEquals(kotlin.math.ln(total), result.logLikelihood, 1e-9)
        for (p in 0 until positions) for (s in 0 until 4) {
            assertEquals(marginal[p][s] / total, result.posteriorStateProbabilities[p][s], 1e-10, "position $p state $s")
        }
    }
}
