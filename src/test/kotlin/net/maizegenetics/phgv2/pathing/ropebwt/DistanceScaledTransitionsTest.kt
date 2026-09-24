package net.maizegenetics.phgv2.pathing.ropebwt

import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test
import kotlin.math.abs
import kotlin.math.exp
import kotlin.math.ln

class DistanceScaledTransitionsTest {

    private val oneMb = 1_000_000.0
    private val defaultProbSwitch = 1e-4

    @Test
    fun stayingIsThePoissonSurvivalOverTheStepDistance() {
        val distances = intArrayOf(0, 1_000, 10_000, 1_000_000)
        val (lnNoSwitch, _) = transitionLogsByDistance(distances, defaultProbSwitch, oneMb, 10)
        val rate = -ln(1.0 - defaultProbSwitch) / oneMb
        for (index in 1 until distances.size) {
            assertEquals(-rate * distances[index], lnNoSwitch[index], 1e-15,
                "step of ${distances[index]}")
        }
        // and at exactly the reference distance, the probability of switching is probSwitch itself
        assertEquals(1.0 - defaultProbSwitch, exp(lnNoSwitch[3]), 1e-12,
            "a 1 Mb step must reproduce --prob-switch exactly")
    }

    @Test
    fun aTenfoldLongerStepCostsTenTimesAsMuch() {
        // The property that makes this worth doing: the penalty tracks sequence, not marker count.
        val (lnNoSwitch, _) = transitionLogsByDistance(intArrayOf(0, 1_000, 10_000), defaultProbSwitch, oneMb, 10)
        assertEquals(10.0, lnNoSwitch[2] / lnNoSwitch[1], 1e-9)
    }

    @Test
    fun switchingIsSharedAmongTheAlternatives() {
        val distances = intArrayOf(0, 500_000)
        for (alternatives in listOf(2, 5, 25)) {
            val (lnNoSwitch, lnSwitch) = transitionLogsByDistance(distances, defaultProbSwitch, oneMb, alternatives)
            val stay = exp(lnNoSwitch[1])
            val perAlternative = exp(lnSwitch[1])
            assertEquals(1.0, stay + perAlternative * (alternatives - 1), 1e-12,
                "staying plus every way of switching must be a distribution, alternatives=$alternatives")
        }
    }

    @Test
    fun withOneAlternativeSwitchingIsImpossibleRatherThanFree() {
        // The divisor would be zero. Left unguarded the quotient becomes infinity and a switch looks
        // infinitely attractive, which is the bug this mirrors in viterbiOptimizedForHaploid.
        val (_, lnSwitch) = transitionLogsByDistance(intArrayOf(0, 1_000), defaultProbSwitch, oneMb, 1)
        assertTrue(lnSwitch[1] < -1e5, "switching must be floored, was ${lnSwitch[1]}")
        assertTrue(lnSwitch[1].isFinite(), "and finite, so it cannot poison a sum")
    }

    @Test
    fun aZeroLengthStepIsNotFreeToSwitchAcross() {
        // Two VCF records at the same position give a step of zero. Charging nothing would let the
        // path change founders for free wherever a site is duplicated.
        val (lnNoSwitch, _) = transitionLogsByDistance(intArrayOf(0, 0, 5), defaultProbSwitch, oneMb, 10)
        assertTrue(lnNoSwitch[1] < 0.0, "a zero step still costs something: ${lnNoSwitch[1]}")
    }

    @Test
    fun thePrecisionHoldsAtTheTinyExponentsThisNormallyProduces() {
        // At 1e-4 per Mb a 1 kb step gives an exponent near -1e-7, where computing 1 - exp(x) by
        // subtraction loses most of its significant digits. expm1 is why this is accurate.
        val (lnNoSwitch, lnSwitch) = transitionLogsByDistance(intArrayOf(0, 1_000), 1e-4, oneMb, 2)
        val naive = ln(1.0 - exp(lnNoSwitch[1]))
        val exact = ln(-kotlin.math.expm1(lnNoSwitch[1]))
        assertEquals(exact, lnSwitch[1], 1e-12)
        // the naive form is accurate enough here to be within a hair, but the exact form is what is
        // used; assert the implementation matches the stable expression rather than the unstable one
        assertTrue(abs(naive - lnSwitch[1]) < 1e-6, "sanity: the two agree to a hair at this scale")
    }

    @Test
    fun elementZeroIsFilledAndNeverNegativeInfinity() {
        // Position 0 has no predecessor so the recursions never read it, but leaving it as a raw 0.0
        // or an infinity invites trouble if something else ever does.
        val (lnNoSwitch, lnSwitch) = transitionLogsByDistance(intArrayOf(0, 1_000, 2_000), defaultProbSwitch, oneMb, 4)
        assertEquals(lnNoSwitch[1], lnNoSwitch[0], 1e-15)
        assertEquals(lnSwitch[1], lnSwitch[0], 1e-15)
        assertTrue(lnSwitch[0].isFinite())
    }

    @Test
    fun badParametersAreRejected() {
        assertThrows(IllegalArgumentException::class.java) {
            transitionLogsByDistance(intArrayOf(), defaultProbSwitch, oneMb, 4)
        }
        for (bad in listOf(0.0, 1.0, -0.1, 1.5)) {
            assertThrows(IllegalArgumentException::class.java, {
                transitionLogsByDistance(intArrayOf(0, 100), bad, oneMb, 4)
            }, "prob-switch $bad must be rejected")
        }
        assertThrows(IllegalArgumentException::class.java) {
            transitionLogsByDistance(intArrayOf(0, 100), defaultProbSwitch, 0.0, 4)
        }
    }

    @Test
    fun stepDistancesAreTheGapsBetweenPositions() {
        val distances = stepDistances(intArrayOf(100, 150, 1_150, 1_150))
        assertArrayEquals(intArrayOf(0, 50, 1_000, 0), distances)
    }

    @Test
    fun nonAscendingPositionsAreRejected() {
        // A negative gap would become a shorter-than-zero step and, unguarded, a cheaper switch.
        val error = assertThrows(IllegalArgumentException::class.java) {
            stepDistances(intArrayOf(100, 50))
        }
        assertTrue(error.message!!.contains("must ascend"), error.message)
    }

    @Test
    fun expectedSwitchesMakesTheParameterInterpretable() {
        // The default, read against a maize genome: 1e-4 per Mb over 2.13 Gb.
        val maize = 2.13e9
        assertEquals(0.213, expectedSwitches(maize, 1e-4, oneMb), 1e-3)
        // and the read-based tuned value, for comparison -- three orders of magnitude smoother
        val ps4gEquivalent = probabilityOfNoSwitch(256.0, 1.0 - 0.999999999, 256.0)
        assertEquals(0.999999999, ps4gEquivalent, 1e-12,
            "prob-switch 1e-9 per 256 bp is the complement of prob-same 0.999999999")
        assertEquals(0.0083, expectedSwitches(maize, 1e-9, 256.0), 1e-4,
            "which is 0.008 switches per maize genome, against roughly 11 in biology")
    }

    @Test
    fun theRecursionsHonourPerPositionTransitions() {
        // The behavioural point: the same emissions and the same number of markers, but a switch is
        // affordable across a long gap and not across a short one.
        //
        // The arithmetic has to be set deliberately, because a switch is only taken when the emission
        // evidence beats its cost. With three sites of advantage E after the boundary, switching wins
        // iff 3E > lnStay - lnSwitch:
        //
        //   10 Mb step at 1e-4 per Mb   cost  6.91 nats  -> 3E = 9 beats it
        //   10 bp step at 1e-12 per Mb  cost 39.1 nats   -> 3E = 9 does not
        //
        // An earlier version of this test used E = 1, where 3 nats loses to both, and the path stayed
        // put in each case for a reason that had nothing to do with the distances.
        val hmm = ViterbiHMM(1.0, 0.99, 0.98)
        val nStates = 2
        val positions = 6
        val advantage = 3.0
        // sites 0-2 favour state 0, sites 3-5 favour state 1
        val emissions = Array(positions) { site ->
            if (site < 3) doubleArrayOf(0.0, -advantage) else doubleArrayOf(-advantage, 0.0)
        }
        val initial = doubleArrayOf(ln(0.5), ln(0.5))

        // Distances: the step into site 3 is enormous, so switching there is cheap.
        val longGapAtBoundary = intArrayOf(0, 10, 10, 10_000_000, 10, 10)
        val (noSwitchA, switchA) = transitionLogsByDistance(longGapAtBoundary, 1e-4, 1e6, nStates)
        val withGap = hmm.viterbiOptimizedForHaploid(nStates, positions, initial,
            { p -> emissions[p] }, noSwitchA, switchA)
        assertEquals(listOf(0, 0, 0, 1, 1, 1), withGap.first.toList(),
            "a long gap at the boundary lets the path follow the emissions")

        // Same emissions, but every step tiny: switching anywhere is expensive, and with only six
        // sites of weak evidence the path should refuse to move.
        val allShort = IntArray(positions) { if (it == 0) 0 else 10 }
        val (noSwitchB, switchB) = transitionLogsByDistance(allShort, 1e-12, 1e6, nStates)
        val withoutGap = hmm.viterbiOptimizedForHaploid(nStates, positions, initial,
            { p -> emissions[p] }, noSwitchB, switchB)
        assertTrue(withoutGap.first.toSet().size == 1,
            "with a severe per-base penalty the same evidence cannot move the path: " +
                    "${withoutGap.first.toList()}")
    }

    @Test
    fun omittingPerPositionTransitionsReproducesTheFixedBehaviour() {
        // Backward compatibility: null must give exactly what the fixed constants gave before, in
        // both recursions, so every existing caller is untouched.
        val hmm = ViterbiHMM(0.0, 0.9999, 0.98)
        val nParents = 3
        val positions = 20
        val random = kotlin.random.Random(4242)
        val emissions = Array(positions) { DoubleArray(nParents * nParents) { random.nextDouble(-4.0, 0.0) } }
        val initial = DoubleArray(nParents * nParents) { -1.0 }

        val withNulls = hmm.viterbiOptimizedForDiploid(nParents, positions, initial, { p -> emissions[p] })
        val explicitFixed = hmm.viterbiOptimizedForDiploid(
            nParents, positions, initial, { p -> emissions[p] },
            DoubleArray(positions) { ln(0.9999) },
            DoubleArray(positions) { ln((1.0 - 0.9999) / (nParents - 1)) }
        )
        assertEquals(withNulls.first.toList(), explicitFixed.first.toList())
        assertEquals(withNulls.second, explicitFixed.second, 1e-9)
    }

    @Test
    fun perPositionTransitionsAreRefusedWhereTheGeneralScanWouldBeNeeded() {
        // At an intermediate coefficient the transition is a full matrix, which would have to be
        // rebuilt per position. Refusing is better than quietly using one step's transition
        // everywhere.
        val hmm = ViterbiHMM(0.5, 0.9999, 0.98)
        val nParents = 2
        val emissions = Array(5) { DoubleArray(4) { -1.0 } }
        val error = assertThrows(IllegalArgumentException::class.java) {
            hmm.findDiploidStatePath(nParents, 5, { p -> emissions[p] }, null,
                DoubleArray(5) { -0.1 }, DoubleArray(5) { -5.0 })
        }
        assertTrue(error.message!!.contains("inbreeding coefficient of 0 or 1"), error.message)
        // and it still works without them
        assertEquals(5, hmm.findDiploidStatePath(nParents, 5, { p -> emissions[p] }).first.size)
    }

    @Test
    fun theSeamPassesPerPositionTransitionsToBothFastPaths() {
        val nParents = 3
        val positions = 12
        val random = kotlin.random.Random(99)
        val emissions = Array(positions) { DoubleArray(nParents * nParents) { random.nextDouble(-4.0, 0.0) } }
        val distances = IntArray(positions) { if (it == 0) 0 else 1_000 * (it + 1) }

        for (f in listOf(0.0, 1.0)) {
            val hmm = ViterbiHMM(f, 0.9999, 0.98)
            val alternatives = if (f == 1.0) nParents else nParents
            val (noSwitch, switch) = transitionLogsByDistance(distances, 1e-3, 1e6, alternatives)
            val scaled = hmm.findDiploidStatePath(nParents, positions, { p -> emissions[p] }, null,
                noSwitch, switch)
            val fixed = hmm.findDiploidStatePath(nParents, positions, { p -> emissions[p] })
            assertEquals(positions, scaled.first.size, "F=$f")
            assertTrue(scaled.second.isFinite(), "F=$f")
            // The two need not agree -- that is the point -- but both must be valid state indices.
            assertTrue(scaled.first.all { it in 0 until nParents * nParents }, "F=$f")
            assertTrue(fixed.first.all { it in 0 until nParents * nParents }, "F=$f")
        }
    }
}
