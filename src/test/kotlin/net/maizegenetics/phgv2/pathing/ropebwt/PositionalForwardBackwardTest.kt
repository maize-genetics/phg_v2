package net.maizegenetics.phgv2.pathing.ropebwt

import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test
import kotlin.math.exp
import kotlin.math.ln

class PositionalForwardBackwardTest {

    /** Dense matrix of [step] at [position], from applying it forward to each unit vector. */
    private fun denseMatrix(step: TransitionStep, position: Int, numberOfStates: Int): Array<DoubleArray> {
        val matrix = Array(numberOfStates) { DoubleArray(numberOfStates) }
        val unit = DoubleArray(numberOfStates)
        val out = DoubleArray(numberOfStates)
        for (from in 0 until numberOfStates) {
            unit.fill(0.0)
            unit[from] = 1.0
            step.apply(position, unit, out, forward = true)
            for (to in 0 until numberOfStates) matrix[from][to] = out[to]
        }
        return matrix
    }

    @Test
    fun testDiploidSwitchStepIsKroneckerProductOfHaploidSteps() {
        val nParents = 3
        val lnNoSwitch = doubleArrayOf(0.0, ln(0.9), ln(0.6))
        val lnSwitch = doubleArrayOf(0.0, ln(0.05), ln(0.2))
        val haploid = SwitchTransitionStep(nParents, false, lnNoSwitch, lnSwitch)
        val diploid = SwitchTransitionStep(nParents, true, lnNoSwitch, lnSwitch)
        val nStates = nParents * nParents

        for (position in 1..2) {
            val t = denseMatrix(haploid, position, nParents)
            for (from in 0 until nParents) {
                assertEquals(1.0, t[from].sum(), 1e-12, "haploid row $from at position $position")
            }
            val vector = DoubleArray(nStates) { 0.1 + it * 0.07 }
            val forwardOut = DoubleArray(nStates)
            val backwardOut = DoubleArray(nStates)
            diploid.apply(position, vector, forwardOut, forward = true)
            diploid.apply(position, vector, backwardOut, forward = false)
            for (c in 0 until nParents) for (d in 0 until nParents) {
                var forwardExpected = 0.0
                var backwardExpected = 0.0
                for (a in 0 until nParents) for (b in 0 until nParents) {
                    forwardExpected += vector[a * nParents + b] * t[a][c] * t[b][d]
                    backwardExpected += t[c][a] * t[d][b] * vector[a * nParents + b]
                }
                assertEquals(forwardExpected, forwardOut[c * nParents + d], 1e-12)
                assertEquals(backwardExpected, backwardOut[c * nParents + d], 1e-12)
            }
        }
    }

    @Test
    fun testPosteriorsAndLikelihoodMatchEnumerationOverAllPaths() {
        // Two parents, diploid: 4 states over 5 positions, with transitions that differ by position and
        // log emissions far below zero, which would underflow if exponentiated without the shift.
        val nParents = 2
        val nStates = 4
        val nPositions = 5
        val lnNoSwitch = DoubleArray(nPositions) { ln(0.95 - 0.1 * it) }
        val lnSwitch = DoubleArray(nPositions) { ln(0.05 + 0.1 * it) }
        val step = SwitchTransitionStep(nParents, true, lnNoSwitch, lnSwitch)
        val lnEmissions = Array(nPositions) { position ->
            DoubleArray(nStates) { state -> -900.0 - 3.0 * ((state + position) % nStates) + 0.5 * state }
        }
        val initial = doubleArrayOf(0.1, 0.2, 0.3, 0.4)

        val result = PositionalForwardBackward(nStates, nPositions, initial, step) { lnEmissions[it].copyOf() }.run()

        // Every path enumerated, with emissions shifted by a constant that is added back to the likelihood.
        val shift = -900.0
        val matrices = (0 until nPositions).map { if (it == 0) null else denseMatrix(step, it, nStates) }
        val marginals = Array(nPositions) { DoubleArray(nStates) }
        var total = 0.0
        val path = IntArray(nPositions)
        var pathCount = 1
        repeat(nPositions) { pathCount *= nStates }
        for (code in 0 until pathCount) {
            var rest = code
            for (position in 0 until nPositions) {
                path[position] = rest % nStates
                rest /= nStates
            }
            var probability = initial[path[0]] * exp(lnEmissions[0][path[0]] - shift)
            for (position in 1 until nPositions) {
                probability *= matrices[position]!![path[position - 1]][path[position]] *
                        exp(lnEmissions[position][path[position]] - shift)
            }
            total += probability
            for (position in 0 until nPositions) marginals[position][path[position]] += probability
        }

        assertEquals(ln(total) + nPositions * shift, result.logLikelihood, 1e-9)
        for (position in 0 until nPositions) for (state in 0 until nStates) {
            assertEquals(marginals[position][state] / total, result.posteriorStateProbabilities[position][state], 1e-12,
                "posterior at position $position, state $state")
        }
    }
}
