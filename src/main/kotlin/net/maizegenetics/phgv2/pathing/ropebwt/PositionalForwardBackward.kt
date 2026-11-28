package net.maizegenetics.phgv2.pathing.ropebwt

import kotlin.math.exp
import kotlin.math.ln

/**
 * One step of an HMM's transition, applied to a vector rather than read off a stored matrix, so
 * that a transition that varies from position to position, or that has structure, never has to be
 * materialised as a full `numberOfStates^2` matrix at every position.
 */
fun interface TransitionStep {
    /**
     * Multiplies [vector] by the transition for the step from `position - 1` into [position] and
     * writes the result to [out]: when [forward] is true, `out[to] = sum over from of vector[from]
     * T(from, to)`; when false, `out[from] = sum over to of T(from, to) vector[to]`. [out] and
     * [vector] are distinct arrays of length numberOfStates.
     */
    fun apply(position: Int, vector: DoubleArray, out: DoubleArray, forward: Boolean)
}

/**
 * Forward-backward algorithm for a Hidden Markov Model defined over a 1-D series of
 * positions rather than time points. The model structure is identical to the usual
 * time-indexed HMM; only the terminology changes (position instead of time).
 *
 * The transition is supplied as a [TransitionStep], which may differ from position to position,
 * as it does when the switch probability is scaled by the distance between positions. The emission
 * depends on both the position and the state, and is supplied as natural log probabilities, since
 * a product over many reads underflows in probability space.
 *
 * Numerical stability:
 * Because the number of positions can be very large (100,000+), the raw forward and
 * backward probabilities would underflow to zero if computed directly. This
 * implementation therefore uses the standard Rabiner scaling approach: at every
 * position the forward vector is renormalized to sum to 1, and the scaling factor
 * used is recorded and reused when computing the backward vector. The log emissions at
 * each position are shifted by their maximum before exponentiating; the shift is a
 * constant across states, so it cancels in the normalization and is added back only to
 * the log-likelihood. This keeps every stored value in a numerically safe range while
 * remaining exactly equivalent (up to floating point rounding) to the unscaled
 * computation.
 *
 * Performance notes:
 * - The emission function is evaluated once per position in each pass; it is
 *   not cached across the whole series because that would require
 *   numberOfPositions * numberOfStates doubles of memory in addition to the forward and
 *   backward tables, which can be prohibitive for very large inputs.
 * - Memory usage is O(numberOfPositions * numberOfStates) for the stored forward and
 *   backward tables (and the returned posterior table). For example, 100,000 positions
 *   with 300 states requires roughly 3 arrays of 30,000,000 doubles each (~240 MB each,
 *   ~720 MB total). Plan memory accordingly for very large problems.
 *
 * @param initialStateProbabilities the state distribution at the first position, in probability space.
 * @param transition the transition step into each position; never called for position 0.
 * @param lnEmissionProbability natural log emission probabilities at a position, one per state.
 */
class PositionalForwardBackward(
    private val numberOfStates: Int,
    private val numberOfPositions: Int,
    private val initialStateProbabilities: DoubleArray,
    private val transition: TransitionStep,
    private val lnEmissionProbability: (Int) -> DoubleArray,
) {

    init {
        require(numberOfStates > 0) { "numberOfStates must be positive" }
        require(numberOfPositions > 0) { "numberOfPositions must be positive" }
        require(initialStateProbabilities.size == numberOfStates) {
            "initialStateProbabilities must have length numberOfStates"
        }
    }

    /**
     * Result of running the forward-backward algorithm.
     *
     * @param posteriorStateProbabilities posteriorStateProbabilities[position][state] is
     *   P(state at this position | entire observed series), i.e. the smoothed marginal
     *   posterior probability.
     * @param logLikelihood the total log-likelihood of the observed series under the model,
     *   log P(observations).
     */
    class Result(
        val posteriorStateProbabilities: Array<DoubleArray>,
        val logLikelihood: Double
    )

    /**
     * Runs the full forward-backward algorithm over all positions and returns the
     * per-position posterior state probabilities together with the overall
     * log-likelihood of the series.
     */
    fun run(): Result {
        val forwardTable = Array(numberOfPositions) { DoubleArray(numberOfStates) }
        val scalingFactors = DoubleArray(numberOfPositions)
        val emissionShifts = DoubleArray(numberOfPositions)

        computeForward(forwardTable, scalingFactors, emissionShifts)

        val backwardTable = Array(numberOfPositions) { DoubleArray(numberOfStates) }
        computeBackward(backwardTable, scalingFactors)

        val posteriorStateProbabilities = computePosteriors(forwardTable, backwardTable)

        var logLikelihood = 0.0
        for (position in 0 until numberOfPositions) {
            logLikelihood += ln(scalingFactors[position]) + emissionShifts[position]
        }

        return Result(posteriorStateProbabilities, logLikelihood)
    }

    /**
     * Fills [emissions] with the emission probabilities at [position], each divided by the largest
     * of them, and returns the log of that divisor.
     */
    private fun scaledEmissions(position: Int, emissions: DoubleArray): Double {
        val lnEmissions = lnEmissionProbability(position)
        require(lnEmissions.size == numberOfStates) { "emission array must have length numberOfStates" }
        val shift = lnEmissions.max()
        require(shift > Double.NEGATIVE_INFINITY) { "every state has zero emission probability at position $position" }
        for (state in 0 until numberOfStates) emissions[state] = exp(lnEmissions[state] - shift)
        return shift
    }

    /**
     * Fills forwardTable[position][state] with the scaled forward probability
     * alpha(position, state) = P(observations up to and including position, state at
     * position | model), normalized to sum to 1 at each position, and records the
     * normalization constant and emission shift used at each position.
     */
    private fun computeForward(forwardTable: Array<DoubleArray>, scalingFactors: DoubleArray,
                               emissionShifts: DoubleArray) {
        val emissions = DoubleArray(numberOfStates)

        // Initialization at the first position.
        val firstPositionForward = forwardTable[0]
        var firstPositionTotal = 0.0
        emissionShifts[0] = scaledEmissions(0, emissions)
        for (state in 0 until numberOfStates) {
            val value = initialStateProbabilities[state] * emissions[state]
            firstPositionForward[state] = value
            firstPositionTotal += value
        }
        scalingFactors[0] = firstPositionTotal
        normalizeInPlace(firstPositionForward, firstPositionTotal)

        // Recursion over the remaining positions.
        for (position in 1 until numberOfPositions) {
            val currentForward = forwardTable[position]
            transition.apply(position, forwardTable[position - 1], currentForward, forward = true)

            emissionShifts[position] = scaledEmissions(position, emissions)
            var currentPositionTotal = 0.0
            for (state in 0 until numberOfStates) {
                currentForward[state] *= emissions[state]
                currentPositionTotal += currentForward[state]
            }

            scalingFactors[position] = currentPositionTotal
            normalizeInPlace(currentForward, currentPositionTotal)
        }
    }

    /**
     * Fills backwardTable[position][state] with the scaled backward probability
     * beta(position, state) = P(observations after position | state at position, model),
     * using the same scaling factors computed during the forward pass so that
     * forward and backward values remain on a consistent scale.
     */
    private fun computeBackward(backwardTable: Array<DoubleArray>, scalingFactors: DoubleArray) {
        val lastPosition = numberOfPositions - 1
        val lastPositionBackward = backwardTable[lastPosition]
        for (state in 0 until numberOfStates) {
            lastPositionBackward[state] = 1.0
        }
        normalizeInPlace(lastPositionBackward, scalingFactors[lastPosition])

        // emission(position + 1, state) * beta(position + 1, state), reused at every position
        val weightedNext = DoubleArray(numberOfStates)
        for (position in lastPosition - 1 downTo 0) {
            val nextBackward = backwardTable[position + 1]
            scaledEmissions(position + 1, weightedNext)
            for (state in 0 until numberOfStates) weightedNext[state] *= nextBackward[state]

            val currentBackward = backwardTable[position]
            transition.apply(position + 1, weightedNext, currentBackward, forward = false)
            normalizeInPlace(currentBackward, scalingFactors[position])
        }
    }

    /**
     * Combines the forward and backward tables into per-position posterior state
     * probabilities, renormalizing explicitly at each position as a safeguard against
     * floating point drift.
     */
    private fun computePosteriors(
        forwardTable: Array<DoubleArray>,
        backwardTable: Array<DoubleArray>
    ): Array<DoubleArray> {
        val posteriorStateProbabilities = Array(numberOfPositions) { DoubleArray(numberOfStates) }

        for (position in 0 until numberOfPositions) {
            val forwardAtPosition = forwardTable[position]
            val backwardAtPosition = backwardTable[position]
            val posteriorAtPosition = posteriorStateProbabilities[position]

            var total = 0.0
            for (state in 0 until numberOfStates) {
                val value = forwardAtPosition[state] * backwardAtPosition[state]
                posteriorAtPosition[state] = value
                total += value
            }
            normalizeInPlace(posteriorAtPosition, total)
        }

        return posteriorStateProbabilities
    }

    private fun normalizeInPlace(values: DoubleArray, total: Double) {
        require(total > 0.0) {
            "Encountered a zero-probability position; check the transition and emission functions."
        }
        val inverseTotal = 1.0 / total
        for (index in values.indices) {
            values[index] *= inverseTotal
        }
    }
}
