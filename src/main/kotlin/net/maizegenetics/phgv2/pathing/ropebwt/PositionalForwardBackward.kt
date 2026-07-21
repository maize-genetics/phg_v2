package net.maizegenetics.phgv2.pathing.ropebwt

import kotlin.math.ln
import kotlin.reflect.KFunction1

/**
 * Forward-backward algorithm for a Hidden Markov Model defined over a 1-D series of
 * positions rather than time points. The model structure is identical to the usual
 * time-indexed HMM; only the terminology changes (position instead of time).
 *
 * The transition probability function is assumed to be homogeneous: it is the same
 * function evaluated at every position, taking only the "from" and "to" state indices.
 * The emission probability function depends on both the position and the state, since
 * observations can differ arbitrarily from position to position.
 *
 * Numerical stability:
 * Because the number of positions can be very large (100,000+), the raw forward and
 * backward probabilities would underflow to zero if computed directly. This
 * implementation therefore uses the standard Rabiner scaling approach: at every
 * position the forward vector is renormalized to sum to 1, and the scaling factor
 * used is recorded and reused when computing the backward vector and the total
 * log-likelihood. This keeps every stored value in a numerically safe range while
 * remaining exactly equivalent (up to floating point rounding) to the unscaled
 * computation.
 *
 * Performance notes:
 * - The transition function is evaluated once per (fromState, toState) pair and cached
 *   in a flat DoubleArray, since it does not depend on position. This avoids repeating
 *   potentially expensive lambda calls numberOfPositions times.
 * - The emission function is evaluated once per (position, state) pair as needed; it is
 *   not cached across the whole series because that would require
 *   numberOfPositions * numberOfStates doubles of memory in addition to the forward and
 *   backward tables, which can be prohibitive for very large inputs.
 * - Inner loops operate on primitive DoubleArray buffers only (no boxing, no per-step
 *   allocation) to keep the O(numberOfPositions * numberOfStates^2) work as fast as
 *   possible on the JVM.
 * - Memory usage is O(numberOfPositions * numberOfStates) for the stored forward and
 *   backward tables (and the returned posterior table). For example, 100,000 positions
 *   with 300 states requires roughly 3 arrays of 30,000,000 doubles each (~240 MB each,
 *   ~720 MB total). Plan memory accordingly for very large problems.
 */
class PositionalForwardBackward(
    private val numberOfStates: Int,
    private val numberOfPositions: Int,
    private val initialStateProbabilities: DoubleArray,
    private val transitionMatrix: DoubleArray,
    private val emissionProbability: KFunction1<Int, DoubleArray>
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

        computeForward(forwardTable, scalingFactors)

        val backwardTable = Array(numberOfPositions) { DoubleArray(numberOfStates) }
        computeBackward(backwardTable, scalingFactors)

        val posteriorStateProbabilities = computePosteriors(forwardTable, backwardTable)

        var logLikelihood = 0.0
        for (position in 0 until numberOfPositions) {
            logLikelihood += ln(scalingFactors[position])
        }

        return Result(posteriorStateProbabilities, logLikelihood)
    }

    /**
     * Fills forwardTable[position][state] with the scaled forward probability
     * alpha(position, state) = P(observations up to and including position, state at
     * position | model), normalized to sum to 1 at each position, and records the
     * normalization constant used at each position in scalingFactors.
     */
    private fun computeForward(forwardTable: Array<DoubleArray>, scalingFactors: DoubleArray) {
        // Initialization at the first position.
        //todo: make sure that position is 0 based and continuous
        val firstPositionForward = forwardTable[0]
        var firstPositionTotal = 0.0
        var emissionArray = emissionProbability(0)
        for (state in 0 until numberOfStates) {
            val value = initialStateProbabilities[state] * emissionArray[state]
            firstPositionForward[state] = value
            firstPositionTotal += value
        }
        scalingFactors[0] = firstPositionTotal
        normalizeInPlace(firstPositionForward, firstPositionTotal)

        // Recursion over the remaining positions.
        for (position in 1 until numberOfPositions) {
            val previousForward = forwardTable[position - 1]
            val currentForward = forwardTable[position]
            var currentPositionTotal = 0.0

            emissionArray = emissionProbability(position)
            for (toState in 0 until numberOfStates) {
                var sumOverPreviousStates = 0.0
                for (fromState in 0 until numberOfStates) {
                    val previousValue = previousForward[fromState]
                    if (previousValue != 0.0) {
                        sumOverPreviousStates += previousValue * transitionMatrix[fromState * numberOfStates + toState]
                    }
                }
                val value = sumOverPreviousStates * emissionArray[toState]
                currentForward[toState] = value
                currentPositionTotal += value
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

        // Reusable buffer holding emission(position + 1, state) values, avoids
        // recomputing emissionProbability multiple times per state pair.

        for (position in lastPosition - 1 downTo 0) {
            val nextBackward = backwardTable[position + 1]
            val currentBackward = backwardTable[position]

            val nextPositionEmissions = emissionProbability(position + 1)

            for (fromState in 0 until numberOfStates) {
                var sumOverNextStates = 0.0
                val rowOffset = fromState * numberOfStates
                for (toState in 0 until numberOfStates) {
                    val nextBackwardValue = nextBackward[toState]
                    if (nextBackwardValue != 0.0) {
                        sumOverNextStates += transitionMatrix[rowOffset + toState] *
                            nextPositionEmissions[toState] * nextBackwardValue
                    }
                }
                currentBackward[fromState] = sumOverNextStates
            }

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

