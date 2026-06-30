package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet

class BackwardForwardHMM(
    private val readMap: Map<Int, MutableList<Ps4gGameteSet>>,
    private val parents: Set<Int>,
    private val pCorrect: Double,
    private val pSameGamete: Double
) {
    /**
     * Hidden Markov Model Forward-Backward Algorithm.
     */

    /**
     * Executes the Forward-Backward algorithm for a given sequence of observations.
     * @param observations List of observation indices corresponding to columns in matrix B.
     * @return Double matrix of size [T x N] containing the posterior probabilities (gamma).
     */
    fun haploidForwardBackward(): Array<DoubleArray> {
        val positionList = readMap.keys.sorted()
        val nPositions = positionList.size
        val nStates = parents.size
        val emissionProbability = HaploidPS4GEmissionProbability(readMap, parents, pCorrect)
        val pNotSameGamete = 1.0 - pSameGamete
        val transitionMatrix = Array(nStates) { DoubleArray(nStates) {pNotSameGamete} }
        for (ndx in 0 until nStates) transitionMatrix[ndx][ndx] = pSameGamete

        //make sure parent indices are sorted  0 to nParents - 1
        val parentIndices = parents.sorted()
        check(parentIndices.last() == nStates - 1) {"gamete indices are not 0..(nGametes - 1) as expected by haploid forward-backward"}

        // 1. Forward Pass (Compute Alpha)
        val alpha = Array(nPositions) { DoubleArray(nStates) }

        // Initialization
        for (i in 0 until nStates) {
            alpha[0][i] = 0.0 + emissionProbability.getLnProbObsGivenState(i, positionList[0])
        }

        // Induction
        for (pos in 1 until nPositions) {
            for (parentIndex in parentIndices) {
                var sum = 0.0
                for (i in 0 until nStates) {
                    sum += alpha[pos - 1][i] * transitionMatrix[i][parentIndex]
                }
                alpha[pos][parentIndex] = sum * emissionProbability.getProbObsGivenState(parentIndex,pos)
            }
        }

        // 2. Backward Pass (Compute Beta)
        val beta = Array(nPositions) { DoubleArray(nStates) }

        // Initialization
        for (i in 0 until nStates) {
            beta[nPositions - 1][i] = 1.0
        }

        // Induction
        for (t in nPositions - 2 downTo 0) {
            for (i in 0 until nStates) {
                var sum = 0.0
                for (j in 0 until nStates) {
                    sum += transitionMatrix[i][j] * emissionProbability.getProbObsGivenState(i,t) * beta[t + 1][j]
                }
                beta[t][i] = sum
            }
        }

        // 3. Compute Posterior Probabilities (Gamma)
        val gamma = Array(nPositions) { DoubleArray(nStates) }
        for (t in 0 until nPositions) {
            var totalProbAtTimeT = 0.0
            for (i in 0 until nStates) {
                gamma[t][i] = alpha[t][i] * beta[t][i]
                totalProbAtTimeT += gamma[t][i]
            }

            // Normalize to make it a valid probability distribution
            if (totalProbAtTimeT > 0.0) {
                for (i in 0 until nStates) {
                    gamma[t][i] /= totalProbAtTimeT
                }
            }
        }

        return gamma
    }


}

