package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.PathFinderHMMPS4G
import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import net.maizegenetics.phgv2.utils.Position
import org.apache.logging.log4j.LogManager
import kotlin.math.ln

class ViterbiHMM(val inbreedingCoefficient: Double, val sameGameteProbability: Double, val probCorrect: Double
) {
    private val myLogger = LogManager.getLogger(ViterbiHMM::class.java)

    fun findHaploidPath(contig: String,
                        gameteIndexMap: Map<Int,String>,
                        readMap: Map<Int, MutableList<Ps4gGameteSet>>,
                        likelyParentSet: Set<Int>): List<Pair<Position, String>>  {

        val nParents = likelyParentSet.size
        val nPositions = readMap.keys.size

        //define emission probability
        val emissionProbabilityCalculator = EmissionProbabilityForViterbiHMM(readMap, likelyParentSet, probCorrect)
        val emissionP = emissionProbabilityCalculator::getLnHaploidEmissionProbabilityArray

        val initProbs = DoubleArray(nParents) {0.0}

        val result = viterbiOptimizedForHaploid(nParents, nPositions, initProbs, emissionP)

        //translate the result
        val positions = readMap.keys.sorted()
        val parentList = likelyParentSet.sorted()
        val resultList = result.first.mapIndexed { index,
                                                   i ->  Pair(
            Position(contig, positions[index]),
            gameteIndexMap[parentList[i]] ?: "none")}
        return resultList

    }

    /**
     *
     * Find the best diploid path for a contig. Inputs are the contig name, a map of gamete index to name, and
     * a map of int position (bin) to a list of Ps4gGameteSets. A Ps4gGameteSet contains an array of gamete indices and the count
     * of the number of read hits to that set of gametes. The output is a list of triples: position, genome1, genome2.
     */
    fun findDiploidPath(contig: String,
                        gameteIndexMap: Map<Int,String>,
                        readMap: Map<Int, MutableList<Ps4gGameteSet>>,
                        likelyParentSet: Set<Int>): List<Triple<Position, String, String>> {
        myLogger.info("Finding diploid path for contig $contig using ViterbiHMM")

        //create the transition matrix
        val nParents = likelyParentSet.size
        val nStates = nParents * nParents
        val nPositions = readMap.keys.size
        val transitionMatrix = DoubleArray(nStates * nStates)
        var ptr = 0
        for (index1 in 0 until nParents) {
            for (index2 in 0 until nParents) {
                for (index3 in 0 until nParents) {
                    for (index4 in 0 until nParents) {
                        transitionMatrix[ptr++] = PathFinderHMMPS4G
                            .DiploidTransitionWithInbreeding(sameGameteProbability, inbreedingCoefficient, nParents)
                            .calculate(Pair(index1, index2), Pair(index3, index4))
                    }
                }
            }
        }

        //emission probabilities (equal -1.0 for testing)
        val emissionProbabilityCalculator = EmissionProbabilityForViterbiHMM(readMap, likelyParentSet, probCorrect)
        val emissionP = emissionProbabilityCalculator::getLnDiploidEmissionProbabilityArray

        //val emissionP = { x: Int -> DoubleArray(nStates) {-1.0} }

        //initial probabilities, use inbreeding coefficent, but 0.0 for testing
        //probability of a homozygote = f, heterozygote = 1-f
        //probability of a specific homozygote = f/nParents heterozygoe = (1-f)/(nParents*nParents - nParents)
        val homozygoteProbabillity = inbreedingCoefficient / nParents
        val heterozygoteProbability = (1.0 - inbreedingCoefficient) / (nParents * nParents - nParents)
        val initProbs = DoubleArray(nParents * nParents) {heterozygoteProbability}
        for (ndx in 0 until nParents) {
            initProbs[ndx * nParents + ndx] = homozygoteProbabillity
        }

        val result = viterbiOptimized(nStates, nPositions, initProbs,
            transitionMatrix, emissionP)

        //translate the result
        val positions = readMap.keys.sorted()
        val parentList = likelyParentSet.sorted()
        val resultList = result.first.mapIndexed { index,
                                                   i ->  Triple(
            Position(contig, positions[index]),
            gameteIndexMap[parentList[i/nParents]] ?: "none",
            gameteIndexMap[parentList[i % nParents]] ?: "none")}
        return resultList
    }

    fun viterbiOptimized(
        stateCount: Int,
        positionCount: Int,
        initialLogProbabilities: DoubleArray,
        transitionLogProbabilities: DoubleArray,
        emissionLogProbabilityFunction: (positionIndex: Int) -> DoubleArray
    ): Pair<IntArray, Double> {
        require(stateCount > 0) { "State count must be positive" }
        require(positionCount > 0) { "State count must be positive" }
        require(initialLogProbabilities.size == stateCount)
        require(transitionLogProbabilities.size == stateCount * stateCount)

        // Rolling buffers: only the previous position's best log-probabilities are
        // needed to compute the current position's, so we avoid storing all positions.
        var previousBestLogProbability = DoubleArray(stateCount)
        var currentBestLogProbability = DoubleArray(stateCount)

        // backPointer must retain full history for backtracing at the end.
        // Flattened as [position * stateCount + stateIndex] for cache-friendly access.
        val backPointer = IntArray(positionCount * stateCount)

        // --- Initialization at the first position ---
        var emissionLogProbabilities = emissionLogProbabilityFunction(0)
        for (stateIndex in 0 until stateCount) {
            previousBestLogProbability[stateIndex] = initialLogProbabilities[stateIndex] + emissionLogProbabilities[stateIndex]
            backPointer[stateIndex] = -1 // no predecessor at position 0
        }

        // --- Recursion across remaining positions ---
        for (position in 1 until positionCount) {
            val backPointerRowOffset = position * stateCount
            emissionLogProbabilities = emissionLogProbabilityFunction(position)
            for (currentStateIndex in 0 until stateCount) {
                val emissionLogProbability =  emissionLogProbabilities[currentStateIndex]

                //here we make the same state the best state so that it will be chosen unless some other is better
                val transitionRowOffsetBase = currentStateIndex // used with previousStateIndex * stateCount below

                var bestPreviousLogProbability = previousBestLogProbability[currentStateIndex] +
                        transitionLogProbabilities[currentStateIndex * stateCount + transitionRowOffsetBase]
                var bestPreviousStateIndex = currentStateIndex

                for (previousStateIndex in 0 until stateCount) {
                    if (previousStateIndex == currentStateIndex) continue  //skip currentStateIndex
                    val transitionLogProbability =
                        transitionLogProbabilities[previousStateIndex * stateCount + transitionRowOffsetBase]

                    val candidateLogProbability =
                        previousBestLogProbability[previousStateIndex] + transitionLogProbability

                    if (candidateLogProbability > bestPreviousLogProbability) {
                        bestPreviousLogProbability = candidateLogProbability
                        bestPreviousStateIndex = previousStateIndex
                    }
                }

                currentBestLogProbability[currentStateIndex] = bestPreviousLogProbability + emissionLogProbability
                backPointer[backPointerRowOffset + currentStateIndex] = bestPreviousStateIndex
            }

            // Swap rolling buffers instead of allocating new arrays each position
            val temporaryBuffer = previousBestLogProbability
            previousBestLogProbability = currentBestLogProbability
            currentBestLogProbability = temporaryBuffer
        }

        // --- Find best final state at the last position ---
        var bestFinalStateIndex = 0
        var bestFinalLogProbability = previousBestLogProbability[0]

        for (stateIndex in 1 until stateCount) {
            if (previousBestLogProbability[stateIndex] > bestFinalLogProbability) {
                bestFinalLogProbability = previousBestLogProbability[stateIndex]
                bestFinalStateIndex = stateIndex
            }
        }

        // --- Backtrace to reconstruct the most likely state index sequence ---
        val bestStateIndicesByPosition = IntArray(positionCount)
        bestStateIndicesByPosition[positionCount - 1] = bestFinalStateIndex

        for (position in positionCount - 1 downTo 1) {
            val currentStateIndex = bestStateIndicesByPosition[position]
            bestStateIndicesByPosition[position - 1] =
                backPointer[position * stateCount + currentStateIndex]
        }
        return bestStateIndicesByPosition to bestFinalLogProbability

    }

    fun viterbiOptimizedForHaploid(
        stateCount: Int,
        positionCount: Int,
        initialLogProbabilities: DoubleArray,
        emissionLogProbabilityFunction: (positionIndex: Int) -> DoubleArray
    ): Pair<IntArray, Double> {
        require(stateCount > 0) { "State count must be positive" }
        require(positionCount > 0) { "State count must be positive" }
        require(initialLogProbabilities.size == stateCount)

        val lnNoSwitch = ln(sameGameteProbability)
        val lnSwitch = ln((1.0 - sameGameteProbability)/(stateCount - 1))
        // Rolling buffers: only the previous position's best log-probabilities are
        // needed to compute the current position's, so we avoid storing all positions.
        var previousBestLogProbability = DoubleArray(stateCount)
        var currentBestLogProbability = DoubleArray(stateCount)

        // backPointer must retain full history for backtracing at the end.
        // Flattened as [position * stateCount + stateIndex] for cache-friendly access.
        val backPointer = IntArray(positionCount * stateCount)

        // --- Initialization at the first position ---
        var emissionLogProbabilities = emissionLogProbabilityFunction(0)
        for (stateIndex in 0 until stateCount) {
            previousBestLogProbability[stateIndex] = initialLogProbabilities[stateIndex] + emissionLogProbabilities[stateIndex]
            backPointer[stateIndex] = -1 // no predecessor at position 0
        }

        // --- Recursion across remaining positions ---
        for (position in 1 until positionCount) {
            val backPointerRowOffset = position * stateCount
            emissionLogProbabilities = emissionLogProbabilityFunction(position)
            for (currentStateIndex in 0 until stateCount) {
                val emissionLogProbability =  emissionLogProbabilities[currentStateIndex]

                //here we make the same state the best state so that it will be chosen unless some other is better
                var samePathLogProbability = previousBestLogProbability[currentStateIndex] + lnNoSwitch

                //get the highest probability previous path
                var bestPreviousPathProbability = previousBestLogProbability[0]
                var bestPreviousStateIndex = 0
                for (previousStateIndex in 1 until stateCount) {
                    if (previousBestLogProbability[previousStateIndex] > bestPreviousPathProbability) {
                        bestPreviousPathProbability = previousBestLogProbability[previousStateIndex]
                        bestPreviousStateIndex = previousStateIndex
                    }
                }

                val switchPathLogProbability = bestPreviousPathProbability + lnSwitch
                if (switchPathLogProbability > samePathLogProbability) {
                    currentBestLogProbability[currentStateIndex] = switchPathLogProbability + emissionLogProbability
                    backPointer[backPointerRowOffset + currentStateIndex] = bestPreviousStateIndex
                } else {
                    currentBestLogProbability[currentStateIndex] = samePathLogProbability + emissionLogProbability
                    backPointer[backPointerRowOffset + currentStateIndex] = currentStateIndex
                }
            }

            // Swap rolling buffers instead of allocating new arrays each position
            val temporaryBuffer = previousBestLogProbability
            previousBestLogProbability = currentBestLogProbability
            currentBestLogProbability = temporaryBuffer
        }

        // --- Find best final state at the last position ---
        var bestFinalStateIndex = 0
        var bestFinalLogProbability = previousBestLogProbability[0]

        for (stateIndex in 1 until stateCount) {
            if (previousBestLogProbability[stateIndex] > bestFinalLogProbability) {
                bestFinalLogProbability = previousBestLogProbability[stateIndex]
                bestFinalStateIndex = stateIndex
            }
        }

        // --- Backtrace to reconstruct the most likely state index sequence ---
        val bestStateIndicesByPosition = IntArray(positionCount)
        bestStateIndicesByPosition[positionCount - 1] = bestFinalStateIndex

        for (position in positionCount - 1 downTo 1) {
            val currentStateIndex = bestStateIndicesByPosition[position]
            bestStateIndicesByPosition[position - 1] =
                backPointer[position * stateCount + currentStateIndex]
        }
        return bestStateIndicesByPosition to bestFinalLogProbability

    }

}