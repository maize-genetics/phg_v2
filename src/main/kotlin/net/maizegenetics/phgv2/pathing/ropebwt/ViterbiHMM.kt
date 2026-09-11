package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import net.maizegenetics.phgv2.utils.Position
import org.apache.logging.log4j.LogManager
import kotlin.math.ln

/**
 * Runs the Viterbi algorithm to infer the most likely path (or pair of paths) of parent gametes
 * through the ordered bins of a single contig, given the reads mapped to each bin.
 *
 * The hidden states are the candidate parent gametes: for a haploid path each state is a single
 * parent; for a diploid path each state is an ordered pair of parents (nParents * nParents states).
 * Emission probabilities come from [EmissionProbabilityForViterbiHMM] and are driven by which
 * gametes the reads in each bin hit. Transition probabilities favor staying on the same
 * gamete(s) between adjacent bins, with a recombination penalty for switching. The inbreeding coefficient
 * affects the transition probabilities. Values > 1 favor transitions to a homozygous state.
 *
 * @param inbreedingCoefficient the inbreeding coefficient (0.0..1.0); used to set diploid initial
 *   state probabilities (probability of a homozygous vs. heterozygous state). Used by diploid path finding only.
 * @param sameGameteProbability the probability that the path stays on the same gamete when moving
 *   from one bin to the next (1 - recombination probability).
 * @param probCorrect the probability that a read maps to the correct haplotype; passed to the
 *   emission probability calculator.
 * @param emissionModel "binomial" (default, historical behaviour) scores only whether a read hits
 *   either founder of the pair. "mixture" models the read as coming from one of the state's two
 *   haplotypes and scores its whole gamete set against a lift-derived sharing table.
 *   "gameteset" classifies each read's site as identical or divergent between the two founders
 *   directly from the gamete set, needing no sharing table. Both alternatives let a homozygous
 *   state be preferred on the evidence rather than through the inbreeding coefficient. Diploid
 *   paths only; the haploid path is unaffected.
 */
class ViterbiHMM(val inbreedingCoefficient: Double, val sameGameteProbability: Double, val probCorrect: Double,
                 val emissionModel: String = "binomial",
                 val sharingTable: SharingTable? = null,
                 val binSize: Int = 256,
                 val sharingClamp: Double = 0.40,
                 val sharingShrink: Double = 0.0,
                 val sharingMatchClamp: Double = 0.40,
                 val presenceTable: PresenceTable? = null,
                 val pavThreshold: Double = 0.02,
                 val pavDamping: Double = 1.0
) {
    private val myLogger = LogManager.getLogger(ViterbiHMM::class.java)

    /**
     * Finds the single best (haploid) path of parent gametes through the bins of a contig.
     *
     * @param contig the name of the contig being imputed.
     * @param gameteIndexMap map of gamete index to gamete (parent) name.
     * @param readMap map of bin position to the list of [Ps4gGameteSet]s observed in that bin; each
     *   [Ps4gGameteSet] holds the gamete indices hit and the count of reads hitting that set.
     * @param likelyParentSet the set of gamete indices to consider as candidate parents (states).
     * @return a list of (position, parent name) pairs, one per bin, ordered by bin position.
     */
    fun findHaploidPath(contig: String,
                        gameteIndexMap: Map<Int,String>,
                        readMap: Map<Int, MutableList<Ps4gGameteSet>>,
                        likelyParentSet: Set<Int>): List<Pair<Position, String>>  {

        val nParents = likelyParentSet.size
        val nPositions = readMap.keys.size

        //define emission probability
        val emissionProbabilityCalculator = EmissionProbabilityForViterbiHMM(readMap, likelyParentSet, probCorrect)
        val emissionP = emissionProbabilityCalculator::getHaploidEmissionProbabilityArray

        val initProbs = DoubleArray(nParents) {0.0}

        val result = viterbiOptimizedForHaploid(nParents, nPositions, initProbs, emissionP)

        //translate the result
        val positions = readMap.keys.sorted()
        val parentList = likelyParentSet.sorted()
        val resultList = result.first.mapIndexed { index, i ->
            //creates a list of Pair(Position, gamete name)
            Pair(Position(contig, positions[index]), gameteIndexMap[parentList[i]] ?: "none")
        }
        return resultList

    }

    /**
     * Finds the best diploid path (a pair of parent gametes per bin) through the bins of a contig.
     *
     * A full transition matrix over the nParents * nParents ordered parent pairs is built from
     * [sameGameteProbability] and [inbreedingCoefficient], and initial state probabilities are set
     * so that homozygous and heterozygous states are weighted by the inbreeding coefficient.
     *
     * @param contig the name of the contig being imputed.
     * @param gameteIndexMap map of gamete index to gamete (parent) name.
     * @param readMap map of bin position to the list of [Ps4gGameteSet]s observed in that bin; each
     *   [Ps4gGameteSet] holds the gamete indices hit and the count of reads hitting that set.
     * @param likelyParentSet the set of gamete indices to consider as candidate parents.
     * @return a list of (position, parent1 name, parent2 name) triples, one per bin, ordered by bin position.
     */
    fun findDiploidPath(contig: String,
                        gameteIndexMap: Map<Int,String>,
                        readMap: Map<Int, MutableList<Ps4gGameteSet>>,
                        likelyParentSet: Set<Int>): List<Triple<Position, String, String>> {

        //create the transition matrix
        val nParents = likelyParentSet.size
        val nStates = nParents * nParents
        val nPositions = readMap.keys.size
        // The general path needs the full matrix; the F = 0 fast path derives its two constants
        // directly and never materialises it (3 MB at twenty-five parents).
        val transitionMatrix = if (inbreedingCoefficient == 0.0) DoubleArray(0) else {
            val matrix = DoubleArray(nStates * nStates)
            val transitionProbabilityCalculator =
                DiploidTransitionProbability(sameGameteProbability, inbreedingCoefficient, nParents)
            var ptr = 0
            for (index1 in 0 until nParents) {
                for (index2 in 0 until nParents) {
                    for (index3 in 0 until nParents) {
                        for (index4 in 0 until nParents) {
                            matrix[ptr++] = transitionProbabilityCalculator
                                .calculateLn(Pair(index1, index2), Pair(index3, index4))
                        }
                    }
                }
            }
            matrix
        }

        //emission probabilities
        val emissionP = if (emissionModel == "gameteset") {
            GameteSetEmissionProbability(readMap, likelyParentSet, probCorrect, presenceTable,
                contig, gameteIndexMap, binSize, pavThreshold,
                pavDamping)::getDiploidEmissionProbabilityArray
        } else if (emissionModel == "mixture") {
            MixtureEmissionProbability(readMap, likelyParentSet, probCorrect, sharingTable, contig,
                gameteIndexMap, binSize, sharingClamp, sharingShrink,
                sharingMatchClamp)::getDiploidEmissionProbabilityArray
        } else {
            EmissionProbabilityForViterbiHMM(readMap, likelyParentSet, probCorrect)::getDiploidEmissionProbabilityArray
        }

        //val emissionP = { x: Int -> DoubleArray(nStates) {-1.0} }

        //initial probabilities, use inbreeding coefficent, but 0.0 for testing
        //probability of a homozygote = f, heterozygote = 1-f
        //probability of a specific homozygote = f/nParents heterozygoe = (1-f)/(nParents*nParents - nParents)
        val homozygoteProbabillity = inbreedingCoefficient / nParents
        val heterozygoteProbability = (1.0 - inbreedingCoefficient) / (nParents * nParents - nParents)
        val initProbs = DoubleArray(nParents * nParents) {lnOrFloor(heterozygoteProbability)}
        for (ndx in 0 until nParents) {
            initProbs[ndx * nParents + ndx] = lnOrFloor(homozygoteProbabillity)
        }

        // At F = 0 the transition is a Kronecker product of two haploid transitions, so the
        // maximisation separates and the recursion drops from O(nParents^4) to O(nParents^2) per
        // position. The paths are identical; only the cost differs. Any other inbreeding
        // coefficient couples the two founders and needs the general scan.
        val result = if (inbreedingCoefficient == 0.0) {
            viterbiOptimizedForDiploid(nParents, nPositions, initProbs, emissionP)
        } else {
            viterbiOptimized(nStates, nPositions, initProbs, transitionMatrix, emissionP)
        }

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

    /**
     * Viterbi for the diploid case when the inbreeding coefficient is zero, in O(nParents^2) per
     * position instead of O(nParents^4).
     *
     * At F = 0 the diploid transition is a Kronecker product of two independent haploid
     * transitions: moving (a, b) -> (c, d) costs `h[a->c] * h[b->d]`, where `h` is `pNoSwitch` when
     * the founder is unchanged and `pSwitch` otherwise. That reproduces `pnn`, `psn` and `pss`
     * exactly. The maximisation therefore separates:
     *
     *     max over (a,b) of [ V(a,b) + h(a->c) + h(b->d) ]
     *         = max over a of [ h(a->c) + max over b of ( V(a,b) + h(b->d) ) ]
     *
     * and because `h` takes only two values, the inner maximum is the row maximum of V unless the
     * row's argmax is `d` itself, in which case it is the row's second maximum. Each stage is one
     * pass over an n x n array, so the whole recursion is quadratic in the parent count rather than
     * quartic. It also never materialises the nStates x nStates transition matrix, which is 3 MB at
     * twenty-five parents.
     *
     * Ties resolve toward staying in the same state, matching [viterbiOptimized]; the two produce
     * identical paths.
     *
     * @return a pair of (best state index per position, log probability of that best path), state
     *   indices being `firstParent * nParents + secondParent` exactly as elsewhere.
     */
    fun viterbiOptimizedForDiploid(
        nParents: Int,
        positionCount: Int,
        initialLogProbabilities: DoubleArray,
        emissionLogProbabilityFunction: (positionIndex: Int) -> DoubleArray
    ): Pair<IntArray, Double> {
        require(nParents > 0) { "Parent count must be positive" }
        require(positionCount > 0) { "Position count must be positive" }
        val stateCount = nParents * nParents
        require(initialLogProbabilities.size == stateCount)

        val lnNoSwitch = ln(sameGameteProbability)
        val lnSwitch = if (nParents > 1) ln((1.0 - sameGameteProbability) / (nParents - 1)) else lnNoSwitch

        var previous = DoubleArray(stateCount)
        var current = DoubleArray(stateCount)
        val backPointer = IntArray(positionCount * stateCount)

        // intermediate: best over the second founder, for every (first founder, target second)
        val partial = DoubleArray(stateCount)
        val partialArg = IntArray(stateCount)

        var emission = emissionLogProbabilityFunction(0)
        for (state in 0 until stateCount) {
            previous[state] = initialLogProbabilities[state] + emission[state]
            backPointer[state] = -1
        }

        for (position in 1 until positionCount) {
            emission = emissionLogProbabilityFunction(position)
            val rowOffset = position * stateCount

            // --- stage one: maximise over the second founder of the previous state ---
            for (first in 0 until nParents) {
                val base = first * nParents
                var bestValue = Double.NEGATIVE_INFINITY; var bestIndex = 0
                var nextValue = Double.NEGATIVE_INFINITY; var nextIndex = -1
                for (second in 0 until nParents) {
                    val v = previous[base + second]
                    if (v > bestValue) {
                        nextValue = bestValue; nextIndex = bestIndex
                        bestValue = v; bestIndex = second
                    } else if (v > nextValue) {
                        nextValue = v; nextIndex = second
                    }
                }
                for (target in 0 until nParents) {
                    val keep = previous[base + target] + lnNoSwitch
                    val otherValue = if (bestIndex == target) nextValue else bestValue
                    val otherIndex = if (bestIndex == target) nextIndex else bestIndex
                    val switch = otherValue + lnSwitch
                    // ties resolve to the lower founder index, matching the full scan's ordering
                    // by state index; there is deliberately no bias toward keeping this founder,
                    // because the full scan only favours staying when the WHOLE state is unchanged
                    if (keep > switch || (keep == switch && target < otherIndex)) {
                        partial[base + target] = keep; partialArg[base + target] = target
                    } else {
                        partial[base + target] = switch; partialArg[base + target] = otherIndex
                    }
                }
            }

            // --- stage two: maximise over the first founder ---
            for (target2 in 0 until nParents) {
                var bestValue = Double.NEGATIVE_INFINITY; var bestIndex = 0
                var nextValue = Double.NEGATIVE_INFINITY; var nextIndex = -1
                for (first in 0 until nParents) {
                    val v = partial[first * nParents + target2]
                    if (v > bestValue) {
                        nextValue = bestValue; nextIndex = bestIndex
                        bestValue = v; bestIndex = first
                    } else if (v > nextValue) {
                        nextValue = v; nextIndex = first
                    }
                }
                for (target1 in 0 until nParents) {
                    val keep = partial[target1 * nParents + target2] + lnNoSwitch
                    val otherValue = if (bestIndex == target1) nextValue else bestValue
                    val otherIndex = if (bestIndex == target1) nextIndex else bestIndex
                    val switch = otherValue + lnSwitch
                    var chosenFirst: Int
                    var best: Double
                    if (keep > switch || (keep == switch && target1 < otherIndex)) {
                        best = keep; chosenFirst = target1
                    } else {
                        best = switch; chosenFirst = otherIndex
                    }
                    var chosenSecond = partialArg[chosenFirst * nParents + target2]
                    // the full scan seeds its search with the identical state, so an unchanged
                    // state wins every tie; reproduce that here rather than inside the stages
                    val state = target1 * nParents + target2
                    val stayWhole = previous[state] + lnNoSwitch + lnNoSwitch
                    if (stayWhole >= best) { best = stayWhole; chosenFirst = target1; chosenSecond = target2 }
                    current[state] = best + emission[state]
                    backPointer[rowOffset + state] = chosenFirst * nParents + chosenSecond
                }
            }

            val swap = previous; previous = current; current = swap
        }

        var bestFinalStateIndex = 0
        var bestFinalLogProbability = previous[0]
        for (state in 1 until stateCount) {
            if (previous[state] > bestFinalLogProbability) {
                bestFinalLogProbability = previous[state]; bestFinalStateIndex = state
            }
        }
        val best = IntArray(positionCount)
        best[positionCount - 1] = bestFinalStateIndex
        for (position in positionCount - 1 downTo 1) {
            best[position - 1] = backPointer[position * stateCount + best[position]]
        }
        return best to bestFinalLogProbability
    }

    /** ln(p), floored so that p == 0.0 yields a large finite penalty rather than -Inf. */
    private fun lnOrFloor(p: Double): Double = if (p <= 0.0) -1.0e6 else ln(p)

    /**
     * General Viterbi implementation that takes an explicit transition matrix. Used for the diploid
     * case where transition probabilities differ per state pair.
     *
     * Memory is kept low by rolling two log-probability buffers (only the previous position is
     * needed to compute the current one); only the back-pointer array retains full history for
     * backtracing. Arrays are flattened as [position * stateCount + stateIndex] for cache-friendly
     * access. Same-state transitions are preferred on ties so the path stays put unless another
     * state is strictly better.
     *
     * @param stateCount the number of hidden states.
     * @param positionCount the number of positions (bins) in the sequence.
     * @param initialLogProbabilities per-state initial log probabilities (length stateCount).
     * @param transitionLogProbabilities flattened stateCount x stateCount matrix of transition log
     *   probabilities, indexed as [previousState * stateCount + currentState].
     * @param emissionLogProbabilityFunction function returning the per-state emission log
     *   probabilities for a given position index.
     * @return a pair of (best state index per position, log probability of that best path).
     */
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

    /**
     * Viterbi implementation specialized for the haploid case. Instead of a full transition matrix,
     * transitions are described by just two log probabilities: staying on the same gamete
     * (ln(sameGameteProbability)) or switching to any of the other states
     * (ln((1 - sameGameteProbability) / (stateCount - 1))). This lets each position be computed by
     * finding the single best previous state once, rather than scanning all state pairs.
     *
     * Memory is kept low by rolling two log-probability buffers; only the back-pointer array retains
     * full history for backtracing. Arrays are flattened as [position * stateCount + stateIndex].
     *
     * @param stateCount the number of hidden states (candidate parents).
     * @param positionCount the number of positions (bins) in the sequence.
     * @param initialLogProbabilities per-state initial log probabilities (length stateCount).
     * @param emissionLogProbabilityFunction function returning the per-state emission log
     *   probabilities for a given position index.
     * @return a pair of (best state index per position, log probability of that best path).
     */
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
                val samePathLogProbability = previousBestLogProbability[currentStateIndex] + lnNoSwitch

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