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
 * Emission probabilities come from [GameteSetEmissionProbability], which classifies each read's
 * site as identical or divergent between the two founders of the candidate pair. A single-path
 * (haploid) imputation is the same recursion at an inbreeding coefficient of 1, where only the
 * homozygous states are reachable and the emission reduces to that model's diagonal. Transition probabilities favor staying on the same
 * gamete(s) between adjacent bins, with a recombination penalty for switching. The inbreeding coefficient
 * affects the transition probabilities. Values > 1 favor transitions to a homozygous state.
 *
 * @param inbreedingCoefficient the inbreeding coefficient (0.0..1.0); used to set diploid initial
 *   state probabilities (probability of a homozygous vs. heterozygous state). Used by diploid path finding only.
 * @param sameGameteProbability the probability that the path stays on the same gamete when moving
 *   from one bin to the next (1 - recombination probability).
 * @param probCorrect the probability that a read maps to the correct haplotype; passed to the
 *   emission probability calculator.
 * @param binSize the bin size the ps4g file was created with, used to convert a bin position to a
 *   reference coordinate when looking a founder up in [presenceTable]. Diploid paths only.
 * @param presenceTable optional per-founder anchor presence, from `phg build-presence-table`. When
 *   supplied, a founder with no alignable sequence in a window stops being read as evidence for
 *   homozygosity. Diploid paths only; the haploid path is unaffected.
 * @param pavThreshold a founder whose anchor presence in a window is at or below this fraction is
 *   treated as absent there.
 * @param pavDamping strength of the presence/absence correction, 0 (off) to 1 (a heterozygous
 *   state involving an absent founder ties with the corresponding homozygous state).
 */
class ViterbiHMM(val inbreedingCoefficient: Double, val sameGameteProbability: Double, val probCorrect: Double,
                 val binSize: Int = 256,
                 val presenceTable: PresenceTable? = null,
                 val pavThreshold: Double = 0.02,
                 val pavDamping: Double = 1.0
) {
    private val myLogger = LogManager.getLogger(ViterbiHMM::class.java)

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

        val nParents = likelyParentSet.size
        val emissionCalculator = GameteSetEmissionProbability(readMap, likelyParentSet, probCorrect,
            presenceTable, contig, gameteIndexMap, binSize, pavThreshold, pavDamping)
        val result = findDiploidStatePath(
            nParents,
            readMap.keys.size,
            emissionCalculator::getDiploidEmissionProbabilityArray,
            emissionCalculator::getHomozygousEmissionProbabilityArray
        )

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
     * Runs the diploid recursion over the `nParents * nParents` ordered pairs of founders and returns
     * the chosen state path with its log probability, states indexed as `first * nParents + second`.
     *
     * Everything here depends only on the state count, the position count and the emission -- not on
     * where the observations came from -- so both a ps4g path and a VCF path use it. What differs
     * between them is how the emission is built and how state indices are turned back into names.
     *
     * The initial state distribution comes from [inbreedingCoefficient]: a homozygous state carries
     * `F / nParents` and a heterozygous one `(1 - F) / (nParents^2 - nParents)`.
     *
     * Both endpoints of the coefficient admit a cheaper recursion than the general scan, and neither
     * changes the path.
     *
     * At `F = 0` the transition is a Kronecker product of two haploid transitions, so the
     * maximisation separates and the recursion drops from `O(nParents^4)` to `O(nParents^2)` per
     * position, never materialising the transition matrix (3 MB at twenty-five parents).
     *
     * At `F = 1` every transition into a heterozygous state has probability zero -- see
     * [DiploidTransitionProbability.probabilityForF1] -- and so does every heterozygous initial
     * state, so only the `nParents` homozygous states are reachable. Among those the transition is
     * `pNoSwitch` to itself and `pSwitch` to any other, which is precisely the haploid transition
     * over `nParents` states, so the haploid recursion solves it exactly on the diagonal of the
     * emission array.
     *
     * Any coefficient strictly between the two couples the two founders and needs the general scan.
     *
     * @param emissionLogProbabilityFunction natural log emission for every ordered pair at a
     *   position, indexed `first * nParents + second`
     * @param homozygousEmissionLogProbabilityFunction the same for the homozygous states alone,
     *   indexed by founder. Used only at `F = 1`, where it saves computing `nParents^2` values to
     *   read `nParents` of them. Omit it and the diagonal is taken from the full array instead,
     *   which is correct but does that wasted work.
     * @param lnNoSwitchByPosition per-position log probability of staying on the same founder,
     *   indexed by the position being entered; see [transitionLogsByDistance]. Null uses the fixed
     *   value from [sameGameteProbability]. Supported at `F = 0` and `F = 1` only.
     * @param lnSwitchByPosition the matching per-position log probability of switching to one
     *   particular other founder.
     */
    fun findDiploidStatePath(
        nParents: Int,
        nPositions: Int,
        emissionLogProbabilityFunction: (positionIndex: Int) -> DoubleArray,
        homozygousEmissionLogProbabilityFunction: ((positionIndex: Int) -> DoubleArray)? = null,
        lnNoSwitchByPosition: DoubleArray? = null,
        lnSwitchByPosition: DoubleArray? = null
    ): Pair<IntArray, Double> {
        val nStates = nParents * nParents

        // A single candidate founder has no heterozygous state at all, and the divisor below would be
        // zero. The value is overwritten before use either way, but computing it would leave an
        // infinity or a NaN sitting in the array in the meantime.
        val homozygoteProbability = inbreedingCoefficient / nParents
        val heterozygoteProbability =
            if (nParents > 1) (1.0 - inbreedingCoefficient) / (nStates - nParents) else 0.0
        val initProbs = DoubleArray(nStates) { lnOrFloor(heterozygoteProbability) }
        for (ndx in 0 until nParents) {
            initProbs[ndx * nParents + ndx] = lnOrFloor(homozygoteProbability)
        }

        return when (inbreedingCoefficient) {
            0.0 -> viterbiOptimizedForDiploid(nParents, nPositions, initProbs,
                emissionLogProbabilityFunction, lnNoSwitchByPosition, lnSwitchByPosition)

            1.0 -> {
                val diagonalInit = DoubleArray(nParents) { initProbs[it * nParents + it] }
                val diagonalEmission = homozygousEmissionLogProbabilityFunction
                    ?: { positionIndex: Int ->
                        val full = emissionLogProbabilityFunction(positionIndex)
                        DoubleArray(nParents) { full[it * nParents + it] }
                    }
                val homozygous = viterbiOptimizedForHaploid(nParents, nPositions, diagonalInit,
                    diagonalEmission, lnNoSwitchByPosition, lnSwitchByPosition)
                // Re-express the chosen founders as diploid state indices (i, i).
                Pair(
                    IntArray(homozygous.first.size) {
                        homozygous.first[it] * nParents + homozygous.first[it]
                    },
                    homozygous.second
                )
            }

            else -> {
                // The general scan reads a precomputed nStates x nStates matrix, so a transition that
                // varies from step to step would mean rebuilding it at every position -- O(nParents^4)
                // of extra work per position, which is the same order as the scan itself. Refusing is
                // better than silently applying the first step's transition everywhere. The two
                // endpoints need only two scalars per step and take it in their stride.
                require(lnNoSwitchByPosition == null && lnSwitchByPosition == null) {
                    "Per-position transitions are supported at an inbreeding coefficient of 0 or 1, " +
                            "not at $inbreedingCoefficient: the general scan would need its " +
                            "transition matrix rebuilt at every position."
                }
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
                viterbiOptimized(nStates, nPositions, initProbs, matrix, emissionLogProbabilityFunction)
            }
        }
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
        emissionLogProbabilityFunction: (positionIndex: Int) -> DoubleArray,
        lnNoSwitchByPosition: DoubleArray? = null,
        lnSwitchByPosition: DoubleArray? = null
    ): Pair<IntArray, Double> {
        require(nParents > 0) { "Parent count must be positive" }
        require(positionCount > 0) { "Position count must be positive" }
        val stateCount = nParents * nParents
        require(initialLogProbabilities.size == stateCount)
        require(lnNoSwitchByPosition == null || lnNoSwitchByPosition.size == positionCount) {
            "lnNoSwitchByPosition must hold one value per position"
        }
        require(lnSwitchByPosition == null || lnSwitchByPosition.size == positionCount) {
            "lnSwitchByPosition must hold one value per position"
        }

        val fixedLnNoSwitch = ln(sameGameteProbability)
        val fixedLnSwitch =
            if (nParents > 1) ln((1.0 - sameGameteProbability) / (nParents - 1)) else fixedLnNoSwitch

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
            // Per-haplotype transition for the step into this position. Both stages below and the
            // whole-state comparison at the end use these, so a varying step size is charged
            // consistently across all three.
            val lnNoSwitch = lnNoSwitchByPosition?.get(position) ?: fixedLnNoSwitch
            val lnSwitch = lnSwitchByPosition?.get(position) ?: fixedLnSwitch

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
        emissionLogProbabilityFunction: (positionIndex: Int) -> DoubleArray,
        lnNoSwitchByPosition: DoubleArray? = null,
        lnSwitchByPosition: DoubleArray? = null
    ): Pair<IntArray, Double> {
        require(stateCount > 0) { "State count must be positive" }
        require(positionCount > 0) { "State count must be positive" }
        require(initialLogProbabilities.size == stateCount)
        require(lnNoSwitchByPosition == null || lnNoSwitchByPosition.size == positionCount) {
            "lnNoSwitchByPosition must hold one value per position"
        }
        require(lnSwitchByPosition == null || lnSwitchByPosition.size == positionCount) {
            "lnSwitchByPosition must hold one value per position"
        }

        val fixedLnNoSwitch = ln(sameGameteProbability)
        // With one state there is nowhere to switch to, and the divisor would be zero: the quotient
        // becomes positive infinity and a switch then looks infinitely attractive, wrecking the
        // score. The probability of switching is genuinely zero, so floor it.
        val fixedLnSwitch = if (stateCount > 1) ln((1.0 - sameGameteProbability) / (stateCount - 1))
                       else lnOrFloor(0.0)
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
            val lnNoSwitch = lnNoSwitchByPosition?.get(position) ?: fixedLnNoSwitch
            val lnSwitch = lnSwitchByPosition?.get(position) ?: fixedLnSwitch
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