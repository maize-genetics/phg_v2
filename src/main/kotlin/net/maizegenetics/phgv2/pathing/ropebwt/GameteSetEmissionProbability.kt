package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import kotlin.math.ln

/**
 * Emission probabilities for the diploid HMM that classify each read's site as identical-by-state
 * or divergent between the two founders of the candidate state, reading that classification off the
 * observed gamete set rather than from any precomputed sharing estimate.
 *
 * For a candidate state (A, B) and an observed gamete set G, split the reads that mapped correctly
 * by whether A and B carry the same sequence at that site:
 *
 *     P(G | AB) = P(G | AB, A=B) P(A=B) + P(G | AB, A!=B) P(A!=B)
 *
 * A correctly mapped read from an A=B site matches both founders, so it appears as `ab`. A
 * correctly mapped read from an A!=B site matches exactly one, with equal chance of either. With
 * `pc` = [probCorrect] and `pe` = 1 - pc, that gives the whole model:
 *
 * | state | G contains | probability |
 * |---|---|---|
 * | (A, A) | A | pc |
 * | (A, A) | not A | pe |
 * | (A, B) | both A and B | pc |
 * | (A, B) | exactly one of A, B | 0.5 * pc |
 * | (A, B) | neither | pe |
 *
 * A homozygous state is just the case where every site is an A=B site, which is why the first two
 * rows carry no factor of one half.
 *
 * Gamete sets combine as the probability of the single observed *sequence* of gamete sets, i.e. the
 * product of their individual probabilities with no multinomial coefficient. The sequence was
 * realised once, so its probability is the relevant quantity; this also matches the probability of
 * the path given no switches.
 *
 * ## What this does not model
 *
 * Presence/absence variation. Where founder B lacks the sequence entirely, every correctly mapped
 * read matches A alone, and the model has no way to distinguish that from a genuine A!=B site --
 * so it reads a PAV block as evidence for the homozygous state. Each such read contributes
 * `ln 0.5` = -0.693 toward homozygous, so a run of roughly `ln(pnn/psn) / 0.693` one-sided reads is
 * enough to force a switch: about 16 at the default `--prob-same` of 0.9999 with six candidate
 * parents. Whether that matters in practice is an empirical question about how long real one-sided
 * blocks are, and whether the recombination prior alone is enough to absorb them.
 *
 * Unlike [MixtureEmissionProbability] this needs no lift-derived sharing table and no clamp, and
 * `probCorrect` is its only parameter.
 */
class GameteSetEmissionProbability(
    val readMap: Map<Int, MutableList<Ps4gGameteSet>>,
    parentSet: Set<Int>,
    val probCorrect: Double
) {
    val parentList = parentSet.sorted()
    val nParents = parentList.size
    val positionList = readMap.keys.sorted()

    private val parentToLocal = HashMap<Int, Int>(nParents * 2).also { map ->
        parentList.forEachIndexed { local, global -> map[global] = local }
    }

    private val lnCorrect = ln(probCorrect)
    private val lnHalfCorrect = ln(0.5 * probCorrect)
    private val lnIncorrect = ln(1.0 - probCorrect)

    /** Reused across positions so no allocation happens per bin. */
    private val inGameteSet = BooleanArray(nParents)
    private val localBuffer = IntArray(nParents)

    /**
     * Natural log emission probabilities at [positionIndex], indexed by ordered pairs of the sorted
     * parent list, as [ViterbiHMM.viterbiOptimized] expects.
     */
    fun getDiploidEmissionProbabilityArray(positionIndex: Int): DoubleArray {
        val gameteSets = readMap[positionList[positionIndex]]!!
        val probabilities = DoubleArray(nParents * nParents)

        for (gameteSet in gameteSets) {
            val count = gameteSet.count
            var size = 0
            for (gamete in gameteSet.gameteIndices) {
                val index = parentToLocal[gamete]
                if (index != null) {
                    localBuffer[size++] = index
                    inGameteSet[index] = true
                }
            }
            var pointer = 0
            for (i in 0 until nParents) {
                val hasFirst = inGameteSet[i]
                for (j in 0 until nParents) {
                    val hasSecond = inGameteSet[j]
                    probabilities[pointer++] += count * when {
                        i == j -> if (hasFirst) lnCorrect else lnIncorrect
                        hasFirst && hasSecond -> lnCorrect          // A = B at this site
                        hasFirst || hasSecond -> lnHalfCorrect      // A != B at this site
                        else -> lnIncorrect
                    }
                }
            }
            for (a in 0 until size) inGameteSet[localBuffer[a]] = false
        }
        return probabilities
    }
}
