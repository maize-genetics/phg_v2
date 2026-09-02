package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import kotlin.math.exp
import kotlin.math.ln

/**
 * Emission probabilities for the diploid HMM, modelled as a mixture over which of the two
 * haplotypes a read came from.
 *
 * Let `r[i][k]` be the probability that a read originating from founder `i` also matches founder
 * `k`, with `r[i][i] = probCorrect`. For an observed gamete set `G` -- the founders a read matched --
 *
 *     P(G | read came from i) = prod over k in G of r[i][k] * prod over k not in G of (1 - r[i][k])
 *
 * and a diploid state (i, j) emits that read as an even mixture of its two haplotypes:
 *
 *     ln P(G | state (i,j)) = ln( 0.5 * P(G | from i) + 0.5 * P(G | from j) )
 *
 * A homozygous state needs no special case: with i == j the mixture collapses to `P(G | from i)`
 * exactly, so one formula covers both.
 *
 * The point of this formulation is that **the partition is the observed gamete set, which does not
 * depend on the state**. Every state scores the same event, so their likelihoods are directly
 * comparable. An earlier four-class formulation (hits both / A only / B only / neither) failed
 * precisely here: it gave heterozygous states more cells than homozygous ones, so the same
 * probability mass was divided further and heterozygous states were penalised by `ln(1/q)` for
 * every read matching both founders -- around 41% of reads inside true heterozygous regions.
 *
 * Scoring against every founder, rather than only the two in the state, also means a read matching
 * some third founder informs the call instead of being lumped into an undifferentiated bucket.
 *
 * `r` is estimated from the reads supplied here, as the symmetric overlap
 * `both / (matches i + matches k - both)`.
 *
 * The conditional form `both / matches i` is the natural estimator but is unusable here: it is
 * contaminated by the sample's own composition. In a sample that is mostly founder X, nearly every
 * read matches X, so `P(matches X | matches anything)` approaches 1 for *every* founder and the
 * matrix stops describing sequence sharing at all. Measured on a 97%-B97 sample, every row gave
 * `r[.][B97] ~= 0.998`, and states pairing B97 with unrelated founders scored as well as B97/B97.
 *
 * The conditional estimate is sound for founders the sample actually carries and biased for those
 * it does not; the symmetric form borrows the estimate from whichever of the pair is better
 * represented, which is the one that is unbiased. It is a pragmatic correction, not a principled
 * one -- deriving `r` from the pangenome (haplotype identity per reference range, or founder-to-
 * founder anchor overlap) would be sample-independent by construction and is the right long-term
 * source.
 *
 * @param eps clamp holding every `r` inside [eps, 1 - eps] so no single read can drive a state to
 *   ln(0) and veto it for the rest of the contig.
 */
class MixtureEmissionProbability(
    val readMap: Map<Int, MutableList<Ps4gGameteSet>>,
    parentSet: Set<Int>,
    val probCorrect: Double,
    val eps: Double = 1e-3
) {
    val parentList = parentSet.sorted()
    val nParents = parentList.size
    val positionList = readMap.keys.sorted()

    private val parentToLocal = HashMap<Int, Int>(nParents * 2).also { map ->
        parentList.forEachIndexed { local, global -> map[global] = local }
    }

    /** ln(1 - r[i][k]) summed over all k: the likelihood of a read from i matching nothing. */
    private val baseline = DoubleArray(nParents)

    /** ln r[i][k] - ln(1 - r[i][k]): the correction applied for each founder actually matched. */
    private val matchDelta = DoubleArray(nParents * nParents)

    private val lnHalf = ln(0.5)
    private val localBuffer = IntArray(nParents)
    private val lnProbabilityFrom = DoubleArray(nParents)

    init {
        estimateSharing()
    }

    /** r[i][k] = symmetric read overlap between founders i and k, accumulated over every position. */
    private fun estimateSharing() {
        val bothTotal = LongArray(nParents * nParents)
        val matchTotal = LongArray(nParents)
        val local = IntArray(nParents)
        for ((_, gameteSets) in readMap) {
            for (gameteSet in gameteSets) {
                val count = gameteSet.count.toLong()
                var size = 0
                for (gamete in gameteSet.gameteIndices) {
                    val index = parentToLocal[gamete]
                    if (index != null) local[size++] = index
                }
                for (a in 0 until size) {
                    val i = local[a]
                    matchTotal[i] += count
                    for (b in 0 until size) bothTotal[i * nParents + local[b]] += count
                }
            }
        }
        for (i in 0 until nParents) {
            var sum = 0.0
            for (k in 0 until nParents) {
                val index = i * nParents + k
                val union = matchTotal[i] + matchTotal[k] - bothTotal[index]
                val raw = when {
                    i == k -> probCorrect
                    union <= 0L -> 0.5
                    else -> bothTotal[index].toDouble() / union
                }
                val r = raw.coerceIn(eps, 1.0 - eps)
                val lnNotR = ln(1.0 - r)
                matchDelta[index] = ln(r) - lnNotR
                sum += lnNotR
            }
            baseline[i] = sum
        }
    }

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
                if (index != null) localBuffer[size++] = index
            }
            // P(G | from i) for every candidate origin, as baseline plus one correction per match
            for (i in 0 until nParents) {
                var value = baseline[i]
                val row = i * nParents
                for (a in 0 until size) value += matchDelta[row + localBuffer[a]]
                lnProbabilityFrom[i] = value
            }
            var pointer = 0
            for (i in 0 until nParents) {
                val fromFirst = lnProbabilityFrom[i]
                for (j in 0 until nParents) {
                    probabilities[pointer++] += count * mixture(fromFirst, lnProbabilityFrom[j])
                }
            }
        }
        return probabilities
    }

    /** ln(0.5 * e^a + 0.5 * e^b), computed stably. Returns exactly `a` when a == b. */
    private fun mixture(a: Double, b: Double): Double {
        val high = if (a > b) a else b
        val low = if (a > b) b else a
        val difference = low - high
        return lnHalf + high + if (difference < -700.0) 0.0 else ln(1.0 + exp(difference))
    }
}
