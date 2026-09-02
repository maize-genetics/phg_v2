package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import org.apache.logging.log4j.LogManager
import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.pow

/**
 * Emission probabilities for the diploid HMM, modelled as a mixture over which of the two
 * haplotypes a read came from.
 *
 * Let `r[i][k]` be the probability that a read originating from founder `i` also matches founder
 * `k`, with `r[i][i] = probCorrect`. For an observed gamete set `G` -- the founders a read matched:
 *
 *     P(G | read came from i) = prod over k in G of r[i][k] * prod over k not in G of (1 - r[i][k])
 *     ln P(G | state (i,j)) = ln( 0.5 * P(G | from i) + 0.5 * P(G | from j) )
 *
 * The partition is the observed gamete set, which does not depend on the state, so every state
 * scores the same event and their likelihoods are comparable. A homozygous state needs no special
 * case: at i == j the mixture collapses exactly to `P(G | from i)`.
 *
 * An earlier four-class formulation (hits both / A only / B only / neither) failed here: it gave
 * heterozygous states more cells than homozygous ones, dividing the same probability mass further
 * and penalising heterozygous states by ln(1/q) for every read matching both founders -- about 41%
 * of reads inside true heterozygous regions.
 *
 * ## Where r comes from
 *
 * Sharing is intensely local: founders are near-identical over IBS blocks and divergent between
 * them. Measured on maize chr1, sharing between a founder pair runs from 0.10 to 0.85 across
 * megabase windows against a genome-wide mean near 0.40, so a single constant per pair discards
 * essentially all of the signal that distinguishes candidate founders.
 *
 * With a [SharingTable] supplied, `r` is therefore local: `r[i][k](w) = s[i][k](w) ^ gamma`, where
 * `s` is the index-derived anchor overlap in the window and `gamma` is a single exponent fitted
 * here from this sample's own reads.
 *
 * The exponent absorbs read length. Longer reads span more sequence and so match a second founder
 * less often, but not uniformly -- they suppress low-sharing windows far more than high-sharing
 * ones, because a long read almost certainly crosses a difference in a divergent window while
 * staying inside an IBS block in a conserved one. Measured across 100/150/250 bp, that distortion
 * is captured by a power law with a scale factor of 1.00 +/- 0.03 and an exponent that is the same
 * for every founder pair (1.30 +/- 0.02 for 100 -> 250 bp). Hence one fitted number per dataset,
 * and a sharing table that needs building only once per index.
 *
 * Without a table, `r` falls back to a single genome-wide symmetric overlap per pair, computed from
 * the reads. That form works but leaves spurious third-founder states; see the notes.
 *
 * @param matchEps clamp for the MATCHED factor, kept separate from [eps].
 *
 * The two factors of this likelihood want opposite things. `ln(1 - r)` in the complement factor
 * diverges as r approaches 1, so it needs a tight clamp. But het-versus-hom discrimination lives
 * entirely in the matched factor, and clamping there compresses it: measured inside true
 * heterozygous intervals, the margin `ln P(AB) - ln P(AA)` falls from 2.05 per read at a clamp of
 * 0.001 to 0.67 at 0.40, while the mixture's fixed `ln 2` cost stays at 0.693. A single clamp
 * therefore appears to trade heterozygote recall for complement stability.
 *
 * That tension turned out to be mostly an artefact of clamping `r[i][i]`, which is `probCorrect`
 * -- a supplied parameter, not an estimate. Clamping rewrote it as 0.60 at a clamp of 0.4, which
 * cost 8 points of heterozygote recall on its own. With the diagonal exempted, equal clamps are
 * best and splitting them is actively worse (97.4% vs 93.2% at e=0), so this parameter defaults
 * to the same value as [eps] and exists only for further study.
 *
 * @param shrink pulls the local sharing toward that pair's contig-wide mean before the exponent is
 *   applied: `s' = (1 - shrink) * s(w) + shrink * mean(s)`. 0 uses the local value as measured, 1
 *   reduces to the global-sharing behaviour.
 *
 * @param eps clamp holding every `r` inside [eps, 1 - eps].
 *
 * The clamp is not just a guard against ln(0). The complement factor `prod (1 - r[i][k])` runs over
 * every founder the read did NOT match -- about twenty of them -- and `ln(1 - r)` diverges as r
 * approaches 1. With a single global sharing value the term is identical for every founder and
 * cancels out of the comparison entirely; with local sharing, which reaches 0.99, it becomes a
 * founder-specific offset spreading 5.7 log units on average and up to 25, swinging by 2.1 between
 * adjacent windows against a recombination penalty of only 10.8. The path then follows the
 * complement term rather than the reads. A clamp bounds ln(1 - r) and keeps that term in proportion.
 * Swept on simulated F2 panels, 0.2-0.4 is a broad optimum: concordance rises from 80.7% at 0.001
 * to 93.9% at 0.2 and plateaus, while imputed breakpoints fall from 1771 to 328. The control that
 * matters is that the same clamp applied to *global* sharing gets steadily worse (77.3% to 67.7%),
 * because clamping destroys the only signal a global value carries. Local sharing does the work;
 * the clamp only stops the complement term from drowning it.
 */
class MixtureEmissionProbability(
    val readMap: Map<Int, MutableList<Ps4gGameteSet>>,
    parentSet: Set<Int>,
    val probCorrect: Double,
    val sharingTable: SharingTable? = null,
    val contig: String = "",
    val gameteIndexMap: Map<Int, String> = emptyMap(),
    val binSize: Int = 256,
    val eps: Double = 0.40,
    val shrink: Double = 0.0,
    val matchEps: Double = 0.40
) {
    val parentList = parentSet.sorted()
    val nParents = parentList.size
    val positionList = readMap.keys.sorted()

    private val myLogger = LogManager.getLogger(MixtureEmissionProbability::class.java)

    private val parentToLocal = HashMap<Int, Int>(nParents * 2).also { map ->
        parentList.forEachIndexed { local, global -> map[global] = local }
    }
    private val parentNames = parentList.map { gameteIndexMap[it] ?: "" }

    /** ln(1 - r[i][k]) summed over k: the likelihood of a read from i matching nothing. */
    private val baseline = DoubleArray(nParents)

    /** ln r[i][k] - ln(1 - r[i][k]): the correction applied for each founder actually matched. */
    private val matchDelta = DoubleArray(nParents * nParents)

    private val lnHalf = ln(0.5)
    private val localBuffer = IntArray(nParents)
    private val lnProbabilityFrom = DoubleArray(nParents)

    private val useTable: Boolean
    private var gamma = 1.0
    private var cachedWindow = Int.MIN_VALUE
    private val contigMeanSharing = DoubleArray(nParents * nParents)

    init {
        val table = sharingTable
        useTable = table != null && contig.isNotEmpty() && table.hasContig(contig) &&
                table.missingTaxa(parentNames.filter { it.isNotEmpty() }).isEmpty() &&
                parentNames.none { it.isEmpty() }
        if (useTable) {
            computeContigMeanSharing(table!!)
            gamma = fitExponent(table)
            myLogger.info("Mixture emission on $contig: local sharing, fitted exponent gamma = %.3f".format(gamma))
        } else {
            if (table != null) myLogger.warn(
                "Mixture emission on $contig: sharing table unusable here, falling back to a " +
                        "genome-wide sharing estimate from the reads"
            )
            setGlobalSharing()
        }
    }

    /** That pair's mean sharing over every window of this contig, the shrinkage target. */
    private fun computeContigMeanSharing(table: SharingTable) {
        val windows = table.windowCount(contig)
        if (windows == 0) return
        for (i in 0 until nParents) {
            for (k in 0 until nParents) {
                var sum = 0.0
                var seen = 0
                for (w in 0 until windows) {
                    val s = table.sharing(contig, w, parentNames[i], parentNames[k])
                    if (!s.isNaN()) { sum += s; seen++ }
                }
                contigMeanSharing[i * nParents + k] = if (seen > 0) sum / seen else 0.5
            }
        }
    }

    private fun windowOf(binPosition: Int) = (binPosition.toLong() * binSize / sharingTable!!.windowSize).toInt()

    /** Read-derived sharing accumulated per window, used only to fit [gamma]. */
    private fun fitExponent(table: SharingTable): Double {
        val windowBoth = HashMap<Int, LongArray>()
        val windowMatch = HashMap<Int, LongArray>()
        val local = IntArray(nParents)
        for ((position, gameteSets) in readMap) {
            val window = windowOf(position)
            val both = windowBoth.getOrPut(window) { LongArray(nParents * nParents) }
            val match = windowMatch.getOrPut(window) { LongArray(nParents) }
            for (gameteSet in gameteSets) {
                val count = gameteSet.count.toLong()
                var size = 0
                for (gamete in gameteSet.gameteIndices) {
                    val index = parentToLocal[gamete]
                    if (index != null) local[size++] = index
                }
                for (a in 0 until size) {
                    match[local[a]] += count
                    for (b in 0 until size) both[local[a] * nParents + local[b]] += count
                }
            }
        }
        // regression through the origin in log-log space: the scale factor is 1 by measurement
        var numerator = 0.0
        var denominator = 0.0
        for ((window, both) in windowBoth) {
            val match = windowMatch[window]!!
            for (i in 0 until nParents) {
                for (k in i + 1 until nParents) {
                    val union = match[i] + match[k] - both[i * nParents + k]
                    if (union < MIN_READS_FOR_FIT) continue
                    val observed = both[i * nParents + k].toDouble() / union
                    val expected = table.sharing(contig, window, parentNames[i], parentNames[k])
                    if (observed <= FIT_FLOOR || observed >= 1.0 - FIT_FLOOR) continue
                    if (expected.isNaN() || expected <= FIT_FLOOR || expected >= 1.0 - FIT_FLOOR) continue
                    val logExpected = ln(expected)
                    numerator += logExpected * ln(observed)
                    denominator += logExpected * logExpected
                }
            }
        }
        return if (denominator <= 0.0) 1.0 else (numerator / denominator).coerceIn(0.25, 4.0)
    }

    /** Fallback: one symmetric overlap per pair for the whole contig. */
    private fun setGlobalSharing() {
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
                    matchTotal[local[a]] += count
                    for (b in 0 until size) bothTotal[local[a] * nParents + local[b]] += count
                }
            }
        }
        val values = DoubleArray(nParents * nParents)
        for (i in 0 until nParents) {
            for (k in 0 until nParents) {
                val index = i * nParents + k
                val union = matchTotal[i] + matchTotal[k] - bothTotal[index]
                values[index] = when {
                    i == k -> probCorrect
                    union <= 0L -> 0.5
                    else -> bothTotal[index].toDouble() / union
                }
            }
        }
        applySharing(values)
    }

    /** Rebuild [baseline] and [matchDelta] for the window containing [binPosition]. */
    private fun setWindowSharing(binPosition: Int) {
        val table = sharingTable!!
        val window = windowOf(binPosition)
        if (window == cachedWindow) return
        cachedWindow = window
        val values = DoubleArray(nParents * nParents)
        for (i in 0 until nParents) {
            for (k in 0 until nParents) {
                val index = i * nParents + k
                values[index] = if (i == k) probCorrect else {
                    val raw = table.sharing(contig, window, parentNames[i], parentNames[k])
                    val s = if (raw.isNaN()) contigMeanSharing[index] else raw
                    val shrunk = (1.0 - shrink) * s + shrink * contigMeanSharing[index]
                    shrunk.coerceIn(1e-6, 1.0 - 1e-6).pow(gamma)
                }
            }
        }
        applySharing(values)
    }

    /**
     * [values] holds unclamped sharing. The complement factor is built from a tightly clamped copy
     * and the matched factor from a loosely clamped one, so stability and discrimination can be set
     * independently.
     */
    private fun applySharing(values: DoubleArray) {
        for (i in 0 until nParents) {
            var sum = 0.0
            for (k in 0 until nParents) {
                val index = i * nParents + k
                // The diagonal is probCorrect, a supplied parameter rather than an estimate, so it
                // is never clamped. Clamping it silently rewrites "a read from i matches i" as
                // 0.70 at a clamp of 0.3, or 0.60 at 0.4, which distorts every state equally in
                // the matched factor and unequally once the two factors are clamped apart.
                val raw = values[index]
                val complement = if (i == k) raw else raw.coerceIn(eps, 1.0 - eps)
                val matched = if (i == k) raw else raw.coerceIn(matchEps, 1.0 - matchEps)
                val lnNotR = ln((1.0 - complement).coerceAtLeast(1e-12))
                matchDelta[index] = ln(matched) - lnNotR
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
        val position = positionList[positionIndex]
        if (useTable) setWindowSharing(position)
        val gameteSets = readMap[position]!!
        val probabilities = DoubleArray(nParents * nParents)

        for (gameteSet in gameteSets) {
            val count = gameteSet.count
            var size = 0
            for (gamete in gameteSet.gameteIndices) {
                val index = parentToLocal[gamete]
                if (index != null) localBuffer[size++] = index
            }
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

    companion object {
        private const val MIN_READS_FOR_FIT = 50L
        private const val FIT_FLOOR = 0.02
    }
}
