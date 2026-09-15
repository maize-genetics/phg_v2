package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import org.apache.logging.log4j.LogManager
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
 * product of their individual probabilities with no multinomial coefficient.
 *
 * ## Presence/absence correction
 *
 * The rule above assumes both founders have sequence at the site. Where one does not -- a deletion,
 * an assembly gap, or a region where its anchors are not unique -- every correctly mapped read
 * matches the other founder alone, and the model reads that as evidence for homozygosity. Each such
 * read contributes `ln 0.5` = -0.693 toward the homozygous state, so a run of roughly
 * `ln(pnn/psn) / 0.693` one-sided reads forces a false switch: about 32 at `--prob-same`
 * 0.999999999 with six candidate parents.
 *
 * Given a [PresenceTable], a founder whose anchor presence in the window is at or below
 * [pavThreshold] is treated as absent, and the heterozygous state then predicts exactly what the
 * homozygous state for the *present* founder predicts:
 *
 * | state (A, B), B absent | G contains | probability |
 * |---|---|---|
 * | | A | pc |
 * | | not A | pe |
 *
 * At full strength the two states tie, which is what the evidence supports: with B absent there is
 * genuinely nothing separating (A, A) from (A, B).
 *
 * [pavDamping] scales the correction: the divergent-site probability becomes `0.5^(1 - damping) *
 * pc` where the founder is flagged absent. At 0 the model is unchanged; at 1 the states tie, as
 * above. Intermediate values reduce the per-read penalty -- and so the rate at which a run of
 * one-sided reads accumulates toward a false switch -- while keeping the heterozygous state
 * strictly below the homozygous one.
 *
 * Scoring is monotonic in [pavDamping] and full strength was best on every benchmark tried: F1
 * arms (false-homozygous sequence down 65-90%), a balanced F2 panel where over-correction would
 * show up and does not overwhelm the gain (+0.35 concordance, +11.6 breakpoint F1, against a
 * measured cost of -0.14 AA recall and -1.47 breakpoint recall), and the maize
 * simulated-validation corpus. An earlier revision of this correction did behave badly at full
 * strength, but that was a bug -- it rewrote the wrong cell of the table, and its tell was a
 * non-monotonic damping sweep -- not a property of the model.
 *
 * State (B, B) is unaffected in all cases and still scores badly, correctly, since a B/B individual
 * could not produce reads where B is absent.
 *
 * This needs no founder-sharing matrix and no clamp; `probCorrect` and the presence threshold are
 * its only parameters.
 */
class GameteSetEmissionProbability(
    val readMap: Map<Int, MutableList<Ps4gGameteSet>>,
    parentSet: Set<Int>,
    val probCorrect: Double,
    val presenceTable: PresenceTable? = null,
    val contig: String = "",
    val gameteIndexMap: Map<Int, String> = emptyMap(),
    val binSize: Int = 256,
    val pavThreshold: Double = 0.02,
    val pavDamping: Double = 1.0
) {
    val parentList = parentSet.sorted()
    val nParents = parentList.size
    val positionList = readMap.keys.sorted()

    private val myLogger = LogManager.getLogger(GameteSetEmissionProbability::class.java)

    private val parentToLocal = HashMap<Int, Int>(nParents * 2).also { map ->
        parentList.forEachIndexed { local, global -> map[global] = local }
    }
    private val parentNames = parentList.map { gameteIndexMap[it] ?: "" }

    private val lnCorrect = ln(probCorrect)
    private val lnHalfCorrect = ln(0.5 * probCorrect)
    private val lnIncorrect = ln(1.0 - probCorrect)

    /** ln(0.5^(1-damping) * pc): the divergent-site value where a founder is flagged absent. */
    private val lnDampedCorrect =
        ln(Math.pow(0.5, 1.0 - pavDamping.coerceIn(0.0, 1.0)) * probCorrect)

    /** Reused across positions so no allocation happens per bin. */
    private val inGameteSet = BooleanArray(nParents)
    private val localBuffer = IntArray(nParents)

    private val useTable: Boolean
    private val absent = BooleanArray(nParents)
    private var cachedWindow = Int.MIN_VALUE
    private var flaggedWindows = 0
    private var totalWindows = 0

    init {
        val table = presenceTable
        useTable = table != null && contig.isNotEmpty() && table.hasContig(contig) &&
                parentNames.none { it.isEmpty() } &&
                table.missingTaxa(parentNames).isEmpty()
        if (table != null && !useTable) myLogger.warn(
            "Presence table unusable on $contig; running without the presence/absence correction"
        )
    }

    private fun setWindow(binPosition: Int) {
        val table = presenceTable!!
        val window = (binPosition.toLong() * binSize / table.windowSize).toInt()
        if (window == cachedWindow) return
        cachedWindow = window
        totalWindows++
        var any = false
        for (i in 0 until nParents) {
            val p = table.presence(contig, window, parentNames[i])
            absent[i] = !p.isNaN() && p <= pavThreshold
            if (absent[i]) any = true
        }
        if (any) flaggedWindows++
    }

    /**
     * Natural log emission probabilities at [positionIndex], indexed by ordered pairs of the sorted
     * parent list, as [ViterbiHMM.viterbiOptimized] expects.
     */
    fun getDiploidEmissionProbabilityArray(positionIndex: Int): DoubleArray {
        val position = positionList[positionIndex]
        if (useTable) setWindow(position)
        val gameteSets = readMap[position]!!
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
                    // Start from the uncorrected rule, then override the single case the
                    // presence/absence correction is about: exactly one founder in G, and the
                    // MISSING one is the one flagged absent. Every other case is untouched, so
                    // damping 0 reproduces the uncorrected model exactly.
                    val base = when {
                        i == j -> if (hasFirst) lnCorrect else lnIncorrect
                        hasFirst && hasSecond -> lnCorrect          // A = B at this site
                        hasFirst || hasSecond -> lnHalfCorrect      // A != B at this site
                        else -> lnIncorrect
                    }
                    val corrected = if (useTable && i != j && (hasFirst != hasSecond) &&
                        ((hasFirst && absent[j] && !absent[i]) || (hasSecond && absent[i] && !absent[j]))
                    ) lnDampedCorrect else base
                    probabilities[pointer++] += count * corrected
                }
            }
            for (a in 0 until size) inGameteSet[localBuffer[a]] = false
        }
        return probabilities
    }

    /** Fraction of visited windows in which at least one candidate founder was flagged absent. */
    fun flaggedWindowFraction(): Double =
        if (totalWindows == 0) 0.0 else flaggedWindows.toDouble() / totalWindows
}
