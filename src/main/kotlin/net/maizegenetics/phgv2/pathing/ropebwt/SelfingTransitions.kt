package net.maizegenetics.phgv2.pathing.ropebwt

import kotlin.math.expm1
import kotlin.math.ln
import kotlin.math.ln1p
import kotlin.math.pow

/**
 * Diploid transitions for a line derived from a biparental F1 by [generations] generations of
 * selfing, S_N = F_(N+1): NAM RILs at S5 are `generations = 5`.
 *
 * The F = 0 diploid model lets the two haplotypes switch independently, so a homozygous breakpoint,
 * (A, A) -> (B, B), costs two switches -- the square of a tiny probability -- while (A, A) -> (A, B)
 * costs one. In a selfed line it is the other way round. Both haplotypes descend from the same
 * crossover once the segment has been fixed, so a homozygous breakpoint is the commonest kind, and a
 * heterozygous segment is short because each further generation of selfing fixes or loses it.
 *
 * The transitions here are the exact two-locus ones. Two loci a recombination fraction `r` apart
 * are followed from the F1 through [generations] rounds of selfing, each a meiosis with one
 * recombination chance between them, and the step from the first locus's ordered genotype to the
 * second's is read off as a conditional probability. To first order in `r`, per unit `r`:
 *
 *     (A, A) -> (B, B)           (2^(N+1) - 2N - 2) / (2^N - 1)    -> 2, the RIL value, as N grows
 *     (A, A) -> (A, B) or (B, A)  N / (2^N - 1) each
 *     (A, B) -> (A, A) or (B, B)  N each
 *     (A, B) -> (B, A)            0 (order r^2)
 *
 * and a locus is heterozygous with probability 2^-N. At N = 5 that is 0.81 homozygous and 0.31
 * heterozygous breakpoints over half a Morgan, where the F = 0 model expects none of the first.
 *
 * ## The recombination fraction is still a smoothing parameter
 *
 * `r` comes from `--prob-switch` and `--prob-switch-distance`, read as the probability of at least
 * one crossover across that distance in one meiosis, and converted to a recombination fraction with
 * Haldane's map. So [generations] fixes the *shape* of the transitions -- how homozygous breakpoints
 * weigh against heterozygous ones and how quickly a heterozygous segment ends -- and the tuned rate
 * still sets how freely the path moves overall. The tuned rate is many orders of magnitude below real
 * recombination for the reason given in [transitionLogsByDistance].
 *
 * ## More than two founders
 *
 * The model is biparental, and `--n-parents 2` is the intended use. With K founders a homozygous
 * state (a, a) shares each kind of change among the K - 1 other founders, which is how
 * [transitionLogsByDistance] shares a switch: (a, a) -> (c, c) takes the biparental (A, A) -> (B, B)
 * probability divided by K - 1, and (a, a) -> (a, c) takes (A, A) -> (A, B) divided by K - 1. It
 * cannot gain two new founders, (a, a) -> (c, d), since a biparental cross offers only one. A
 * heterozygous state (a, b) has already named its pair, so it moves only among (a, a), (b, b) and
 * (b, a), with the biparental probabilities unshared. Each row still sums to one, and the initial
 * probabilities stay the stationary distribution; sharing a heterozygote's changes among all the
 * founders would break that, draining homozygotes into heterozygous states of a cross that does not
 * exist.
 *
 * The transition is a dense `K^2 x K^2` matrix, scanned in full at every position, so this is for
 * small founder sets; see [MAX_PARENTS].
 *
 * @param generations generations of selfing after the F1, at least 1
 * @param nParents number of candidate founders, K
 */
class SelfingTransitions(val generations: Int, val nParents: Int) {

    companion object {
        /** Largest founder set accepted: the recursion costs `nParents^4` per position. */
        const val MAX_PARENTS = 8

        /** ln(0), floored like the other recursions so an impossible step is finite. */
        private const val LN_ZERO = -1.0e6

        /** Biparental ordered genotypes, indexed `first * 2 + second` with A = 0 and B = 1. */
        const val AA = 0
        const val AB = 1
        const val BA = 2
        const val BB = 3

        /**
         * Recombination fraction across [distance] when [probSwitch] is the probability of at least
         * one crossover across [referenceDistance]: crossovers are a Poisson process at
         * `-ln(1 - probSwitch) / referenceDistance` per unit, and Haldane's map gives
         * `r = (1 - exp(-2 * rate * distance)) / 2`. Computed with expm1, since the exponent is
         * normally far below the double-precision step at 1.
         */
        fun recombinationFraction(distance: Double, probSwitch: Double, referenceDistance: Double): Double {
            val morgans = -ln1p(-probSwitch) / referenceDistance * distance
            return -expm1(-2.0 * morgans) / 2.0
        }
    }

    init {
        require(generations >= 1) { "selfing generations must be at least 1, not $generations" }
        require(nParents in 1..MAX_PARENTS) {
            "the selfing model scans a dense nParents^4 transition at every position, so it takes " +
                    "1 to $MAX_PARENTS parents, not $nParents; restrict them with --n-parents"
        }
    }

    private val nStates = nParents * nParents

    /** Probability that a locus is heterozygous: 2^-N. */
    val heterozygosity: Double = 0.5.pow(generations)

    /**
     * For each (from, to) pair of K-founder states, flattened `from * nStates + to`: the biparental
     * transition it maps to, or -1 if it has none.
     */
    private val biparentalIndex = IntArray(nStates * nStates)

    /** ln of how many K-founder targets share that biparental transition. */
    private val lnShare = DoubleArray(nStates * nStates)

    init {
        val lnAlternatives = if (nParents > 1) ln((nParents - 1).toDouble()) else 0.0
        for (from in 0 until nStates) {
            val a = from / nParents
            val b = from % nParents
            for (to in 0 until nStates) {
                val c = to / nParents
                val d = to % nParents
                val cell = from * nStates + to
                val firstChanged = c != a
                val secondChanged = d != b
                if (a == b) {
                    if (firstChanged && secondChanged && c != d) {
                        biparentalIndex[cell] = -1
                    } else {
                        // Relabel a as A; any other founder is B.
                        val target = (if (firstChanged) 2 else 0) + (if (secondChanged) 1 else 0)
                        biparentalIndex[cell] = AA * 4 + target
                        lnShare[cell] = if (firstChanged || secondChanged) lnAlternatives else 0.0
                    }
                } else if ((c == a || c == b) && (d == a || d == b)) {
                    // Relabel a as A and b as B.
                    val target = (if (c == a) 0 else 2) + (if (d == a) 0 else 1)
                    biparentalIndex[cell] = AB * 4 + target
                    lnShare[cell] = 0.0
                } else {
                    biparentalIndex[cell] = -1
                }
            }
        }
    }

    /**
     * Log initial probabilities over the K-founder states, indexed `first * nParents + second`: the
     * single-locus genotype frequencies after [generations] generations, 1 - 2^-N homozygous and
     * 2^-N heterozygous, shared equally among the states of each kind. This is also the stationary
     * distribution of the transitions.
     */
    fun initialLogProbabilities(): DoubleArray = DoubleArray(nStates) { state ->
        if (state / nParents == state % nParents) ln((1.0 - heterozygosity) / nParents)
        else ln(heterozygosity / (nParents * (nParents - 1)))
    }

    /**
     * The biparental two-locus transition across recombination fraction [r], as 16 probabilities
     * flattened `from * 4 + to` over [AA], [AB], [BA] and [BB]: the probability of each genotype at
     * a locus given the genotype at a locus a distance `r` before it.
     *
     * Two loci are followed through F1 -> [generations] selfings in one plant line. The two parents
     * are founder A (0) and founder B (1).
     * - A two-locus haplotype is the founders it carries at the two loci, indexed
     *   `(founder at locus 1) * 2 + (founder at locus 2)`: 0 = A-A, 1 = A-B, 2 = B-A, 3 = B-B.
     * - An individual is an ordered pair of haplotypes (haplotype 1, haplotype 2), indexed
     *   `h1 * 4 + h2`, so there are 16 kinds. `individuals[i]` is the probability that the plant in
     *   the current generation is of kind `i`. The pair is ordered because the HMM states are: AB
     *   (haplotype 1 from A, haplotype 2 from B) is a different state from BA.
     *
     * The F1 is the A-A haplotype paired with the B-B haplotype. Each selfing draws two independent
     * gametes from the plant, one for each of the offspring's two haplotypes. A gamete is either of the plant's
     * own haplotypes with probability (1 - r) / 2 each, or either recombinant with probability r / 2
     * each. A recombinant takes locus 1 from one haplotype and locus 2 from the other.
     *
     * After the last generation, each individual is read as a genotype at each locus:
     * - the genotype at locus 1 is (haplotype 1's founder at locus 1, haplotype 2's founder at locus 1);
     * - the genotype at locus 2 is the same at locus 2.
     * Summing gives the joint distribution of the two genotypes. Dividing each row by its total gives
     * the transition.
     *
     * The leaving probabilities are sums of products with no cancellation, so they stay accurate when
     * `r` is far below the double-precision step at 1, as a tuned `--prob-switch` makes it. Staying is
     * formed as one minus their sum.
     */
    fun biparentalTransitions(r: Double): DoubleArray {
        require(r in 0.0..0.5) { "recombination fraction must be between 0 and 0.5, not $r" }
        fun atLocus1(haplotype: Int) = haplotype / 2
        fun atLocus2(haplotype: Int) = haplotype % 2
        fun haplotype(locus1: Int, locus2: Int) = locus1 * 2 + locus2

        var individuals = DoubleArray(16)
        individuals[haplotype(0, 0) * 4 + haplotype(1, 1)] = 1.0  // the F1
        val gametes = DoubleArray(4)  // gamete haplotype -> probability, for one individual
        repeat(generations) {
            val next = DoubleArray(16)
            for (individual in 0 until 16) {
                val p = individuals[individual]
                if (p == 0.0) continue
                val h1 = individual / 4
                val h2 = individual % 4
                gametes.fill(0.0)
                gametes[h1] += (1.0 - r) / 2.0
                gametes[h2] += (1.0 - r) / 2.0
                gametes[haplotype(atLocus1(h1), atLocus2(h2))] += r / 2.0
                gametes[haplotype(atLocus1(h2), atLocus2(h1))] += r / 2.0
                // the offspring takes one gamete as haplotype 1 and an independent one as haplotype 2
                for (g1 in 0 until 4) for (g2 in 0 until 4) next[g1 * 4 + g2] += p * gametes[g1] * gametes[g2]
            }
            individuals = next
        }

        // Joint distribution, flattened `(genotype at locus 1) * 4 + (genotype at locus 2)`, each
        // genotype being (haplotype 1's founder) * 2 + (haplotype 2's founder), which is AA, AB, BA or BB.
        val joint = DoubleArray(16)
        for (individual in 0 until 16) {
            val h1 = individual / 4
            val h2 = individual % 4
            val genotype1 = atLocus1(h1) * 2 + atLocus1(h2)
            val genotype2 = atLocus2(h1) * 2 + atLocus2(h2)
            joint[genotype1 * 4 + genotype2] += individuals[individual]
        }

        // transition from -> to = P(locus 2 is `to` | locus 1 is `from`)
        val transitions = DoubleArray(16)
        for (from in 0 until 4) {
            val marginal = (0 until 4).sumOf { joint[from * 4 + it] }
            var leaving = 0.0
            for (to in 0 until 4) {
                if (to == from) continue
                transitions[from * 4 + to] = joint[from * 4 + to] / marginal
                leaving += transitions[from * 4 + to]
            }
            transitions[from * 4 + from] = 1.0 - leaving
        }
        return transitions
    }

    /**
     * ln of [biparentalTransitions], with staying taken as ln1p of minus the leaving total so it
     * keeps its precision when leaving is below the double-precision step at 1.
     */
    fun biparentalLogTransitions(r: Double): DoubleArray {
        val p = biparentalTransitions(r)
        val logs = DoubleArray(16)
        for (from in 0 until 4) {
            var leaving = 0.0
            for (to in 0 until 4) {
                if (to == from) continue
                leaving += p[from * 4 + to]
                logs[from * 4 + to] = if (p[from * 4 + to] > 0.0) ln(p[from * 4 + to]) else LN_ZERO
            }
            logs[from * 4 + from] = ln1p(-leaving)
        }
        return logs
    }

    /**
     * Fills [out] with the K-founder log transition matrix, flattened `from * nStates + to`, built
     * from the biparental [biparentalLogs] as described on the class.
     */
    fun fillLogMatrix(biparentalLogs: DoubleArray, out: DoubleArray) {
        // With one founder there is one state and nowhere to go.
        if (nParents == 1) {
            out[0] = 0.0
            return
        }
        for (cell in out.indices) {
            val index = biparentalIndex[cell]
            out[cell] = if (index < 0) LN_ZERO else biparentalLogs[index] - lnShare[cell]
        }
    }

    /** The K-founder log transition matrix across recombination fraction [r]. */
    fun logMatrix(r: Double): DoubleArray =
        DoubleArray(nStates * nStates).also { fillLogMatrix(biparentalLogTransitions(r), it) }

    /**
     * Per-position K-founder log transition matrices for steps of [distances], as from
     * [stepDistances], for the general Viterbi scan. Element `i` governs the step into position `i`.
     *
     * The biparental transitions are cached by distance -- binned positions repeat the same few step
     * lengths -- and the returned function expands one into a single reused buffer, so the matrix it
     * returns is valid only until the next call.
     */
    fun logMatricesByDistance(
        distances: IntArray,
        probSwitch: Double,
        referenceDistance: Double
    ): (positionIndex: Int) -> DoubleArray {
        require(probSwitch > 0.0 && probSwitch < 1.0) { "prob-switch must be between 0 and 1 exclusive" }
        require(referenceDistance > 0.0) { "reference distance must be positive" }
        val cache = HashMap<Int, DoubleArray>()
        // As in transitionLogsByDistance, a zero-length step is floored at one unit.
        val byPosition = Array(distances.size) { position ->
            val distance = distances[position].coerceAtLeast(1)
            cache.getOrPut(distance) {
                biparentalLogTransitions(recombinationFraction(distance.toDouble(), probSwitch, referenceDistance))
            }
        }
        val buffer = DoubleArray(nStates * nStates)
        return { position -> fillLogMatrix(byPosition[position], buffer); buffer }
    }
}
