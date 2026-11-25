package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.ContigSites.Companion.MAX_ALLELES_PER_SITE
import net.maizegenetics.phgv2.pathing.ropebwt.ContigSites.Companion.MISSING
import kotlin.math.ln

/**
 * Emission probabilities for a diploid HMM whose observations are **genotypes** rather than reads.
 *
 * A ps4g emission scores read mappings and multiplies its per-read term by a read count. A VCF
 * genotype is one call derived from an unknown depth, so there is no honest count to multiply by;
 * this model therefore scores one observation per site and nothing else. That is the whole reason it
 * exists rather than routing a VCF through a ps4g file.
 *
 * ## The model
 *
 * For a candidate founder pair the panel says what each founder carries, so the genotypes that pair
 * could produce are known by **Mendelian sampling**: each founder contributes one gamete drawn
 * uniformly from its own alleles, and the sample's genotype is the unordered pair of the two gametes.
 * With `pc` = [probCorrect] and `G` the observed genotype,
 *
 *     P(G | state) = pc * P(G | the pair's Mendelian sampling) + (1 - pc)
 *
 * and
 *
 *     P(G | (A,B)) = (1/|A|)(1/|B|) * #{ (alpha, beta) : alpha in A, beta in B, {alpha, beta} == G }
 *
 * where `|A|` counts A's alleles. Under the intended assumption -- a panel that is haploid or mostly
 * homozygous -- each founder has one allele, the pair's genotype is determined, and the inner
 * probability is 1 or 0:
 *
 * | state | predicts | `P(G \| state)` |
 * |---|---|---|
 * | (A, A) | {a, a} | 1 if `G = {a,a}`, else 0 |
 * | (A, B) | {a, b} | 1 if `G = {a,b}`, else 0 |
 *
 * Heterozygous panel sites need no special case, because the same formula covers them: with A
 * heterozygous for {a1, a2} and B homozygous, state (A, B) predicts {a1, b} and {a2, b} at one half
 * each. State (A, A) draws **two** gametes from A, which is selfing A -- {a1,a1}, {a1,a2}, {a2,a2} at
 * a quarter, a half, a quarter -- and is exactly right for an inbred descendant of a founder that was
 * not fully inbred itself.
 *
 * At `pc = 0.98` the resulting log emissions are
 *
 * | `P(G \| state)` | emission | ln |
 * |---|---|---|
 * | 1 | 1.0000 | 0.0000 |
 * | 1/2 | 0.5100 | -0.6733 |
 * | 1/4 | 0.2650 | -1.3280 |
 * | 0 | 0.0200 | -3.9120 |
 *
 * The `(1 - pc)` floor is not normalised over genotypes and does not need to be: it is the same for
 * every state at a position, so it cancels in every Viterbi comparison while keeping a mismatch
 * finite rather than negative infinity. [GameteSetEmissionProbability] uses the same device.
 *
 * ## Absent data
 *
 * A **missing sample genotype** leaves every state at 0.0, so the site neither favours nor penalises
 * any path -- it simply carries no information.
 *
 * A **missing founder allele** is scored as a founder about which nothing is known: the average over
 * the founders that *are* called at the site. A gap in the panel is evidence neither for nor against
 * a founder, and the two obvious constants each read it as one or the other. A mismatch would count it
 * against the founder. A perfect match -- ln 1, which is what this model once used -- counts it *for*
 * the founder: a founder missing at many sites then tied with the true founder at every one of them
 * and beat every called founder that mismatched, so the path was drawn toward founders with the most
 * missing data, and composition then copied their no-calls into the output.
 *
 * Averaging over the called founders is the same as giving the missing founder each allele at the
 * frequency `f` the called founders carry it, a heterozygous founder contributing half of each:
 *
 * | state | `P(G \| state)` before the `pc` mixing |
 * |---|---|
 * | (M, B), M missing | `(1/\|B\|) * sum over beta in B of f(the allele completing G with beta)` |
 * | (M, M), the same founder | the mean over called founders C of `P(G \| (C, C))` |
 * | (M, N), two missing founders | `f(g)^2` for `G = {g, g}`, `2 f(g1) f(g2)` for `G = {g1, g2}` |
 *
 * At a site where every called founder carries the same allele this is 1 or 0, exactly as for a called
 * founder. Where no founder at all is called, every state stays at 0.0, since there is nothing to
 * average and the site cannot tell any two states apart.
 *
 * A sample allele **no founder carries** needs no handling: every state mismatches, every state gets
 * `ln(1 - pc)`, and the site again discriminates between nothing. The same is true of a site that is
 * monomorphic across the panel, where every state predicts the same genotype. Both cost a little time
 * and could be filtered by the caller, but neither distorts the path.
 *
 * ## Cost
 *
 * `O(nFounders^2)` per site with at most four allele comparisons per state, matching the state space
 * the Viterbi has to scan anyway. Where every founder is homozygous -- the common case -- only one
 * comparison of the four is reached. Each founder's alleles are read once per site, and a site with a
 * missing founder costs `O(nFounders)` more for the averages. A faster form is available if a dense
 * panel ever makes this the bottleneck: precompute per site which founders carry each of the sample's
 * two alleles, reducing the inner work to two array reads, at the cost of a second code path to keep
 * in step with this one.
 *
 * @param sites the panel and sample genotypes for one contig, allele-index encoded
 * @param sampleIndex which sample of [ContigSites.sampleNames] this instance scores
 * @param probCorrect the probability that a genotype call is right
 */
class VcfGenotypeEmissionProbability(
    val sites: ContigSites,
    val sampleIndex: Int,
    val probCorrect: Double
) {
    private val nFounders = sites.nFounders

    /**
     * `ln(pc * matches/total + (1 - pc))` for every ratio the Mendelian count can produce.
     * With at most two alleles per founder, `total` is 1, 2 or 4, so indexing by `[total][matches]`
     * covers every case and keeps `ln` out of the inner loop.
     */
    private val lnEmission: Array<DoubleArray> = Array(5) { total ->
        DoubleArray(total + 1) { matches ->
            ln(probCorrect * matches.toDouble() / total.coerceAtLeast(1) + (1.0 - probCorrect))
        }
    }

    /** Each founder's distinct alleles at the site in hand, and how many: 0 when uncalled. */
    private val gameteCount = IntArray(nFounders)
    private val gamete1 = IntArray(nFounders)
    private val gamete2 = IntArray(nFounders)

    /** Allele frequency among the called founders at the site in hand, indexed by allele. */
    private val frequency = DoubleArray(MAX_ALLELES_PER_SITE)

    /** For a missing founder M: the log emission of (M, founder), and of (M, M). */
    private val lnWithMissing = DoubleArray(nFounders)

    /**
     * Natural log emission probabilities at [positionIndex], indexed by ordered pairs of founders as
     * `state = first * nFounders + second`, which is what [ViterbiHMM.viterbiOptimized] and its fast
     * paths expect.
     */
    fun getDiploidEmissionProbabilityArray(positionIndex: Int): DoubleArray {
        val probabilities = DoubleArray(nFounders * nFounders)

        val observed1 = sites.sampleAllele1(positionIndex, sampleIndex)
        val observed2 = sites.sampleAllele2(positionIndex, sampleIndex)
        // A missing genotype says nothing about any state, so every state stays at ln(1) = 0.
        if (observed1 == MISSING || observed2 == MISSING) return probabilities
        val o1 = observed1.toInt()
        val o2 = observed2.toInt()

        var missingFounders = 0
        for (founder in 0 until nFounders) {
            gameteCount[founder] = gametesOf(positionIndex, founder)
            if (gameteCount[founder] == 0) missingFounders++
        }
        // No founder is called, so there is nothing to average and every state is alike.
        if (missingFounders == nFounders) return probabilities
        var lnBothMissing = 0.0
        var lnSameMissing = 0.0
        if (missingFounders > 0) {
            val (both, same) = prepareMissing(o1, o2)
            lnBothMissing = ln(probCorrect * both + (1.0 - probCorrect))
            lnSameMissing = ln(probCorrect * same + (1.0 - probCorrect))
        }

        var pointer = 0
        for (first in 0 until nFounders) {
            val firstCount = gameteCount[first]
            for (second in 0 until nFounders) {
                val secondCount = gameteCount[second]
                if (firstCount == 0 || secondCount == 0) {
                    probabilities[pointer++] = when {
                        first == second -> lnSameMissing
                        firstCount == 0 && secondCount == 0 -> lnBothMissing
                        firstCount == 0 -> lnWithMissing[second]
                        else -> lnWithMissing[first]
                    }
                    continue
                }
                probabilities[pointer++] = lnEmission[firstCount * secondCount][matches(first, second, o1, o2)]
            }
        }
        return probabilities
    }

    /**
     * How many of the `|A| x |B|` gamete combinations of founders [first] and [second] give the observed
     * genotype {[o1], [o2]}. Unordered, since the sample VCF is unphased: {alpha, beta} and
     * {beta, alpha} are the same observation.
     */
    private fun matches(first: Int, second: Int, o1: Int, o2: Int): Int {
        var matches = 0
        for (a in 0 until gameteCount[first]) {
            val alpha = if (a == 0) gamete1[first] else gamete2[first]
            for (b in 0 until gameteCount[second]) {
                val beta = if (b == 0) gamete1[second] else gamete2[second]
                if ((alpha == o1 && beta == o2) || (alpha == o2 && beta == o1)) matches++
            }
        }
        return matches
    }

    /**
     * Fills [frequency] from the called founders and [lnWithMissing] from it, for the site in hand,
     * and returns the probabilities before the `pc` mixing for two different missing founders and for
     * one missing founder paired with itself. See the class notes for the three forms.
     */
    private fun prepareMissing(o1: Int, o2: Int): Pair<Double, Double> {
        var called = 0
        var sameSum = 0.0
        for (founder in 0 until nFounders) {
            val count = gameteCount[founder]
            if (count == 0) continue
            called++
            frequency[gamete1[founder]] += 1.0 / count
            if (count == 2) frequency[gamete2[founder]] += 0.5
            sameSum += matches(founder, founder, o1, o2).toDouble() / (count * count)
        }
        val f1 = frequency[o1] / called
        val f2 = frequency[o2] / called

        for (founder in 0 until nFounders) {
            val count = gameteCount[founder]
            if (count == 0) continue
            // the missing founder supplies whichever allele completes the genotype with beta
            var sum = 0.0
            for (b in 0 until count) {
                val beta = if (b == 0) gamete1[founder] else gamete2[founder]
                if (beta == o1) sum += f2
                if (beta == o2 && o1 != o2) sum += f1
            }
            lnWithMissing[founder] = ln(probCorrect * sum / count + (1.0 - probCorrect))
        }

        // reset only what was touched, so the next site starts from zero
        for (founder in 0 until nFounders) {
            val count = gameteCount[founder]
            if (count == 0) continue
            frequency[gamete1[founder]] = 0.0
            if (count == 2) frequency[gamete2[founder]] = 0.0
        }
        val both = if (o1 == o2) f1 * f1 else 2.0 * f1 * f2
        return both to sameSum / called
    }

    /**
     * Writes [founder]'s distinct alleles at [site] into [gamete1] and [gamete2] and returns how many
     * there are: 0 when the founder has no call, 1 when it is haploid or homozygous, 2 when
     * heterozygous.
     *
     * Collapsing a homozygous call to one allele is what makes the Mendelian count come out right.
     * Were both halves of `a/a` kept, state (A, B) would enumerate four combinations rather than two
     * and `matches/total` would still be correct, but state (A, A) with A homozygous would report
     * 4/4 rather than 1/1 -- the same probability, yet reached through a `total` of 4, which would
     * then differ from the heterozygous case in a way the lookup table could not distinguish.
     */
    private fun gametesOf(site: Int, founder: Int): Int {
        val allele1 = sites.founderAllele1(site, founder)
        if (allele1 == MISSING) return 0
        val allele2 = sites.founderAllele2(site, founder)
        gamete1[founder] = allele1.toInt()
        if (allele2 == MISSING || allele2 == allele1) return 1
        gamete2[founder] = allele2.toInt()
        return 2
    }
}
