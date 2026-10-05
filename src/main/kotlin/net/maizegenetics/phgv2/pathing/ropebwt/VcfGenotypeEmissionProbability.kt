package net.maizegenetics.phgv2.pathing.ropebwt

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
 * A **missing founder allele** leaves every state involving that founder at 0.0. The pair cannot be
 * evaluated, and the alternative -- scoring it as a mismatch -- would read a gap in the panel as
 * evidence against a founder, which it is not.
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
 * comparison of the four is reached. A faster form is available if a dense panel ever makes this the
 * bottleneck: precompute per site which founders carry each of the sample's two alleles, reducing the
 * inner work to two array reads, at the cost of a second code path to keep in step with this one.
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

    /** Reused across positions so the per-site call allocates only its result. */
    private val founderGametes = IntArray(2)
    private val otherGametes = IntArray(2)

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

        var pointer = 0
        for (first in 0 until nFounders) {
            val firstCount = gametesOf(positionIndex, first, founderGametes)
            for (second in 0 until nFounders) {
                if (firstCount == 0) {
                    // This founder has no call here, so the pair cannot be evaluated. Leaving it at
                    // 0.0 keeps a gap in the panel from counting as evidence against the founder.
                    pointer++
                    continue
                }
                val secondCount = gametesOf(positionIndex, second, otherGametes)
                if (secondCount == 0) {
                    pointer++
                    continue
                }
                var matches = 0
                for (a in 0 until firstCount) {
                    val alpha = founderGametes[a]
                    for (b in 0 until secondCount) {
                        val beta = otherGametes[b]
                        // Unordered comparison: the sample VCF is unphased, so {alpha, beta} and
                        // {beta, alpha} are the same observation.
                        if ((alpha == observed1.toInt() && beta == observed2.toInt()) ||
                            (alpha == observed2.toInt() && beta == observed1.toInt())
                        ) matches++
                    }
                }
                probabilities[pointer++] = lnEmission[firstCount * secondCount][matches]
            }
        }
        return probabilities
    }

    /**
     * Writes [founder]'s distinct alleles at [site] into [into] and returns how many there are: 0
     * when the founder has no call, 1 when it is haploid or homozygous, 2 when heterozygous.
     *
     * Collapsing a homozygous call to one allele is what makes the Mendelian count come out right.
     * Were both halves of `a/a` kept, state (A, B) would enumerate four combinations rather than two
     * and `matches/total` would still be correct, but state (A, A) with A homozygous would report
     * 4/4 rather than 1/1 -- the same probability, yet reached through a `total` of 4, which would
     * then differ from the heterozygous case in a way the lookup table could not distinguish.
     */
    private fun gametesOf(site: Int, founder: Int, into: IntArray): Int {
        val allele1 = sites.founderAllele1(site, founder)
        if (allele1 == MISSING) return 0
        val allele2 = sites.founderAllele2(site, founder)
        into[0] = allele1.toInt()
        if (allele2 == MISSING || allele2 == allele1) return 1
        into[1] = allele2.toInt()
        return 2
    }
}
