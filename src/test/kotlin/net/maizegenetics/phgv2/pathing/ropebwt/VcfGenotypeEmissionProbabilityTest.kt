package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.ContigSites.Companion.MISSING
import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test
import kotlin.math.ln
import kotlin.random.Random

class VcfGenotypeEmissionProbabilityTest {

    private val probCorrect = 0.98
    private val match = 0.0
    private val mismatch = ln(1.0 - probCorrect)
    private val half = ln(0.5 * probCorrect + (1.0 - probCorrect))
    private val quarter = ln(0.25 * probCorrect + (1.0 - probCorrect))

    /**
     * Builds one site. Each founder and each sample is given as a pair of allele indices, with
     * [MISSING] for a no-call; a homozygous or haploid call repeats its allele.
     */
    private fun oneSite(
        founders: List<Pair<Int, Int>>,
        sample: Pair<Int, Int>,
        position: Int = 100
    ): VcfGenotypeEmissionProbability {
        val founderAlleles = ByteArray(founders.size * 2)
        founders.forEachIndexed { index, (first, second) ->
            founderAlleles[index * 2] = first.toByte()
            founderAlleles[index * 2 + 1] = second.toByte()
        }
        val sites = ContigSites(
            founderNames = founders.indices.map { "founder$it" },
            sampleNames = listOf("sample"),
            positions = intArrayOf(position),
            founderAlleles = founderAlleles,
            sampleAlleles = byteArrayOf(sample.first.toByte(), sample.second.toByte())
        )
        return VcfGenotypeEmissionProbability(sites, 0, probCorrect)
    }

    private fun DoubleArray.state(first: Int, second: Int, nFounders: Int) = this[first * nFounders + second]

    // ---------------------------------------------------------------- the homozygous-panel table

    @Test
    fun aHomozygousPairPredictsExactlyOneGenotype() {
        // founder0 = allele 0, founder1 = allele 1. Sample is heterozygous 0/1, so only the pairs
        // that combine one of each explain it.
        val p = oneSite(listOf(0 to 0, 1 to 1), 0 to 1).getDiploidEmissionProbabilityArray(0)
        assertEquals(match, p.state(0, 1, 2), 1e-12, "(f0, f1) predicts {0,1}")
        assertEquals(match, p.state(1, 0, 2), 1e-12, "the reversed pair predicts the same, unphased")
        assertEquals(mismatch, p.state(0, 0, 2), 1e-12, "(f0, f0) predicts {0,0}")
        assertEquals(mismatch, p.state(1, 1, 2), 1e-12, "(f1, f1) predicts {1,1}")
    }

    @Test
    fun aHomozygousSampleIsExplainedOnlyByTheMatchingHomozygousPair() {
        val p = oneSite(listOf(0 to 0, 1 to 1), 1 to 1).getDiploidEmissionProbabilityArray(0)
        assertEquals(match, p.state(1, 1, 2), 1e-12)
        assertEquals(mismatch, p.state(0, 0, 2), 1e-12)
        assertEquals(mismatch, p.state(0, 1, 2), 1e-12)
        assertEquals(mismatch, p.state(1, 0, 2), 1e-12)
    }

    @Test
    fun twoFoundersCarryingTheSameAlleleAreIndistinguishable() {
        // A site monomorphic across the panel discriminates between nothing, which is correct rather
        // than a defect -- it just costs time.
        val p = oneSite(listOf(0 to 0, 0 to 0), 0 to 0).getDiploidEmissionProbabilityArray(0)
        assertTrue(p.all { it == match }, "every state predicts {0,0}: ${p.toList()}")
    }

    // ---------------------------------------------------------------- Mendelian sampling

    @Test
    fun aHeterozygousFounderSplitsItsGameteBetweenTwoAlleles() {
        // founder0 is het {0,1}, founder1 is homozygous 0. State (f0, f1) predicts {0,0} and {1,0},
        // each at one half.
        val founders = listOf(0 to 1, 0 to 0)
        val onZeroZero = oneSite(founders, 0 to 0).getDiploidEmissionProbabilityArray(0)
        assertEquals(half, onZeroZero.state(0, 1, 2), 1e-12, "{0,0} is one of the two outcomes")
        val onZeroOne = oneSite(founders, 0 to 1).getDiploidEmissionProbabilityArray(0)
        assertEquals(half, onZeroOne.state(0, 1, 2), 1e-12, "{0,1} is the other")
        val onOneOne = oneSite(founders, 1 to 1).getDiploidEmissionProbabilityArray(0)
        assertEquals(mismatch, onOneOne.state(0, 1, 2), 1e-12, "{1,1} needs two copies of allele 1")
    }

    @Test
    fun aHomozygousStateOnAHeterozygousFounderIsSelfing() {
        // State (f0, f0) with f0 het {0,1} draws two independent gametes from f0, so it predicts
        // {0,0}, {0,1}, {1,1} at a quarter, a half, a quarter. This is the row most easily got
        // wrong: collapsing it to "predicts {0,1}" would be treating the founder as a single gamete.
        val founders = listOf(0 to 1)
        assertEquals(quarter, oneSite(founders, 0 to 0).getDiploidEmissionProbabilityArray(0)[0], 1e-12)
        assertEquals(half, oneSite(founders, 0 to 1).getDiploidEmissionProbabilityArray(0)[0], 1e-12)
        assertEquals(quarter, oneSite(founders, 1 to 1).getDiploidEmissionProbabilityArray(0)[0], 1e-12)
    }

    @Test
    fun aHomozygousStateOnAHomozygousFounderIsCertain() {
        // The counterpart of the row above: with f0 homozygous there is only one gamete, so the
        // selfed genotype is determined. The gamete collapse is what keeps this at 1/1 rather than
        // arriving at 4/4 through a different denominator.
        val p = oneSite(listOf(0 to 0), 0 to 0).getDiploidEmissionProbabilityArray(0)
        assertEquals(match, p[0], 1e-12)
        assertEquals(mismatch, oneSite(listOf(0 to 0), 0 to 1).getDiploidEmissionProbabilityArray(0)[0], 1e-12)
    }

    @Test
    fun bothFoundersHeterozygousGivesQuarters() {
        // f0 het {0,1}, f1 het {0,1}. State (f0, f1) has four equally likely gamete combinations:
        // {0,0} once, {0,1} twice, {1,1} once.
        val founders = listOf(0 to 1, 0 to 1)
        assertEquals(quarter, oneSite(founders, 0 to 0).getDiploidEmissionProbabilityArray(0).state(0, 1, 2), 1e-12)
        assertEquals(half, oneSite(founders, 0 to 1).getDiploidEmissionProbabilityArray(0).state(0, 1, 2), 1e-12)
        assertEquals(quarter, oneSite(founders, 1 to 1).getDiploidEmissionProbabilityArray(0).state(0, 1, 2), 1e-12)
    }

    @Test
    fun theMendelianCountMatchesABruteForceEnumeration() {
        // An independent implementation over random sites, including heterozygous founders and
        // multi-allelic sites, so the optimised gamete collapse cannot quietly disagree with the
        // definition it is meant to implement.
        val random = Random(20260924)
        val nFounders = 4
        val nAlleles = 3
        repeat(300) {
            val founders = (0 until nFounders).map {
                if (random.nextInt(4) == 0) MISSING.toInt() to MISSING.toInt()
                else random.nextInt(nAlleles) to random.nextInt(nAlleles)
            }
            val sample = random.nextInt(nAlleles) to random.nextInt(nAlleles)
            val actual = oneSite(founders, sample).getDiploidEmissionProbabilityArray(0)

            for (i in 0 until nFounders) for (j in 0 until nFounders) {
                val gi = gametesByHand(founders[i])
                val gj = gametesByHand(founders[j])
                val expected = if (gi.isEmpty() || gj.isEmpty()) 0.0 else {
                    var matches = 0
                    for (a in gi) for (b in gj) {
                        if ((a == sample.first && b == sample.second) ||
                            (a == sample.second && b == sample.first)) matches++
                    }
                    ln(probCorrect * matches.toDouble() / (gi.size * gj.size) + (1.0 - probCorrect))
                }
                assertEquals(expected, actual.state(i, j, nFounders), 1e-12,
                    "founders=$founders sample=$sample state=($i,$j)")
            }
        }
    }

    /** Distinct gametes of a founder, written out longhand for the brute-force check. */
    private fun gametesByHand(founder: Pair<Int, Int>): List<Int> = when {
        founder.first == MISSING.toInt() -> emptyList()
        founder.second == MISSING.toInt() || founder.second == founder.first -> listOf(founder.first)
        else -> listOf(founder.first, founder.second)
    }

    // ---------------------------------------------------------------- absent data

    @Test
    fun aMissingSampleGenotypeLeavesEveryStateNeutral() {
        val p = oneSite(listOf(0 to 0, 1 to 1), MISSING.toInt() to MISSING.toInt())
            .getDiploidEmissionProbabilityArray(0)
        assertTrue(p.all { it == 0.0 }, "the site must carry no information: ${p.toList()}")
    }

    @Test
    fun aMissingFounderLeavesOnlyItsOwnPairsNeutral() {
        // founder1 has no call. Pairs involving it cannot be evaluated, so they stay at 0.0 -- but
        // the pairs that do not involve it must still be scored normally, or a single gap in the
        // panel would wipe out the whole site.
        val p = oneSite(listOf(0 to 0, MISSING.toInt() to MISSING.toInt(), 1 to 1), 0 to 1)
            .getDiploidEmissionProbabilityArray(0)
        assertEquals(0.0, p.state(0, 1, 3), "pair with the missing founder is neutral")
        assertEquals(0.0, p.state(1, 0, 3))
        assertEquals(0.0, p.state(1, 1, 3))
        assertEquals(0.0, p.state(1, 2, 3))
        assertEquals(match, p.state(0, 2, 3), 1e-12, "founders 0 and 2 still explain {0,1}")
        assertEquals(mismatch, p.state(0, 0, 3), 1e-12, "and are still penalised where they do not")
    }

    @Test
    fun missingIsNotTreatedAsAMismatch() {
        // The distinction that matters: a gap in the panel must score higher than a founder that is
        // present and wrong, otherwise absence becomes evidence against.
        val p = oneSite(listOf(MISSING.toInt() to MISSING.toInt(), 0 to 0), 1 to 1)
            .getDiploidEmissionProbabilityArray(0)
        assertTrue(p.state(0, 0, 2) > p.state(1, 1, 2),
            "missing (${p.state(0, 0, 2)}) must beat present-and-wrong (${p.state(1, 1, 2)})")
    }

    @Test
    fun anAlleleNoFounderCarriesDiscriminatesBetweenNothing() {
        // Sample carries allele 2, which no founder has. Every state mismatches, so the site is
        // uninformative rather than misleading -- no handling needed, but worth pinning.
        val p = oneSite(listOf(0 to 0, 1 to 1), 2 to 2).getDiploidEmissionProbabilityArray(0)
        assertTrue(p.all { it == mismatch }, "every state equally wrong: ${p.toList()}")
    }

    // ---------------------------------------------------------------- shape and multiple sites

    @Test
    fun theArrayIsIndexedByOrderedPairsAsTheViterbiExpects() {
        val p = oneSite(listOf(0 to 0, 1 to 1, 0 to 0), 0 to 1).getDiploidEmissionProbabilityArray(0)
        assertEquals(9, p.size, "nFounders squared")
        // founders 0 and 2 both carry allele 0, so (0,1) (1,0) (2,1) (1,2) all explain {0,1}
        assertEquals(match, p.state(2, 1, 3), 1e-12)
        assertEquals(match, p.state(1, 2, 3), 1e-12)
        assertEquals(mismatch, p.state(0, 2, 3), 1e-12, "two copies of allele 0 cannot give {0,1}")
    }

    @Test
    fun eachSiteIsScoredFromItsOwnGenotypes() {
        // Two sites where the answer flips, to catch an index applied to the wrong row.
        val sites = ContigSites(
            founderNames = listOf("f0", "f1"),
            sampleNames = listOf("s"),
            positions = intArrayOf(10, 20),
            // site 0: f0 = 0/0, f1 = 1/1    site 1: the founders swap, f0 = 1/1, f1 = 0/0
            founderAlleles = byteArrayOf(0, 0, 1, 1, 1, 1, 0, 0),
            // the sample is 0/0 at both sites, so only the founders move and the answer must flip
            sampleAlleles = byteArrayOf(0, 0, 0, 0)
        )
        val emission = VcfGenotypeEmissionProbability(sites, 0, probCorrect)
        val site0 = emission.getDiploidEmissionProbabilityArray(0)
        assertEquals(match, site0.state(0, 0, 2), 1e-12, "site 0: f0 carries the sample's allele")
        assertEquals(mismatch, site0.state(1, 1, 2), 1e-12, "site 0: f1 does not")
        val site1 = emission.getDiploidEmissionProbabilityArray(1)
        assertEquals(mismatch, site1.state(0, 0, 2), 1e-12, "site 1: f0 no longer carries it")
        assertEquals(match, site1.state(1, 1, 2), 1e-12, "site 1: f1 now does")
    }

    @Test
    fun eachSampleIsScoredFromItsOwnGenotypes() {
        val sites = ContigSites(
            founderNames = listOf("f0", "f1"),
            sampleNames = listOf("s0", "s1"),
            positions = intArrayOf(10),
            founderAlleles = byteArrayOf(0, 0, 1, 1),
            // s0 is 0/0, s1 is 1/1
            sampleAlleles = byteArrayOf(0, 0, 1, 1)
        )
        val first = VcfGenotypeEmissionProbability(sites, 0, probCorrect).getDiploidEmissionProbabilityArray(0)
        val second = VcfGenotypeEmissionProbability(sites, 1, probCorrect).getDiploidEmissionProbabilityArray(0)
        assertEquals(match, first.state(0, 0, 2), 1e-12)
        assertEquals(mismatch, first.state(1, 1, 2), 1e-12)
        assertEquals(mismatch, second.state(0, 0, 2), 1e-12)
        assertEquals(match, second.state(1, 1, 2), 1e-12)
    }

    @Test
    fun repeatedCallsForOnePositionAgree() {
        // The class reuses scratch arrays across calls, so a leak between calls would show here.
        val emission = oneSite(listOf(0 to 1, 1 to 1, MISSING.toInt() to MISSING.toInt()), 0 to 1)
        val first = emission.getDiploidEmissionProbabilityArray(0)
        val second = emission.getDiploidEmissionProbabilityArray(0)
        assertArrayEquals(first, second, 1e-12)
    }

    @Test
    fun probCorrectSetsTheMismatchPenaltyAndNothingElse() {
        for (pc in listOf(0.9, 0.98, 0.999)) {
            val sites = ContigSites(listOf("f0", "f1"), listOf("s"), intArrayOf(10),
                byteArrayOf(0, 0, 1, 1), byteArrayOf(0, 1))
            val p = VcfGenotypeEmissionProbability(sites, 0, pc).getDiploidEmissionProbabilityArray(0)
            assertEquals(0.0, p.state(0, 1, 2), 1e-12, "a match is ln(1) whatever pc is")
            assertEquals(ln(1.0 - pc), p.state(0, 0, 2), 1e-12, "a mismatch is ln(1 - pc)")
        }
    }

    // ---------------------------------------------------------------- the container's own checks

    @Test
    fun contigSitesRejectsArraysThatDoNotMatchItsDimensions() {
        val tooShort = assertThrows(IllegalArgumentException::class.java) {
            ContigSites(listOf("f0", "f1"), listOf("s"), intArrayOf(10, 20),
                byteArrayOf(0, 0, 1, 1), byteArrayOf(0, 0, 1, 1))
        }
        assertTrue(tooShort.message!!.contains("founderAlleles"), tooShort.message)

        val sampleTooShort = assertThrows(IllegalArgumentException::class.java) {
            ContigSites(listOf("f0"), listOf("s0", "s1"), intArrayOf(10),
                byteArrayOf(0, 0), byteArrayOf(0, 0))
        }
        assertTrue(sampleTooShort.message!!.contains("sampleAlleles"), sampleTooShort.message)
    }

    @Test
    fun contigSitesAccessorsReadTheRightCell() {
        val sites = ContigSites(
            founderNames = listOf("f0", "f1"),
            sampleNames = listOf("s0", "s1"),
            positions = intArrayOf(10, 20),
            // site 0: f0 = 1/2, f1 = 3/4 ; site 1: f0 = 5/6, f1 = 7/8
            founderAlleles = byteArrayOf(1, 2, 3, 4, 5, 6, 7, 8),
            // site 0: s0 = 9/10, s1 = 11/12 ; site 1: s0 = 13/14, s1 = 15/16
            sampleAlleles = byteArrayOf(9, 10, 11, 12, 13, 14, 15, 16)
        )
        assertEquals(1.toByte(), sites.founderAllele1(0, 0)); assertEquals(2.toByte(), sites.founderAllele2(0, 0))
        assertEquals(3.toByte(), sites.founderAllele1(0, 1)); assertEquals(4.toByte(), sites.founderAllele2(0, 1))
        assertEquals(5.toByte(), sites.founderAllele1(1, 0)); assertEquals(6.toByte(), sites.founderAllele2(1, 0))
        assertEquals(8.toByte(), sites.founderAllele2(1, 1))
        assertEquals(9.toByte(), sites.sampleAllele1(0, 0)); assertEquals(12.toByte(), sites.sampleAllele2(0, 1))
        assertEquals(13.toByte(), sites.sampleAllele1(1, 0)); assertEquals(16.toByte(), sites.sampleAllele2(1, 1))
    }
}
