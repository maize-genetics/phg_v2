package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import org.junit.jupiter.api.Test
import kotlin.math.ln
import kotlin.test.assertEquals

class GameteSetEmissionProbabilityTest {

    private val probCorrect = 0.98
    private val lnPc = ln(probCorrect)
    private val lnHalfPc = ln(0.5 * probCorrect)
    private val lnPe = ln(1.0 - probCorrect)

    /** Parents 10, 20, 30 map to local indices 0, 1, 2. */
    private fun emission(vararg sets: Ps4gGameteSet) =
        GameteSetEmissionProbability(mapOf(100 to sets.toMutableList()), setOf(10, 20, 30), probCorrect)
            .getDiploidEmissionProbabilityArray(0)

    @Test
    fun eachCaseOfTheModelTable() {
        // gamete set {10, 20}: parents 0 and 1 present, parent 2 absent
        val p = emission(Ps4gGameteSet(intArrayOf(10, 20), 1))
        val n = 3
        // homozygous, founder present
        assertEquals(lnPc, p[0 * n + 0], 1e-12)
        assertEquals(lnPc, p[1 * n + 1], 1e-12)
        // homozygous, founder absent
        assertEquals(lnPe, p[2 * n + 2], 1e-12)
        // heterozygous, both founders in the gamete set: an A = B site
        assertEquals(lnPc, p[0 * n + 1], 1e-12)
        assertEquals(lnPc, p[1 * n + 0], 1e-12)
        // heterozygous, exactly one founder in the gamete set: an A != B site
        assertEquals(lnHalfPc, p[0 * n + 2], 1e-12)
        assertEquals(lnHalfPc, p[2 * n + 0], 1e-12)
        assertEquals(lnHalfPc, p[1 * n + 2], 1e-12)
    }

    @Test
    fun heterozygousStateWithNeitherFounderPresent() {
        // gamete set {30}: for the pair (0, 1) neither founder is present
        val p = emission(Ps4gGameteSet(intArrayOf(30), 1))
        assertEquals(lnPe, p[0 * 3 + 1], 1e-12)
        assertEquals(lnPc, p[2 * 3 + 2], 1e-12)
    }

    @Test
    fun readCountsMultiply() {
        val one = emission(Ps4gGameteSet(intArrayOf(10, 20), 1))
        val three = emission(Ps4gGameteSet(intArrayOf(10, 20), 3))
        for (i in one.indices) assertEquals(3 * one[i], three[i], 1e-12)
    }

    @Test
    fun gameteSetsInTheSameBinAdd() {
        val p = emission(Ps4gGameteSet(intArrayOf(10, 20), 1), Ps4gGameteSet(intArrayOf(30), 1))
        // pair (0,1): first set is an A = B site (lnPc), second has neither founder (lnPe)
        assertEquals(lnPc + lnPe, p[0 * 3 + 1], 1e-12)
        // pair (2,2): first set lacks parent 2 (lnPe), second contains it (lnPc)
        assertEquals(lnPe + lnPc, p[2 * 3 + 2], 1e-12)
    }

    @Test
    fun discriminationMatchesTheAnalyticPrediction() {
        val n = 3
        // an "ab" read cannot separate the heterozygous state from the homozygous one
        val ab = emission(Ps4gGameteSet(intArrayOf(10, 20), 1))
        assertEquals(0.0, ab[0 * n + 1] - ab[0 * n + 0], 1e-12)
        // an "a only" read favours the homozygous state by ln 2
        val aOnly = emission(Ps4gGameteSet(intArrayOf(10), 1))
        assertEquals(ln(0.5), aOnly[0 * n + 1] - aOnly[0 * n + 0], 1e-12)
        // a "b only" read favours the heterozygous state by ln(0.5 * pc / pe)
        val bOnly = emission(Ps4gGameteSet(intArrayOf(20), 1))
        assertEquals(ln(0.5 * probCorrect / (1.0 - probCorrect)), bOnly[0 * n + 1] - bOnly[0 * n + 0], 1e-12)
    }
}
