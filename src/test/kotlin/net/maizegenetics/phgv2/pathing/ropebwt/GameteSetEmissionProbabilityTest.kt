package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import org.junit.jupiter.api.Test
import kotlin.math.ln
import kotlin.test.assertEquals

class GameteSetEmissionProbabilityTest {

    /** Writes a presence table holding one chromosome and the given per-taxon presence values. */
    private fun presenceTable(vararg presence: Pair<String, Float>): PresenceTable {
        val file = java.io.File.createTempFile("presence", ".bin").also { it.deleteOnExit() }
        java.io.DataOutputStream(file.outputStream().buffered()).use { out ->
            fun int(v: Int) = out.writeInt(Integer.reverseBytes(v))
            fun str(s: String) { int(s.length); out.write(s.toByteArray()) }
            out.write("PAVPRES".toByteArray()); out.write(1)
            int(50_000); int(presence.size)
            presence.forEach { str(it.first) }
            int(1); str("chr1"); int(1)
            val buf = java.nio.ByteBuffer.allocate(presence.size * 4)
                .order(java.nio.ByteOrder.LITTLE_ENDIAN)
            presence.forEach { buf.putFloat(it.second) }
            out.write(buf.array())
        }
        return PresenceTable(file)
    }


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

    /**
     * With founder P2 absent, state (P0, P2) must predict exactly what state (P0, P0) predicts:
     * the two tie, so the transition prior decides rather than the absent founder's missing reads
     * being read as evidence for homozygosity.
     */
    @Test
    fun anAbsentFounderTiesTheHeterozygousStateWithTheHomozygousOne() {
        val table = presenceTable("P0" to 0.9f, "P1" to 0.9f, "P2" to 0.0f)
        val names = mapOf(10 to "P0", 20 to "P1", 30 to "P2")
        // gamete set {P0}: without the correction, (0,2) would score ln(0.5*pc) as an A != B site
        val p = GameteSetEmissionProbability(
            mapOf(100 to mutableListOf(Ps4gGameteSet(intArrayOf(10), 1))),
            setOf(10, 20, 30), probCorrect, table, "chr1", names, 256, 0.05
        ).getDiploidEmissionProbabilityArray(0)
        val n = 3
        assertEquals(lnPc, p[0 * n + 0], 1e-12)
        assertEquals(lnPc, p[0 * n + 2], 1e-12)   // tied with the homozygous state, not penalised
        assertEquals(lnPc, p[2 * n + 0], 1e-12)   // and symmetric in the ordered pair
        // the pair not involving the absent founder keeps the ordinary rule
        assertEquals(lnHalfPc, p[0 * n + 1], 1e-12)
        // the absent founder's own homozygous state is still scored badly, correctly
        assertEquals(lnPe, p[2 * n + 2], 1e-12)
    }

    /**
     * The regression this test exists for: with damping 0 the correction must reproduce the
     * uncorrected model EXACTLY, even where a founder is flagged absent. An earlier version
     * failed here, because it also rewrote the case where the flagged founder is in the gamete
     * set but the other is not, which is common -- absence is a window-level call and a typical
     * read matches many founders regardless.
     */
    @Test
    fun zeroDampingReproducesTheUncorrectedModelEvenWhenFoundersAreFlagged() {
        val table = presenceTable("P0" to 0.0f, "P1" to 0.9f, "P2" to 0.0f)
        val names = mapOf(10 to "P0", 20 to "P1", 30 to "P2")
        for (gs in listOf(intArrayOf(10), intArrayOf(20), intArrayOf(30),
                          intArrayOf(10, 20), intArrayOf(20, 30), intArrayOf(10, 30))) {
            val damped = GameteSetEmissionProbability(
                mapOf(100 to mutableListOf(Ps4gGameteSet(gs, 1))),
                setOf(10, 20, 30), probCorrect, table, "chr1", names, 256, 0.05, 0.0
            ).getDiploidEmissionProbabilityArray(0)
            val plain = emission(Ps4gGameteSet(gs, 1))
            for (i in plain.indices) assertEquals(plain[i], damped[i], 1e-12,
                "gamete set ${gs.toList()} state index $i")
        }
    }

    @Test
    fun anAbsentFounderOnlyChangesTheOneSidedCase() {
        val table = presenceTable("P0" to 0.9f, "P1" to 0.9f, "P2" to 0.0f)
        val names = mapOf(10 to "P0", 20 to "P1", 30 to "P2")
        // gamete set {P0, P2}: P2 is flagged absent but IS present in G, so the pair (0,2) is an
        // ordinary A = B site and must be untouched
        val p = GameteSetEmissionProbability(
            mapOf(100 to mutableListOf(Ps4gGameteSet(intArrayOf(10, 30), 1))),
            setOf(10, 20, 30), probCorrect, table, "chr1", names, 256, 0.05, 1.0
        ).getDiploidEmissionProbabilityArray(0)
        assertEquals(lnPc, p[0 * 3 + 2], 1e-12)
        // and (1,2): P1 not in G, P2 in G and flagged absent -> still an ordinary divergent site
        assertEquals(lnHalfPc, p[1 * 3 + 2], 1e-12)
    }

    @Test
    fun aPresentFounderIsUnaffectedByTheTable() {
        val table = presenceTable("P0" to 0.9f, "P1" to 0.9f, "P2" to 0.9f)
        val names = mapOf(10 to "P0", 20 to "P1", 30 to "P2")
        val withTable = GameteSetEmissionProbability(
            mapOf(100 to mutableListOf(Ps4gGameteSet(intArrayOf(10), 1))),
            setOf(10, 20, 30), probCorrect, table, "chr1", names, 256, 0.05
        ).getDiploidEmissionProbabilityArray(0)
        val without = emission(Ps4gGameteSet(intArrayOf(10), 1))
        for (i in without.indices) assertEquals(without[i], withTable[i], 1e-12)
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
