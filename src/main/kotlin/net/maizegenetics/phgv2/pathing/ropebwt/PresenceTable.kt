package net.maizegenetics.phgv2.pathing.ropebwt

import java.io.DataInputStream
import java.io.File
import java.nio.ByteBuffer
import java.nio.ByteOrder

/**
 * Per-founder anchor presence, precomputed once per pangenome index from a `ropebwt3 lift` file.
 *
 * For each reference window and each founder, this holds the fraction of that window's reference
 * anchor grid the founder possesses. A value near zero means the founder has no alignable sequence
 * there, so no read can support it -- whether through deletion, assembly gap, or loss of anchor
 * uniqueness. All three have the same consequence for imputation.
 *
 * This is deliberately not the same quantity as [SharingTable]. Sharing answers "how much sequence
 * do these two founders hold in common here", which conflates one founder being absent with both
 * being present but divergent. Only the first means reads cannot support that founder, and only the
 * first should trigger a presence/absence correction.
 *
 * Validated against independently called gVCF deletions of 5 kb or more: windows with presence at
 * or below 0.05 are genuinely deleted 87-89% of the time. Note the base rate -- NAM founders lack
 * roughly 37% of the reference at that scale, so absence is common rather than exceptional.
 *
 * Built by `scripts/build_presence_table.py`. Binary, little-endian:
 *   magic "PAVPRES" + version byte, int32 windowSize, int32 nTaxa,
 *   per taxon (int32 nameLen, UTF-8 name), int32 nChrom,
 *   per chrom (int32 nameLen, UTF-8 name, int32 nWindows),
 *   then float32[nChrom][nWindows][nTaxa].
 */
class PresenceTable(file: File) {
    val windowSize: Int
    val nTaxa: Int
    private val taxonIndex: Map<String, Int>
    private val chromWindows = LinkedHashMap<String, Int>()
    private val chromData = HashMap<String, FloatArray>()

    init {
        DataInputStream(file.inputStream().buffered(1 shl 20)).use { input ->
            val magic = ByteArray(8).also { input.readFully(it) }
            require(String(magic, Charsets.ISO_8859_1).startsWith("PAVPRES")) {
                "${file.name} is not a presence table"
            }
            fun int() = Integer.reverseBytes(input.readInt())
            fun str() = ByteArray(int()).also { input.readFully(it) }.toString(Charsets.UTF_8)
            windowSize = int()
            nTaxa = int()
            taxonIndex = (0 until nTaxa).associate { str() to it }
            repeat(int()) { chromWindows[str()] = int() }
            for ((chrom, windows) in chromWindows) {
                val cells = windows * nTaxa
                val bytes = ByteArray(cells * 4).also { input.readFully(it) }
                val buffer = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN).asFloatBuffer()
                chromData[chrom] = FloatArray(cells).also { buffer.get(it) }
            }
        }
    }

    fun hasContig(contig: String) = chromData.containsKey(contig)

    fun missingTaxa(names: Collection<String>) = names.filterNot { taxonIndex.containsKey(it) }

    fun windowCount(contig: String) = chromWindows[contig] ?: 0

    /** Fraction of [window]'s anchor grid that [taxon] possesses, or NaN if either is unknown. */
    fun presence(contig: String, window: Int, taxon: String): Double {
        val data = chromData[contig] ?: return Double.NaN
        val taxonSlot = taxonIndex[taxon] ?: return Double.NaN
        val windows = chromWindows[contig]!!
        return data[window.coerceIn(0, windows - 1) * nTaxa + taxonSlot].toDouble()
    }
}
