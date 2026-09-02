package net.maizegenetics.phgv2.pathing.ropebwt

import java.io.DataInputStream
import java.io.File
import java.nio.ByteBuffer
import java.nio.ByteOrder

/**
 * Windowed founder-to-founder sequence sharing, precomputed once per pangenome index from a
 * `ropebwt3 lift` file and read here as a lookup.
 *
 * Sharing is the local Jaccard overlap of two founders' anchor sets: of the reference anchors
 * present in either founder, the fraction present in both. Because it is derived from the index
 * rather than from any sample's reads, it is available for founders the sample does not carry --
 * which is exactly where a sample-derived estimate breaks down.
 *
 * Built by `scripts/build_sharing_table.py`. Binary, little-endian:
 *   magic "PSHARE" + version byte + pad, int32 windowSize, int32 nTaxa,
 *   per taxon (int32 nameLen, UTF-8 name), int32 nChrom,
 *   per chrom (int32 nameLen, UTF-8 name, int32 nWindows),
 *   then float32[nChrom][nWindows][nTaxa][nTaxa] row-major.
 */
class SharingTable(file: File) {
    val windowSize: Int
    val nTaxa: Int
    private val taxonIndex: Map<String, Int>
    private val chromWindows = LinkedHashMap<String, Int>()
    private val chromData = HashMap<String, FloatArray>()

    init {
        DataInputStream(file.inputStream().buffered(1 shl 20)).use { input ->
            val magic = ByteArray(8).also { input.readFully(it) }
            require(String(magic, Charsets.ISO_8859_1).startsWith("PSHARE")) {
                "${file.name} is not a sharing table"
            }
            fun int() = Integer.reverseBytes(input.readInt())
            fun str() = ByteArray(int()).also { input.readFully(it) }.toString(Charsets.UTF_8)
            windowSize = int()
            nTaxa = int()
            taxonIndex = (0 until nTaxa).associate { str() to it }
            repeat(int()) { chromWindows[str()] = int() }
            for ((chrom, windows) in chromWindows) {
                val cells = windows * nTaxa * nTaxa
                val bytes = ByteArray(cells * 4).also { input.readFully(it) }
                val buffer = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN).asFloatBuffer()
                chromData[chrom] = FloatArray(cells).also { buffer.get(it) }
            }
        }
    }

    fun hasContig(contig: String) = chromData.containsKey(contig)

    fun missingTaxa(names: Collection<String>) = names.filterNot { taxonIndex.containsKey(it) }

    fun windowCount(contig: String) = chromWindows[contig] ?: 0

    /** Sharing between two founders in [window] of [contig], or NaN if either is unknown. */
    fun sharing(contig: String, window: Int, taxonA: String, taxonB: String): Double {
        val data = chromData[contig] ?: return Double.NaN
        val a = taxonIndex[taxonA] ?: return Double.NaN
        val b = taxonIndex[taxonB] ?: return Double.NaN
        val windows = chromWindows[contig]!!
        val w = window.coerceIn(0, windows - 1)
        return data[(w * nTaxa + a) * nTaxa + b].toDouble()
    }

    fun windowOf(referencePosition: Long) = (referencePosition / windowSize).toInt()
}
