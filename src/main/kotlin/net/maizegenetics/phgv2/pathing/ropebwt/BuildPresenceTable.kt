package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.cli.logCommand
import org.apache.logging.log4j.LogManager
import java.io.BufferedInputStream
import java.io.DataOutputStream
import java.io.File
import java.io.InputStream
import java.nio.ByteBuffer
import java.nio.ByteOrder
import java.util.zip.GZIPInputStream

/**
 * Builds the per-founder anchor presence table that [GameteSetEmissionProbability] uses for its
 * presence/absence correction, from a `ropebwt3 lift` file.
 *
 * A pairwise sharing measure -- "how much sequence do founders i and j hold in common here" --
 * conflates two situations that are not alike: founder j is *absent* (a presence/absence variant),
 * or both are present but divergent. Only the first means no read can support j, and only the first
 * should trigger a correction. This table therefore records, for each reference window and each
 * founder, the fraction of that window's reference anchor grid the founder possesses. Near zero
 * means the founder has no alignable sequence there -- deletion, assembly gap, or loss of anchor
 * uniqueness alike, all of which have the same consequence for imputation.
 *
 * ## Inputs
 *
 * `--lift-file` is the binary map `ropebwt3 lift` writes, little-endian:
 *
 *     magic "LIFT" + version byte, int64 nSeq, int64 nPoints, int64 offset[nSeq + 1],
 *     then nPoints x (int64 carrierPos, int64 refPos, int32 carrierSeqId, int32 refSeqId)
 *
 * `offset` partitions the point array by carrier sequence: the points of sequence `i` are
 * `[offset[i], offset[i + 1])`. Points are stored in carrier order, so this reads them as one
 * forward stream rather than seeking.
 *
 * `--length-file` is the index's `.fmd.len.gz`, one name/length line per sequence, tab separated,
 * in sequence-id order. A sequence is named `taxon_contig`, so the taxon is the part before the
 * first underscore, and the reference's own chromosomes are those named `refPrefix_chr*`.
 *
 * ## Output
 *
 * Little-endian, and exactly what [PresenceTable] reads:
 *
 *     magic "PAVPRES" + version byte, int32 windowSize, int32 nTaxa,
 *     per taxon (int32 nameLen, UTF-8 name), int32 nChrom,
 *     per chrom (int32 nameLen, UTF-8 name, int32 nWindows),
 *     then float32[nChrom][nWindows][nTaxa]
 *
 * Contig names are written with the reference prefix stripped (`B73_chr1` becomes `chr1`) so they
 * match the contig names in a ps4g file. A reference chromosome carrying no lift points is written
 * with a window count of zero and contributes no float block.
 *
 * The reference taxon is a special case: the lift maps *other* founders onto the reference, so the
 * reference has no points of its own. It holds every anchor by construction and is written as 1.0
 * in every window that has any anchor at all.
 *
 * ## Memory
 *
 * Reference positions are held in memory, grouped by (chromosome, taxon), so that duplicates
 * arising from several carrier contigs of one founder hitting the same reference position can be
 * removed before counting. That is 8 bytes per point: roughly 0.6 GB for a maize index with 72
 * million points. Run with a heap to match.
 */
class BuildPresenceTable : CliktCommand(help = "Build a founder anchor-presence table from a ropebwt3 lift file") {

    private val myLogger = LogManager.getLogger(BuildPresenceTable::class.java)

    val liftFile by option(help = "The binary map written by `ropebwt3 lift`. Required parameter.")
        .required()
        .validate { require(File(it).exists()) { "$it is not a valid file" } }

    val lengthFile by option(help = "The index's gzipped sequence-length file (.fmd.len.gz), one " +
            "tab-separated name and length per sequence, in sequence-id order. Required parameter.")
        .required()
        .validate { require(File(it).exists()) { "$it is not a valid file" } }

    val referencePrefix by option(help = "The taxon name of the reference, used to find its " +
            "chromosomes in the length file as <reference-prefix>_chr*. Required parameter.")
        .required()

    val windowSize by option(help = "Window size in base pairs. Default = 50000")
        .int()
        .default(50_000)
        .validate { require(it > 0) { "window-size must be positive" } }

    val outputFile by option(help = "The presence table to write. Required parameter.")
        .required()

    override fun run() {
        logCommand(this)
        val blocks = buildBlocks(File(liftFile), File(lengthFile), referencePrefix, windowSize)
        writeTable(File(outputFile), windowSize, blocks.taxa, blocks.chroms)
        myLogger.info("Wrote $outputFile (${File(outputFile).length()} bytes)")
    }

    /** One reference chromosome's table: [presence] is [nWindows] x nTaxa, row major, or null when empty. */
    data class ChromBlock(val name: String, val nWindows: Int, val presence: FloatArray?)

    data class Blocks(val taxa: List<String>, val chroms: List<ChromBlock>)

    /**
     * Reads [lift] and returns the presence fractions per reference chromosome.
     *
     * Kept separate from [run] so it can be exercised against a synthetic lift file without going
     * through the CLI.
     */
    fun buildBlocks(lift: File, lengths: File, refPrefix: String, window: Int): Blocks {
        val names = ArrayList<String>()
        val seqLengths = ArrayList<Long>()
        GZIPInputStream(lengths.inputStream().buffered()).bufferedReader().forEachLine { line ->
            if (line.isBlank()) return@forEachLine
            val parts = line.split('\t')
            names.add(parts[0])
            seqLengths.add(parts[1].trim().toLong())
        }

        val taxa = names.map { it.substringBefore("_") }.distinct().sorted()
        val taxonOf = taxa.withIndex().associate { (index, name) -> name to index }
        val seqTaxon = IntArray(names.size) { taxonOf[names[it].substringBefore("_")]!! }

        // Reference chromosomes, ordered numerically by the digits after "chr" so chr2 precedes
        // chr10; anything non-numeric sorts last, matching the builder this replaces.
        val refChroms = names.indices
            .filter { names[it].startsWith(refPrefix + "_chr") }
            .sortedBy { names[it].substringAfter("chr").toIntOrNull() ?: 999 }
        val refSlot = HashMap<Int, Int>(refChroms.size * 2)
        refChroms.forEachIndexed { slot, seqId -> refSlot[seqId] = slot }

        // refPositions[chromSlot][taxon] accumulates every reference position a founder's anchors
        // land on; duplicates are removed per (chromosome, taxon) below.
        val refPositions = Array(refChroms.size) { Array(taxa.size) { LongList() } }

        BufferedInputStream(lift.inputStream(), 1 shl 22).use { input ->
            val header = ByteArray(21)
            readFully(input, header)
            require(String(header, 0, 4, Charsets.ISO_8859_1) == "LIFT") {
                "${lift.name} is not a ropebwt3 lift file"
            }
            val headerBuffer = ByteBuffer.wrap(header).order(ByteOrder.LITTLE_ENDIAN)
            val nSeq = headerBuffer.getLong(5)
            val nPoints = headerBuffer.getLong(13)
            require(nSeq == names.size.toLong()) {
                "lift has $nSeq sequences but ${lengths.name} lists ${names.size}"
            }

            val offsetBytes = ByteArray(8 * (nSeq + 1).toInt())
            readFully(input, offsetBytes)
            val offsets = ByteBuffer.wrap(offsetBytes).order(ByteOrder.LITTLE_ENDIAN).asLongBuffer()

            // Points are contiguous in carrier order, so walking the offset boundaries forward
            // gives each point's carrier sequence without seeking.
            var carrier = 0
            var carrierEnd = offsets.get(1)
            val chunk = ByteArray(POINT_BYTES * POINTS_PER_CHUNK)
            var pointIndex = 0L
            while (pointIndex < nPoints) {
                val want = minOf(POINTS_PER_CHUNK.toLong(), nPoints - pointIndex).toInt()
                readFully(input, chunk, want * POINT_BYTES)
                val buffer = ByteBuffer.wrap(chunk).order(ByteOrder.LITTLE_ENDIAN)
                for (i in 0 until want) {
                    while (pointIndex >= carrierEnd && carrier < nSeq - 1) {
                        carrier++
                        carrierEnd = offsets.get(carrier + 1)
                    }
                    val base = i * POINT_BYTES
                    val refPos = buffer.getLong(base + 8)
                    val refSeq = buffer.getInt(base + 20)
                    val slot = refSlot[refSeq]
                    if (slot != null) refPositions[slot][seqTaxon[carrier]].add(refPos)
                    pointIndex++
                }
            }
        }

        val refTaxonSlot = taxonOf[refPrefix]
        val blocks = refChroms.mapIndexed { slot, seqId ->
            val contig = names[seqId].substringAfter("_")
            val perTaxon = refPositions[slot].map { it.sortedDistinct() }
            val anchorCount = perTaxon.sumOf { it.size }
            if (anchorCount == 0) return@mapIndexed ChromBlock(contig, 0, null)

            val nWindows = (seqLengths[seqId] / window).toInt() + 1
            // The grid is every reference position any founder anchors to, deduplicated across
            // founders: the denominator each founder's own count is taken against.
            val grid = LongList(anchorCount).also { all -> perTaxon.forEach { all.addAll(it) } }
                .sortedDistinct()
            val total = binCounts(grid, window, nWindows)

            val presence = FloatArray(nWindows * taxa.size)
            for (taxon in taxa.indices) {
                if (taxon == refTaxonSlot) {
                    // The reference holds every anchor by construction: the lift maps the other
                    // founders onto it, so it carries no points of its own.
                    for (w in 0 until nWindows) {
                        presence[w * taxa.size + taxon] = if (total[w] > 0) 1.0f else 0.0f
                    }
                    continue
                }
                val positions = perTaxon[taxon]
                if (positions.isEmpty()) continue
                val counts = binCounts(positions, window, nWindows)
                for (w in 0 until nWindows) {
                    if (total[w] > 0) presence[w * taxa.size + taxon] = counts[w].toFloat() / total[w]
                }
            }
            myLogger.info("${names[seqId]}: ${grid.size} anchors, $nWindows windows")
            ChromBlock(contig, nWindows, presence)
        }
        myLogger.info("taxa=${taxa.size} window=$window")
        return Blocks(taxa, blocks)
    }

    /** Counts how many of [sorted] fall in each of [nWindows] windows of [window] bases. */
    private fun binCounts(sorted: LongArray, window: Int, nWindows: Int): IntArray {
        val counts = IntArray(nWindows)
        for (position in sorted) {
            val w = (position / window).toInt()
            if (w in 0 until nWindows) counts[w]++
        }
        return counts
    }

    fun writeTable(out: File, window: Int, taxa: List<String>, chroms: List<ChromBlock>) {
        DataOutputStream(out.outputStream().buffered(1 shl 20)).use { output ->
            output.write("PAVPRES".toByteArray(Charsets.ISO_8859_1))
            output.write(1)
            fun int(value: Int) = output.writeInt(Integer.reverseBytes(value))
            fun str(value: String) = value.toByteArray(Charsets.UTF_8).let { int(it.size); output.write(it) }
            int(window)
            int(taxa.size)
            taxa.forEach { str(it) }
            int(chroms.size)
            chroms.forEach { str(it.name); int(it.nWindows) }
            for (chrom in chroms) {
                val presence = chrom.presence ?: continue
                val bytes = ByteArray(presence.size * 4)
                ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN).asFloatBuffer().put(presence)
                output.write(bytes)
            }
        }
    }

    private fun readFully(input: InputStream, buffer: ByteArray, length: Int = buffer.size) {
        var read = 0
        while (read < length) {
            val n = input.read(buffer, read, length - read)
            require(n > 0) { "unexpected end of lift file" }
            read += n
        }
    }

    /** A growable primitive long array: the point arrays are far too large to box. */
    class LongList(capacity: Int = 16) {
        private var values = LongArray(maxOf(capacity, 4))
        var size = 0
            private set

        fun add(value: Long) {
            if (size == values.size) values = values.copyOf(values.size * 2)
            values[size++] = value
        }

        fun addAll(other: LongArray) {
            var needed = values.size
            while (needed < size + other.size) needed *= 2
            if (needed != values.size) values = values.copyOf(needed)
            other.copyInto(values, size)
            size += other.size
        }

        /** Sorted copy with duplicates removed. */
        fun sortedDistinct(): LongArray {
            if (size == 0) return LongArray(0)
            val copy = values.copyOf(size)
            copy.sort()
            var kept = 1
            for (i in 1 until copy.size) if (copy[i] != copy[kept - 1]) copy[kept++] = copy[i]
            return copy.copyOf(kept)
        }
    }

    companion object {
        private const val POINT_BYTES = 24
        private const val POINTS_PER_CHUNK = 1 shl 16
    }
}
