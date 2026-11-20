package net.maizegenetics.phgv2.pathing

import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.getBufferedReader
import net.maizegenetics.phgv2.utils.getBufferedWriter
import java.io.BufferedWriter
import java.io.Closeable
import java.io.File
import java.io.FileInputStream
import java.io.InputStream
import java.util.Arrays
import java.util.zip.GZIPInputStream

/**
 * One row of a VCF haplotype index: at [contig]:[position], the reference-panel allele
 * [allele] is carried by exactly the haplotypes in [hapIds] within [refRange]. [hapIds]
 * is expected to be deduplicated and sorted.
 *
 * Note: [position] and [refRange] both refer to the start of the underlying VCF record.
 * A record whose span crosses a reference range boundary (e.g. a multi-base indel) is
 * attributed entirely to the range containing its start position - it is not split
 * across ranges.
 */
data class VcfHaplotypeIndexEntry(
    val contig: String,
    val position: Int,
    val allele: String,
    val refRange: ReferenceRange,
    val hapIds: List<String>
) {
    val positionKey: Position get() = Position(contig, position)
}

/**
 * All the rows for a single panel position, collapsed into an allele -> hapIds lookup.
 * This is the shape a downstream consumer (e.g. a readMapping builder) wants: given an
 * observed allele at a position, get the haplotype ids it implies.
 */
data class VcfHaplotypeIndexPosition(
    val position: Position,
    val refRange: ReferenceRange,
    val alleleToHapIds: Map<String, List<String>>
)

/**
 * Reader/writer and shared allele-keying logic for the VCF haplotype index produced by
 * [BuildVcfHaplotypeIndex]. Any code that consumes this index (e.g. a future readMapping
 * builder) should use [alleleKey] on its own VCF alleles and [readIndex] to load the file -
 * that is the entire contract between producer and consumer.
 */
object VcfHaplotypeIndexUtils {

    const val MAGIC = "#PHG_VCF_HAPLOTYPE_INDEX"
    const val VERSION = "1.0"
    const val COLUMN_HEADER = "#chrom\tpos\tallele\trefRange\thapIds"

    /**
     * Returns the string used to key an [allele] in the index. Uses [Allele.getDisplayString],
     * which preserves symbolic alleles (e.g. "<DEL>") verbatim. Do NOT use
     * [Allele.getBaseString] here - htsjdk returns an empty string from getBaseString() for
     * symbolic alleles (it only returns actual bases, and symbolic alleles have none), which
     * would silently collapse every symbolic allele at a site into a single "" key.
     *
     * Throws [IllegalArgumentException] for a NO_CALL allele - callers must filter those out
     * before computing an allele key, since a no-call carries no allele information to index.
     */
    fun alleleKey(allele: Allele): String {
        require(!allele.isNoCall) { "alleleKey: NO_CALL alleles must be filtered out before calling alleleKey" }
        return allele.displayString
    }

    /**
     * Formats a single index row as a tab-delimited line (no trailing newline).
     */
    fun formatEntry(entry: VcfHaplotypeIndexEntry): String {
        return "${entry.contig}\t${entry.position}\t${entry.allele}\t${entry.refRange}\t${entry.hapIds.joinToString(",")}"
    }

    /**
     * Parses a single tab-delimited data row (not a header/comment line) into a
     * [VcfHaplotypeIndexEntry]. Throws [IllegalArgumentException] if the line does not have
     * exactly 5 tab-delimited fields.
     */
    fun parseEntry(line: String): VcfHaplotypeIndexEntry {
        val fields = line.split("\t")
        require(fields.size == 5) { "Malformed VCF haplotype index line (expected 5 tab-delimited fields, got ${fields.size}): $line" }
        return VcfHaplotypeIndexEntry(
            contig = fields[0],
            position = fields[1].toInt(),
            allele = fields[2],
            refRange = parseRefRange(fields[3]),
            hapIds = fields[4].split(",")
        )
    }

    /**
     * Parses a ReferenceRange formatted as "contig:start-end" (i.e. [ReferenceRange.toString]).
     *
     * Deliberately does NOT use [ReferenceRange.parse]: that function splits on the first "-"
     * in the string, which mangles contig names that themselves contain a dash (e.g.
     * "scaffold-1:100-200"). This parser splits on the LAST ":" and the LAST "-" instead, so
     * dashes and colons inside the contig name are preserved correctly.
     */
    fun parseRefRange(field: String): ReferenceRange {
        val contig = field.substringBeforeLast(":")
        val coords = field.substringAfterLast(":")
        val start = coords.substringBeforeLast("-")
        val end = coords.substringAfterLast("-")
        require(contig.isNotEmpty() && start.isNotEmpty() && start != coords) { "Malformed refRange field: $field" }
        return ReferenceRange(contig, start.toInt(), end.toInt())
    }

    /**
     * Writes the comment-line header block for a VCF haplotype index file, ending with the
     * [COLUMN_HEADER] line. [graphChecksum] is optional (pass "" to omit it) but recommended:
     * it lets a downstream reader confirm the index was built from the expected HaplotypeGraph
     * before trusting its haplotype ids.
     */
    fun writeHeader(
        writer: BufferedWriter,
        commandLine: String = "",
        hvcfDir: String = "",
        panelVcf: String = "",
        graphChecksum: String = ""
    ) {
        writer.write("$MAGIC\n")
        writer.write("#version=$VERSION\n")
        if (commandLine.isNotEmpty()) writer.write("#PHG Command: $commandLine\n")
        if (hvcfDir.isNotEmpty()) writer.write("#hvcfDir=$hvcfDir\n")
        if (panelVcf.isNotEmpty()) writer.write("#panelVcf=$panelVcf\n")
        if (graphChecksum.isNotEmpty()) writer.write("#graphChecksum=$graphChecksum\n")
        writer.write("$COLUMN_HEADER\n")
    }

    /**
     * Writes a full VCF haplotype index file: header block followed by one line per entry,
     * in [entries] order. If [outputFile] ends in ".gz" the file is gzip-compressed.
     */
    fun writeIndex(
        outputFile: String,
        entries: List<VcfHaplotypeIndexEntry>,
        commandLine: String = "",
        hvcfDir: String = "",
        panelVcf: String = "",
        graphChecksum: String = ""
    ) {
        getBufferedWriter(outputFile).use { writer ->
            writeHeader(writer, commandLine, hvcfDir, panelVcf, graphChecksum)
            entries.forEach { entry ->
                writer.write(formatEntry(entry))
                writer.write("\n")
            }
        }
    }

    /**
     * Reads every data row of a VCF haplotype index file, in file order. Lines starting with
     * "#" are treated as header/comment lines and skipped, except that the exact
     * [COLUMN_HEADER] line must appear before any data row - this is what distinguishes a
     * valid index file from an arbitrary tab-delimited file that merely starts with "#" lines.
     *
     * Throws [IllegalArgumentException] if [COLUMN_HEADER] is never seen, or if a data row is
     * malformed.
     */
    fun readEntries(inputFile: String): List<VcfHaplotypeIndexEntry> {
        val entries = mutableListOf<VcfHaplotypeIndexEntry>()
        var sawColumnHeader = false
        getBufferedReader(inputFile).useLines { lines ->
            lines.forEach { line ->
                when {
                    line.isBlank() -> {}
                    line == COLUMN_HEADER -> sawColumnHeader = true
                    line.startsWith("#") -> {}
                    else -> {
                        require(sawColumnHeader) {
                            "$inputFile is not a valid VCF haplotype index: found a data row before the '$COLUMN_HEADER' line"
                        }
                        entries.add(parseEntry(line))
                    }
                }
            }
        }
        require(sawColumnHeader) { "$inputFile is not a valid VCF haplotype index: missing the '$COLUMN_HEADER' column header line" }
        return entries
    }

    /**
     * Reads a VCF haplotype index file and groups its rows by position, in the shape a
     * downstream consumer wants: for a given [Position], the allele -> hapIds lookup needed
     * to score an observed allele.
     *
     * WARNING: this materializes the ENTIRE index in memory, at several times its on-disk size.
     * Only use it for small indexes (tests, debugging). Production consumers should stream the
     * index with [VcfHaplotypeIndexStreamReader] instead.
     *
     * If the index file contains more than one row for the same (position, allele) - e.g. from
     * duplicate records in the source panel VCF - their hapIds are unioned rather than the
     * later row silently overwriting the earlier one.
     */
    fun readIndex(inputFile: String): Map<Position, VcfHaplotypeIndexPosition> {
        val alleleHapIdsByPosition = LinkedHashMap<Position, MutableMap<String, MutableSet<String>>>()
        val rangeByPosition = HashMap<Position, ReferenceRange>()

        readEntries(inputFile).forEach { entry ->
            val key = entry.positionKey
            rangeByPosition.putIfAbsent(key, entry.refRange)
            alleleHapIdsByPosition
                .getOrPut(key) { linkedMapOf() }
                .getOrPut(entry.allele) { sortedSetOf() }
                .addAll(entry.hapIds)
        }

        return alleleHapIdsByPosition.mapValues { (position, alleleMap) ->
            VcfHaplotypeIndexPosition(
                position = position,
                refRange = rangeByPosition.getValue(position),
                alleleToHapIds = alleleMap.mapValues { it.value.toList() }
            )
        }
    }

    /**
     * Returns the "#graphChecksum=..." value recorded in the index's header block, or null if
     * the index was written without one.
     */
    fun readGraphChecksum(inputFile: String): String? {
        getBufferedReader(inputFile).useLines { lines ->
            for (line in lines) {
                if (line == COLUMN_HEADER) return null
                if (line.startsWith("#graphChecksum=")) return line.substringAfter("#graphChecksum=")
            }
        }
        return null
    }

    /**
     * Makes one fast pass over a "#"-headed, tab-delimited text file whose first column is a
     * contig name (a VCF, or a VCF haplotype index) and returns its contigs in the order their
     * records appear. Works on plain text or gzip/bgzip (".gz") files; ".bcf" files are read
     * through htsjdk instead.
     *
     * Only the first column of each line is ever decoded - the rest of the line is skipped at
     * the byte level - so this is cheap even for very wide VCFs. Throws
     * [IllegalArgumentException] if a contig's records are not contiguous (e.g. chr1, chr2, chr1),
     * since no single contig order could then describe the file.
     */
    fun scanContigOrder(inputFile: String): List<String> {
        if (inputFile.endsWith(".bcf")) return scanContigOrderWithHtsjdk(inputFile)

        val order = mutableListOf<String>()
        val seen = HashSet<String>()

        var field = ByteArray(256)
        var fieldLen = 0
        var last = ByteArray(0)

        fun finishField() {
            if (fieldLen == last.size && Arrays.equals(field, 0, fieldLen, last, 0, last.size)) return
            val contig = String(field, 0, fieldLen, Charsets.UTF_8)
            require(seen.add(contig)) {
                "$inputFile is not sorted: records for contig '$contig' appear in more than one " +
                    "non-contiguous block. Sort the file so each contig's records are together."
            }
            order.add(contig)
            last = field.copyOf(fieldLen)
        }

        //0 = at the start of a line, 1 = reading the first (contig) column, 2 = skipping to end of line
        var state = 0
        val fileStream = FileInputStream(inputFile)
        val stream: InputStream = if (inputFile.endsWith(".gz")) GZIPInputStream(fileStream, 1 shl 16) else fileStream
        stream.use { input ->
            val buffer = ByteArray(1 shl 16)
            while (true) {
                val n = input.read(buffer)
                if (n < 0) break
                for (i in 0 until n) {
                    val b = buffer[i]
                    when (state) {
                        0 -> when (b) {
                            HASH, NEWLINE_LF, NEWLINE_CR -> if (b == HASH) state = 2
                            else -> {
                                fieldLen = 0
                                field[fieldLen++] = b
                                state = 1
                            }
                        }
                        1 -> when (b) {
                            TAB -> { finishField(); state = 2 }
                            NEWLINE_LF -> { finishField(); state = 0 }
                            NEWLINE_CR -> {}
                            else -> {
                                if (fieldLen == field.size) field = field.copyOf(field.size * 2)
                                field[fieldLen++] = b
                            }
                        }
                        else -> if (b == NEWLINE_LF) state = 0
                    }
                }
            }
        }
        if (state == 1) finishField()

        return order
    }

    private fun scanContigOrderWithHtsjdk(inputFile: String): List<String> {
        val order = mutableListOf<String>()
        val seen = HashSet<String>()
        VCFFileReader(File(inputFile), false).use { reader ->
            reader.forEach { record ->
                if (order.isNotEmpty() && order.last() == record.contig) return@forEach
                require(seen.add(record.contig)) {
                    "$inputFile is not sorted: records for contig '${record.contig}' appear in more " +
                        "than one non-contiguous block. Sort the file so each contig's records are together."
                }
                order.add(record.contig)
            }
        }
        return order
    }

    private const val HASH = '#'.code.toByte()
    private const val TAB = '\t'.code.toByte()
    private const val NEWLINE_LF = '\n'.code.toByte()
    private const val NEWLINE_CR = '\r'.code.toByte()

    /**
     * Verifies that the contigs shared by the index and another file (e.g. the impute VCF) appear
     * in the SAME relative order in both. Contigs present in only one file are ignored - they are
     * simply skipped while streaming. Different tools sort chromosomes differently (lexicographic
     * "1,10,2" vs natural "1,2,10", "chr" prefixes, scaffold placement), and a streaming merge of
     * two files sorted differently would silently drop evidence, so this fails loudly instead.
     */
    fun requireMatchingContigOrder(
        indexContigOrder: List<String>,
        otherContigOrder: List<String>,
        indexFile: String,
        otherFile: String
    ) {
        val indexContigs = indexContigOrder.toHashSet()
        val otherContigs = otherContigOrder.toHashSet()
        val sharedInIndexOrder = indexContigOrder.filter { it in otherContigs }
        val sharedInOtherOrder = otherContigOrder.filter { it in indexContigs }
        if (sharedInIndexOrder == sharedInOtherOrder) return

        val firstDiff = sharedInIndexOrder.indices.first { sharedInIndexOrder[it] != sharedInOtherOrder[it] }
        fun preview(contigs: List<String>) =
            if (contigs.size <= 40) contigs.toString() else "${contigs.take(40)} ... (${contigs.size} total)"

        throw IllegalArgumentException(
            "Contig order differs between the index ($indexFile) and $otherFile. Both files must list " +
                "their shared contigs in the same order so they can be streamed together.\n" +
                "  First difference at shared contig #${firstDiff + 1}: index has '${sharedInIndexOrder[firstDiff]}', " +
                "$otherFile has '${sharedInOtherOrder[firstDiff]}'.\n" +
                "  index order:  ${preview(sharedInIndexOrder)}\n" +
                "  $otherFile order: ${preview(sharedInOtherOrder)}\n" +
                "Re-sort $otherFile so its contigs follow the index's order, or rebuild the index from a " +
                "panel VCF sorted the same way as $otherFile."
        )
    }
}

/**
 * Streams a VCF haplotype index one position at a time, so a consumer never holds more than a
 * single position's rows in memory (contrast [VcfHaplotypeIndexUtils.readIndex], which loads it all).
 *
 * Consecutive rows for the same position are merged exactly like [VcfHaplotypeIndexUtils.readIndex]
 * does (hapIds unioned per allele). hapId strings are pooled per reference range, so every
 * [VcfHaplotypeIndexPosition] in a range shares the same String instances.
 *
 * Throws [IllegalArgumentException] if positions within a contig ever go backwards, or a contig's
 * rows are not contiguous - the streaming merge in [BuildReadMappingFromVcf] depends on both.
 */
class VcfHaplotypeIndexStreamReader(private val inputFile: String) : Closeable {

    private val reader = getBufferedReader(inputFile, 1 shl 16)

    private var pending: VcfHaplotypeIndexEntry? = null

    private val seenContigs = HashSet<String>()
    private var lastContig: String? = null
    private var lastPosition = Int.MIN_VALUE

    private var currentRange: ReferenceRange? = null
    private val hapIdPool = HashMap<String, String>()

    init {
        try {
            var sawColumnHeader = false
            while (!sawColumnHeader) {
                val line = reader.readLine() ?: break
                when {
                    line.isBlank() -> {}
                    line == VcfHaplotypeIndexUtils.COLUMN_HEADER -> sawColumnHeader = true
                    line.startsWith("#") -> {}
                    else -> throw IllegalArgumentException(
                        "$inputFile is not a valid VCF haplotype index: found a data row before the " +
                            "'${VcfHaplotypeIndexUtils.COLUMN_HEADER}' line"
                    )
                }
            }
            require(sawColumnHeader) {
                "$inputFile is not a valid VCF haplotype index: missing the '${VcfHaplotypeIndexUtils.COLUMN_HEADER}' column header line"
            }
            pending = readNextEntry()
        } catch (e: Exception) {
            reader.close()
            throw e
        }
    }

    /**
     * Returns the next position in file order, with all of its rows merged, or null at end of file.
     */
    fun next(): VcfHaplotypeIndexPosition? {
        val first = pending ?: return null
        val alleleToHapIds = LinkedHashMap<String, List<String>>()
        alleleToHapIds[first.allele] = first.hapIds

        var nextEntry = readNextEntry()
        while (nextEntry != null && nextEntry.position == first.position && nextEntry.contig == first.contig) {
            val existing = alleleToHapIds[nextEntry.allele]
            alleleToHapIds[nextEntry.allele] =
                if (existing == null) nextEntry.hapIds else (existing + nextEntry.hapIds).toSortedSet().toList()
            nextEntry = readNextEntry()
        }
        pending = nextEntry

        return VcfHaplotypeIndexPosition(first.positionKey, first.refRange, alleleToHapIds)
    }

    private fun readNextEntry(): VcfHaplotypeIndexEntry? {
        while (true) {
            val line = reader.readLine() ?: return null
            if (line.isBlank() || line.startsWith("#")) continue

            val raw = VcfHaplotypeIndexUtils.parseEntry(line)
            val contig = checkSorted(raw.contig, raw.position)

            if (raw.refRange != currentRange) {
                currentRange = raw.refRange
                hapIdPool.clear()
            }
            return raw.copy(
                contig = contig,
                refRange = currentRange!!,
                hapIds = sortedDistinct(raw.hapIds).map { hapIdPool.getOrPut(it) { it } }
            )
        }
    }

    /**
     * Validates sort order and returns a canonical (shared) instance of [contig].
     */
    private fun checkSorted(contig: String, position: Int): String {
        val previous = lastContig
        if (previous != null && previous == contig) {
            require(position >= lastPosition) {
                "$inputFile is not sorted by position: $contig:$position appears after $contig:$lastPosition"
            }
            lastPosition = position
            return previous
        }
        require(seenContigs.add(contig)) {
            "$inputFile is not sorted: rows for contig '$contig' appear in more than one non-contiguous block"
        }
        lastContig = contig
        lastPosition = position
        return contig
    }

    /**
     * The index builder already writes hapIds sorted and deduplicated; only pay for sorting when
     * a row violates that, so hapId-set equality matches [VcfHaplotypeIndexUtils.readIndex].
     */
    private fun sortedDistinct(hapIds: List<String>): List<String> {
        for (i in 1 until hapIds.size) {
            if (hapIds[i - 1] >= hapIds[i]) return hapIds.toSortedSet().toList()
        }
        return hapIds
    }

    override fun close() = reader.close()
}
