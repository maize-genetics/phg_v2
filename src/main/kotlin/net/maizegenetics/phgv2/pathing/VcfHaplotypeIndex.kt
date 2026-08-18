package net.maizegenetics.phgv2.pathing

import htsjdk.variant.variantcontext.Allele
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.getBufferedReader
import net.maizegenetics.phgv2.utils.getBufferedWriter
import java.io.BufferedWriter

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
}
