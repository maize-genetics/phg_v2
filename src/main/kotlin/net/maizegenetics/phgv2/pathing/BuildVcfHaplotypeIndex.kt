package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.headerCommand
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.utils.getBufferedWriter
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.TreeMap

/**
 * Resolves a contig + 1-based position to the [ReferenceRange] that contains it.
 *
 * Backed by one [TreeMap] per contig, keyed on the range start, built once from the graph.
 *
 * Splitting by contig first is deliberate: a single TreeMap<Position, ReferenceRange> across all
 * contigs (as used e.g. by Hvcf2Vcf.buildPositionToRefRangeMap) is unsound, because
 * [net.maizegenetics.phgv2.utils.Position.compareTo] is not a total order over mixed
 * numeric/non-numeric contig names - it strips a "chr" prefix and compares contigs numerically
 * when both parse as integers, else falls back to a plain string compare. With contigs
 * {"9", "10", "chrUn"}, "10" > "9" numerically but "10" < "chrUn" while "9" > "chrUn"
 * lexicographically, which is not transitive. Splitting by contig first means each per-contig
 * TreeMap is keyed on plain Ints, which is a genuine total order, so no cross-contig re-check is
 * needed after a lookup and the input VCF does not need to be sorted in any particular way.
 */
class RefRangeLocator(rangesByContig: Map<String, List<ReferenceRange>>) {

    constructor(graph: HaplotypeGraph) : this(graph.rangesByContig())

    private val rangesByStart: Map<String, TreeMap<Int, ReferenceRange>> =
        rangesByContig.mapValues { (_, ranges) ->
            TreeMap<Int, ReferenceRange>().apply { ranges.forEach { range -> put(range.start, range) } }
        }

    /**
     * Returns the ReferenceRange containing [contig]:[position] (1-based, inclusive), or null if
     * [contig] is not in the graph or no range covers [position].
     */
    fun findRange(contig: String, position: Int): ReferenceRange? {
        val range = rangesByStart[contig]?.floorEntry(position)?.value ?: return null
        return if (position <= range.end) range else null
    }

    /**
     * Returns adjacent range pairs (within a contig) that overlap. PHG reference ranges are
     * expected to never overlap; when they do, [findRange] resolves a position within the
     * overlap to the range with the greater start.
     */
    fun detectOverlaps(): List<Pair<ReferenceRange, ReferenceRange>> {
        return rangesByStart.values.flatMap { treeMap ->
            treeMap.values.toList().zipWithNext().filter { (a, b) -> b.start <= a.end }
        }
    }
}

/**
 * Caches [HaplotypeGraph.sampleGameteToHaplotypeId] for the most recently requested range.
 *
 * A panel VCF is coordinate-sorted in the overwhelming majority of cases, so consecutive records
 * almost always request the same range as the previous one; a one-slot cache captures nearly all
 * of that benefit at O(numberOfGametes) memory. Correctness does not depend on the VCF being
 * sorted - a cache miss simply costs a fresh (and correct) lookup.
 */
class HapIdsByRangeCache(private val graph: HaplotypeGraph) {

    private var cachedRange: ReferenceRange? = null
    private var cachedMap: Map<SampleGamete, String> = emptyMap()

    fun hapIds(range: ReferenceRange): Map<SampleGamete, String> {
        if (range != cachedRange) {
            cachedMap = graph.sampleGameteToHaplotypeId(range)
            cachedRange = range
        }
        return cachedMap
    }
}

/**
 * Running tallies for a single build-vcf-haplotype-index run.
 */
data class VcfHaplotypeIndexStats(
    var totalPanelRecords: Int = 0,
    var positionsIndexed: Int = 0,
    var positionsWithNoRefRange: Int = 0,
    var noCallAlleles: Int = 0,
    var gametesWithNoHaplotypeInRange: Int = 0,
    var gametesFromSamplesNotInGraph: Int = 0,
    var allelesWithNoHapIds: Int = 0,
    var rowsWritten: Int = 0,
    val samplesNotInGraph: MutableSet<String> = sortedSetOf()
)

/**
 * Builds a tab-delimited index mapping each allele observed in a reference-panel VCF to the PHG
 * haplotype ids that carry it.
 *
 * For every record in the panel VCF, the ReferenceRange covering its start position is located in
 * a HaplotypeGraph built from --hvcf-dir. Then, for each sample/gamete's allele at that position,
 * the haplotype id that sample/gamete carries in that range is looked up. Rows are grouped by
 * (position, allele) with a deduplicated, sorted list of the haplotype ids observed for that
 * allele.
 *
 * A record's span is not considered - only its start position is used to locate the reference
 * range, matching the existing convention in ConvertVcf2Ps4gFile. A multi-base record whose start
 * falls in one range but whose end reaches into the next is therefore attributed entirely to the
 * range containing its start.
 *
 * This is step 1 of a two-step VCF-to-readMapping pipeline: step 2 (not implemented here) will
 * read this index (via [VcfHaplotypeIndexUtils.readIndex]) together with a VCF to impute, and
 * emit one readMapping file per sample.
 */
class BuildVcfHaplotypeIndex : CliktCommand(
    help = "Build a tab-delimited index mapping each reference-panel VCF allele to the PHG haplotype ids that carry it."
) {

    private val myLogger = LogManager.getLogger(BuildVcfHaplotypeIndex::class.java)

    val hvcfDir by option(help = "Full path to a directory of hVCF files used to build the HaplotypeGraph. Required parameter.")
        .required()
        .validate { require(File(it).isDirectory) { "$it is not a valid directory." } }

    val panelVcf by option(help = "Full path to the reference-panel VCF file. Required parameter.")
        .required()
        .validate { require(File(it).isFile) { "$it is not a valid file" } }

    val outputFile by option(help = "Full path to the output index file. Ending the name in .gz will gzip the output. Required parameter.")
        .required()

    override fun run() {
        logCommand(this)

        myLogger.info("Building HaplotypeGraph from $hvcfDir")
        val graph = HaplotypeGraph(hvcfDir)

        val locator = RefRangeLocator(graph)
        val overlaps = locator.detectOverlaps()
        if (overlaps.isNotEmpty()) {
            myLogger.warn(
                "Found ${overlaps.size} overlapping reference range pair(s) in the HaplotypeGraph; " +
                    "a position falling in an overlap resolves to the range with the greater start. " +
                    "First few: ${overlaps.take(10)}"
            )
        }

        val stats = buildAndWriteIndex(graph, locator, panelVcf, outputFile, headerCommand(this))
        logSummary(graph, stats)
    }

    /**
     * Streams [panelVcfFile], writing one index row per (position, allele) with at least one
     * resolved haplotype id to [outputFileName]. Memory use is bounded by the size of a single
     * reference range's sample-to-haplotype map, not by the number of rows in the panel VCF.
     */
    fun buildAndWriteIndex(
        graph: HaplotypeGraph,
        locator: RefRangeLocator,
        panelVcfFile: String,
        outputFileName: String,
        commandLine: String
    ): VcfHaplotypeIndexStats {

        val graphSampleNames = graph.samples().toSet()
        val cache = HapIdsByRangeCache(graph)
        val stats = VcfHaplotypeIndexStats()

        getBufferedWriter(outputFileName).use { writer ->
            VcfHaplotypeIndexUtils.writeHeader(writer, commandLine, hvcfDir, panelVcfFile, graph.checksum)

            var lastPositionKey: Pair<String, Int>? = null
            VCFFileReader(File(panelVcfFile), false).use { reader ->
                reader.forEach { record ->
                    val key = record.contig to record.start
                    if (key == lastPositionKey) {
                        myLogger.warn(
                            "Duplicate panel VCF record at ${record.contig}:${record.start}; rows for " +
                                "this position will be merged when the index is read."
                        )
                    }
                    lastPositionKey = key

                    processRecord(record, locator, graphSampleNames, cache::hapIds, stats).forEach { entry ->
                        writer.write(VcfHaplotypeIndexUtils.formatEntry(entry))
                        writer.write("\n")
                        stats.rowsWritten++
                    }
                }
            }
        }

        return stats
    }

    /**
     * Converts one panel VCF record into index rows. Returns an empty list when no reference
     * range covers the record's start position, or when every observed allele at the position
     * resolves to zero haplotype ids. [stats] is mutated in place with counts for every edge
     * case encountered.
     *
     * [hapIdLookup] takes the range containing the record and returns the graph's
     * SampleGamete -> haplotype id map for that range (see [HapIdsByRangeCache]); it is passed as
     * a lambda so this function can be unit tested without constructing a HaplotypeGraph.
     */
    fun processRecord(
        record: VariantContext,
        locator: RefRangeLocator,
        graphSampleNames: Set<String>,
        hapIdLookup: (ReferenceRange) -> Map<SampleGamete, String>,
        stats: VcfHaplotypeIndexStats
    ): List<VcfHaplotypeIndexEntry> {

        stats.totalPanelRecords++

        val range = locator.findRange(record.contig, record.start)
        if (range == null) {
            stats.positionsWithNoRefRange++
            return emptyList()
        }
        stats.positionsIndexed++

        val gameteToHapId = hapIdLookup(range)

        //sortedMapOf/sortedSetOf give deterministic row order and deduplicated, sorted hapIds
        //without any extra work.
        val alleleToHapIds = sortedMapOf<String, MutableSet<String>>()
        val observedAlleles = sortedSetOf<String>()

        record.genotypes.forEach { genotype: Genotype ->
            genotype.alleles.forEachIndexed { alleleIndex, allele: Allele ->
                if (allele.isNoCall) {
                    stats.noCallAlleles++
                    return@forEachIndexed
                }

                val alleleKey = VcfHaplotypeIndexUtils.alleleKey(allele)
                observedAlleles.add(alleleKey)

                val sampleGamete = SampleGamete(genotype.sampleName, alleleIndex)
                val hapId = gameteToHapId[sampleGamete]
                if (hapId == null) {
                    if (sampleGamete.name in graphSampleNames) {
                        stats.gametesWithNoHaplotypeInRange++
                    } else {
                        stats.gametesFromSamplesNotInGraph++
                        stats.samplesNotInGraph.add(sampleGamete.name)
                    }
                    return@forEachIndexed
                }

                alleleToHapIds.getOrPut(alleleKey) { sortedSetOf() }.add(hapId)
            }
        }

        stats.allelesWithNoHapIds += (observedAlleles - alleleToHapIds.keys).size

        return alleleToHapIds.map { (allele, hapIds) ->
            VcfHaplotypeIndexEntry(record.contig, record.start, allele, range, hapIds.toList())
        }
    }

    fun logSummary(graph: HaplotypeGraph, stats: VcfHaplotypeIndexStats) {
        myLogger.info("build-vcf-haplotype-index summary:")
        myLogger.info("  panel VCF records read: ${stats.totalPanelRecords}")
        myLogger.info("  positions indexed: ${stats.positionsIndexed}")
        myLogger.info("  positions with no covering reference range: ${stats.positionsWithNoRefRange}")
        myLogger.info("  index rows written: ${stats.rowsWritten}")
        myLogger.info("  no-call alleles skipped: ${stats.noCallAlleles}")
        myLogger.info("  gametes with no haplotype in their range: ${stats.gametesWithNoHaplotypeInRange}")
        myLogger.info("  gametes from samples not in the graph: ${stats.gametesFromSamplesNotInGraph}")
        myLogger.info("  alleles observed but resolving to no haplotype: ${stats.allelesWithNoHapIds}")

        if (stats.samplesNotInGraph.isNotEmpty()) {
            myLogger.warn(
                "${stats.samplesNotInGraph.size} panel VCF sample(s) are not in the HaplotypeGraph and " +
                    "contributed no haplotype ids: ${stats.samplesNotInGraph.take(20)}"
            )
        }
        if (stats.totalPanelRecords > 0 && stats.positionsIndexed == 0) {
            myLogger.error(
                "No panel VCF position fell inside any reference range. This usually means the panel VCF " +
                    "and the hVCFs use different contig names. Graph contigs: ${graph.contigs}"
            )
        }
    }
}
