package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.long
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap
import it.unimi.dsi.fastutil.objects.Reference2IntOpenHashMap
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.utils.getBufferedWriter
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * Running tallies for a single build-read-mapping-from-vcf run. Counters are Longs: with
 * thousands of samples and millions of sites, per-allele counts easily exceed Int.MAX_VALUE.
 */
data class ReadMappingStats(
    var totalImputeRecords: Long = 0,
    var positionsMatched: Long = 0,
    var positionsNotInIndex: Long = 0,
    var noCallAlleles: Long = 0,
    var allelesNotInIndex: Long = 0,
    var gameteHits: Long = 0,
    var totalEvidenceUnits: Long = 0,
    var adFallbackCount: Long = 0,
    var zeroDepthAlleleCalls: Long = 0,
    var flushes: Long = 0,
    val samplesWithNoEvidence: MutableList<String> = mutableListOf()
)

/**
 * Builds one readMapping file per sample from a VCF, using the index produced by
 * [BuildVcfHaplotypeIndex] to translate each of the VCF's observed alleles into the set of PHG
 * haplotype ids that carry it.
 *
 * This is step 2 of the VCF-to-readMapping pipeline. Step 1 ([BuildVcfHaplotypeIndex]) does the
 * expensive HaplotypeGraph work once, up front; this step needs only the flat index file it
 * produces - no HaplotypeGraph is built or needed here.
 *
 * For every genotype call in --impute-vcf, each called allele (one per gamete - e.g. 2 for a
 * diploid site) is treated as one independent piece of evidence, exactly like a single sequencing
 * read: it contributes a +1 (or, with --use-allele-depth, its AD value) to the count for whichever
 * set of haplotypes the index says carries that allele at that position. Both alleles of a
 * heterozygous call land in the same per-sample readMapping map, since - like a real read - a
 * called allele doesn't reveal which parental chromosome it came from.
 *
 * Memory: the index and the impute VCF are streamed together in a single sorted merge pass, so
 * neither is ever held in memory. Because every hapId set belongs to exactly one reference range
 * and both files are sorted, a sample's counts for a range are final once the merge moves past
 * that range; counts are buffered in compact int maps and appended to the per-sample files at a
 * range boundary whenever more than --max-buffered-entries (sample, hapId set) counts are held.
 * Peak memory is therefore set by that threshold, not by index size, sample count, or genome size.
 *
 * The merge requires both files to list their shared contigs in the same order. Before any
 * evidence is counted, each file gets one fast first-column-only pass to read its contig order,
 * and the run fails with a descriptive error if the orders disagree.
 *
 * Output files are gzipped and named <sampleName>_readMapping.txt.gz (a suffix FindPaths accepts), one per
 * sample in --impute-vcf - including samples with zero resolvable evidence, so a keyfile-driven
 * downstream step can assume every sample's file exists. A pathKeyFile.txt is also written to
 * --output-dir so it can be passed directly to `find-paths --path-keyfile`.
 */
class BuildReadMappingFromVcf : CliktCommand(
    help = "Build one readMapping file per sample from a VCF, using a VCF haplotype index produced by build-vcf-haplotype-index."
) {

    private val myLogger = LogManager.getLogger(BuildReadMappingFromVcf::class.java)

    val indexFile by option(help = "Full path to the index file produced by build-vcf-haplotype-index. Required parameter.")
        .required()
        .validate { require(File(it).isFile) { "$it is not a valid file" } }

    val imputeVcf by option(help = "Full path to the VCF file to impute readMapping evidence from. Its contigs must be in the same order as the index's. Required parameter.")
        .required()
        .validate { require(File(it).isFile) { "$it is not a valid file" } }

    val outputDir by option("-o", "--output-dir", help = "Directory to write one gzipped <sampleName>_readMapping.txt.gz file per VCF sample. Created if it does not exist. Required parameter.")
        .required()

    val useAlleleDepth by option("--use-allele-depth", help = "If set, weight each called allele's count by its AD (allele depth) FORMAT value instead of a flat 1. Falls back to 1 when AD is absent or malformed for a given genotype.")
        .flag(default = false)

    val maxBufferedEntries by option(
        "--max-buffered-entries",
        help = "Number of (sample, hapId set) counts to hold in memory before appending them to the output files " +
            "(at the next reference range boundary). Each entry costs roughly 10-20 bytes. Default: $DEFAULT_MAX_BUFFERED_ENTRIES"
    )
        .long()
        .default(DEFAULT_MAX_BUFFERED_ENTRIES)
        .validate { require(it > 0) { "--max-buffered-entries must be positive" } }

    companion object {
        const val DEFAULT_MAX_BUFFERED_ENTRIES = 100_000_000L

        /** Output files are gzipped; FindPaths accepts the "_readMapping.txt.gz" suffix. */
        const val READ_MAPPING_SUFFIX = "_readMapping.txt.gz"

        /** Full path of [sampleName]'s readMapping file in [outputDir] (no trailing slash). */
        fun readMappingFileName(outputDir: String, sampleName: String) = "$outputDir/$sampleName$READ_MAPPING_SUFFIX"
    }

    override fun run() {
        logCommand(this)

        //trim a trailing slash so "$outputDirPath/..." below never produces a double slash,
        //regardless of whether the user passed one
        val outputDirPath = outputDir.trimEnd('/')
        File(outputDirPath).mkdirs()

        val (sampleNames, stats) = buildReadMappings(indexFile, imputeVcf, outputDirPath, useAlleleDepth, maxBufferedEntries)
        writePathKeyFile(outputDirPath, sampleNames)

        logSummary(stats, sampleNames.size)
    }

    /**
     * Returns the count to add for one called [allele]. Flat 1 unless [useAlleleDepth] is set AND
     * the genotype has a usable AD value for this exact allele; falls back to 1 (counted in
     * [ReadMappingStats.adFallbackCount]) if AD is absent or the allele's index falls outside the
     * AD array. An AD value of 0 returns 0 (counted in [ReadMappingStats.zeroDepthAlleleCalls])
     * rather than falling back to 1 - zero reported depth is real information, not missing data.
     */
    fun alleleWeight(
        genotype: Genotype,
        allele: Allele,
        record: VariantContext,
        useAlleleDepth: Boolean,
        stats: ReadMappingStats
    ): Int {
        if (!useAlleleDepth) return 1
        if (!genotype.hasAD()) {
            stats.adFallbackCount++
            return 1
        }
        val ad = genotype.ad
        val alleleIndex = record.getAlleleIndex(allele)
        if (alleleIndex < 0 || alleleIndex >= ad.size) {
            stats.adFallbackCount++
            return 1
        }
        val depth = ad[alleleIndex]
        if (depth <= 0) {
            stats.zeroDepthAlleleCalls++
            return 0
        }
        return depth
    }

    /**
     * Converts one impute-VCF [record], already matched to its [indexPosition], into evidence.
     * [onEvidence] is called once per counted allele with the genotype's sample name, the hapId
     * set the allele resolves to, and its weight. (Keyed by name, not position: htsjdk does not
     * guarantee record.genotypes iterates in header order.) [stats] is mutated in place with a counter for every edge case encountered.
     */
    fun processRecord(
        record: VariantContext,
        indexPosition: VcfHaplotypeIndexPosition,
        useAlleleDepth: Boolean,
        stats: ReadMappingStats,
        onEvidence: (sampleName: String, hapIds: List<String>, weight: Int) -> Unit
    ) {
        record.genotypes.forEach { genotype: Genotype ->
            genotype.alleles.forEach { allele: Allele ->
                if (allele.isNoCall) {
                    stats.noCallAlleles++
                    return@forEach
                }

                val hapIds = indexPosition.alleleToHapIds[VcfHaplotypeIndexUtils.alleleKey(allele)]
                if (hapIds == null) {
                    stats.allelesNotInIndex++
                    return@forEach
                }

                val weight = alleleWeight(genotype, allele, record, useAlleleDepth, stats)
                if (weight <= 0) return@forEach

                stats.gameteHits++
                stats.totalEvidenceUnits += weight
                onEvidence(genotype.sampleName, hapIds, weight)
            }
        }
    }

    /**
     * Streams [indexFile] and [imputeVcfFile] together in one sorted merge pass and writes one
     * <sampleName>_readMapping.txt.gz per VCF sample into [outputDir] (no trailing slash). Returns
     * the VCF's sample names (header order) and the run's stats.
     *
     * Fails before counting anything if the two files list their shared contigs in different
     * orders, or if either file's records for a contig are not contiguous; fails mid-run if
     * positions within a contig go backwards in either file.
     */
    fun buildReadMappings(
        indexFile: String,
        imputeVcfFile: String,
        outputDir: String,
        useAlleleDepth: Boolean,
        maxBufferedEntries: Long = DEFAULT_MAX_BUFFERED_ENTRIES
    ): Pair<List<String>, ReadMappingStats> {

        val indexContigOrder = verifyContigOrder(indexFile, imputeVcfFile)
        val indexContigRank = HashMap<String, Int>()
        indexContigOrder.forEachIndexed { rank, contig -> indexContigRank[contig] = rank }

        val stats = ReadMappingStats()
        val sampleNames: List<String>

        VCFFileReader(File(imputeVcfFile), false).use { vcfReader ->
            sampleNames = vcfReader.header.sampleNamesInOrder.toList()
            val writer = ReadMappingFileWriter(outputDir, sampleNames, imputeVcfFile)
            val accumulator = ReadMappingAccumulator(sampleNames.size)
            val sampleIndexByName = HashMap<String, Int>()
            sampleNames.forEachIndexed { i, name -> sampleIndexByName[name] = i }
            val addEvidence = { sampleName: String, hapIds: List<String>, weight: Int ->
                accumulator.add(sampleIndexByName.getValue(sampleName), hapIds, weight)
            }

            VcfHaplotypeIndexStreamReader(indexFile).use { indexReader ->
                var cursor = indexReader.next()
                var currentRange: ReferenceRange? = null
                var lastVcfContig: String? = null
                var lastVcfPosition = Int.MIN_VALUE

                vcfReader.forEach { record ->
                    stats.totalImputeRecords++

                    val contig = record.contig
                    val start = record.start
                    if (contig == lastVcfContig) {
                        require(start >= lastVcfPosition) {
                            "$imputeVcfFile is not sorted by position: $contig:$start appears after $contig:$lastVcfPosition"
                        }
                    } else {
                        lastVcfContig = contig
                    }
                    lastVcfPosition = start

                    val recordRank = indexContigRank[contig]
                    if (recordRank == null) {
                        stats.positionsNotInIndex++
                        return@forEach
                    }

                    //advance the index past everything strictly before this record
                    while (cursor != null) {
                        val cursorRank = indexContigRank.getValue(cursor.position.contig)
                        if (cursorRank > recordRank || (cursorRank == recordRank && cursor.position.position >= start)) break
                        cursor = indexReader.next()
                    }

                    val indexPosition = cursor
                    if (indexPosition == null || indexPosition.position.position != start || indexPosition.position.contig != contig) {
                        stats.positionsNotInIndex++
                        return@forEach
                    }
                    stats.positionsMatched++

                    //Only flush at a reference range boundary: a hapId set never spans ranges, so
                    //this guarantees no set is ever written twice to the same sample's file.
                    if (indexPosition.refRange != currentRange) {
                        if (accumulator.bufferedEntries >= maxBufferedEntries) {
                            accumulator.flush(writer)
                            stats.flushes++
                        }
                        currentRange = indexPosition.refRange
                    }

                    processRecord(record, indexPosition, useAlleleDepth, stats, addEvidence)
                }
            }

            accumulator.flush(writer)
            stats.flushes++

            accumulator.sampleHasEvidence.forEachIndexed { i, hasEvidence ->
                if (!hasEvidence) stats.samplesWithNoEvidence.add(sampleNames[i])
            }
        }

        return sampleNames to stats
    }

    /**
     * One fast first-column pass over each file to read its contig order, then a check that the
     * shared contigs are in the same order in both. Returns the index's contig order.
     */
    private fun verifyContigOrder(indexFile: String, imputeVcfFile: String): List<String> {
        val startTime = System.currentTimeMillis()
        val indexContigOrder = VcfHaplotypeIndexUtils.scanContigOrder(indexFile)
        val vcfContigOrder = VcfHaplotypeIndexUtils.scanContigOrder(imputeVcfFile)
        myLogger.info(
            "Contig order check: index has ${indexContigOrder.size} contig(s), impute VCF has " +
                "${vcfContigOrder.size} (scanned in ${(System.currentTimeMillis() - startTime) / 1000}s)"
        )
        VcfHaplotypeIndexUtils.requireMatchingContigOrder(indexContigOrder, vcfContigOrder, indexFile, imputeVcfFile)

        val shared = indexContigOrder.toHashSet().apply { retainAll(vcfContigOrder.toSet()) }
        if (shared.size < vcfContigOrder.size) {
            myLogger.warn(
                "${vcfContigOrder.size - shared.size} impute VCF contig(s) are not in the index and will be skipped: " +
                    vcfContigOrder.filter { it !in shared }.take(20)
            )
        }
        return indexContigOrder
    }

    /**
     * Writes a pathKeyFile.txt (sampleName\tfilename per row) to [outputDir], so the directory
     * this command wrote to can be passed directly as `find-paths --path-keyfile`. [outputDir]
     * must not have a trailing slash.
     */
    fun writePathKeyFile(outputDir: String, sampleNames: Collection<String>) {
        getBufferedWriter("$outputDir/pathKeyFile.txt").use { writer ->
            writer.write("sampleName\tfilename\n")
            sampleNames.forEach { writer.write("$it\t${readMappingFileName(outputDir, it)}\n") }
        }
    }

    fun logSummary(stats: ReadMappingStats, numSamples: Int) {
        myLogger.info("build-read-mapping-from-vcf summary:")
        myLogger.info("  impute VCF records read: ${stats.totalImputeRecords}")
        myLogger.info("  positions matched in index: ${stats.positionsMatched}")
        myLogger.info("  positions not in index: ${stats.positionsNotInIndex}")
        myLogger.info("  no-call alleles skipped: ${stats.noCallAlleles}")
        myLogger.info("  alleles not found in index: ${stats.allelesNotInIndex}")
        myLogger.info("  gamete hits (evidence events): ${stats.gameteHits}")
        myLogger.info("  total evidence units written: ${stats.totalEvidenceUnits}")
        myLogger.info("  AD fallback to flat count: ${stats.adFallbackCount}")
        myLogger.info("  zero-depth allele calls skipped: ${stats.zeroDepthAlleleCalls}")
        myLogger.info("  buffer flushes to disk: ${stats.flushes}")
        myLogger.info("  samples with readMapping files written: $numSamples")

        if (stats.samplesWithNoEvidence.isNotEmpty()) {
            myLogger.warn(
                "${stats.samplesWithNoEvidence.size} sample(s) had zero resolvable evidence: " +
                    stats.samplesWithNoEvidence.take(50)
            )
        }

        if (stats.totalImputeRecords > 0 && stats.positionsMatched == 0L) {
            myLogger.error(
                "No impute-VCF position matched the index. This usually means the index and the " +
                    "impute VCF use different contig names."
            )
        }
    }
}

/**
 * Buffers per-sample counts keyed by a compact int id per distinct hapId set, instead of a
 * Map<List<String>, Int> per sample. Each hapId set is stored once no matter how many samples
 * hit it, and each (sample, set) count costs a couple of ints in a primitive hash map.
 */
class ReadMappingAccumulator(numSamples: Int) {

    //content-based interning: one id per distinct hapId set
    private val setIdByContent = Object2IntOpenHashMap<List<String>>().apply { defaultReturnValue(-1) }

    //identity cache in front of it: every sample at a given position passes the SAME List
    //instance, so this avoids re-hashing a (potentially long) hapId list for every sample
    private val setIdByIdentity = Reference2IntOpenHashMap<List<String>>().apply { defaultReturnValue(-1) }

    private val hapIdSets = ArrayList<List<String>>()
    private val countsBySample = Array(numSamples) { Int2IntOpenHashMap() }

    val sampleHasEvidence = BooleanArray(numSamples)

    /** Number of (sample, hapId set) counts currently held. */
    var bufferedEntries = 0L
        private set

    fun add(sampleIndex: Int, hapIds: List<String>, weight: Int) {
        var setId = setIdByIdentity.getInt(hapIds)
        if (setId < 0) {
            setId = setIdByContent.getInt(hapIds)
            if (setId < 0) {
                setId = hapIdSets.size
                hapIdSets.add(hapIds)
                setIdByContent.put(hapIds, setId)
            }
            setIdByIdentity.put(hapIds, setId)
        }
        //addTo returns the previous value (0 if absent); weights are always > 0
        if (countsBySample[sampleIndex].addTo(setId, weight) == 0) bufferedEntries++
        sampleHasEvidence[sampleIndex] = true
    }

    /**
     * Appends every buffered count to its sample's file and clears the buffer. Must only be called
     * at a reference range boundary (see [BuildReadMappingFromVcf.buildReadMappings]).
     */
    fun flush(writer: ReadMappingFileWriter) {
        val joined = arrayOfNulls<String>(hapIdSets.size)
        countsBySample.forEachIndexed { sampleIndex, counts ->
            if (counts.isEmpty()) return@forEachIndexed
            writer.append(sampleIndex) { out ->
                val iterator = counts.int2IntEntrySet().fastIterator()
                while (iterator.hasNext()) {
                    val entry = iterator.next()
                    val setId = entry.intKey
                    val hapIdString = joined[setId] ?: hapIdSets[setId].joinToString(",").also { joined[setId] = it }
                    out.write(hapIdString)
                    out.write("\t")
                    out.write(entry.intValue.toString())
                    out.write("\n")
                }
            }
            counts.clear()
        }
        setIdByContent.clear()
        setIdByIdentity.clear()
        hapIdSets.clear()
        bufferedEntries = 0
    }
}

/**
 * Owns the per-sample readMapping files: creates each one with its header up front (so every
 * sample gets a file even with zero evidence), then lets [ReadMappingAccumulator.flush] append
 * count rows. Files are opened only for the duration of a flush, so thousands of samples never
 * need thousands of simultaneously open file handles.
 *
 * Files are gzipped (the ".gz" suffix makes getBufferedWriter compress). Each append adds a new
 * gzip member to the file; concatenated members are a valid gzip stream per RFC 1952, and
 * getBufferedReader (java.util.zip.GZIPInputStream) reads all of them, so the file reads back as
 * one continuous text file.
 *
 * The header matches [AlignmentUtils.exportReadMapping], so [AlignmentUtils.importReadMapping]
 * and FindPaths read these files unchanged.
 */
class ReadMappingFileWriter(outputDir: String, sampleNames: List<String>, sourceFile: String) {

    private val files = sampleNames.map { File(BuildReadMappingFromVcf.readMappingFileName(outputDir, it)) }

    init {
        File(outputDir).mkdirs()
        sampleNames.forEachIndexed { i, sampleName ->
            getBufferedWriter(files[i]).use { out ->
                out.write("#sampleName=$sampleName\n")
                out.write("#filename1=$sourceFile\n")
                out.write("HapIds\tcount\n")
            }
        }
    }

    fun append(sampleIndex: Int, block: (java.io.BufferedWriter) -> Unit) {
        getBufferedWriter(files[sampleIndex], append = true).use(block)
    }
}
