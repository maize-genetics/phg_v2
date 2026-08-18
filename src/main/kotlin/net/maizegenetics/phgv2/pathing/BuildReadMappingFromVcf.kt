package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.getBufferedWriter
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * Running tallies for a single build-read-mapping-from-vcf run.
 */
data class ReadMappingStats(
    var totalImputeRecords: Int = 0,
    var positionsMatched: Int = 0,
    var positionsNotInIndex: Int = 0,
    var noCallAlleles: Int = 0,
    var allelesNotInIndex: Int = 0,
    var gameteHits: Int = 0,
    var totalEvidenceUnits: Long = 0,
    var adFallbackCount: Int = 0,
    var zeroDepthAlleleCalls: Int = 0
)

/**
 * Builds one readMapping file per sample from a VCF, using the index produced by
 * [BuildVcfHaplotypeIndex] to translate each of the VCF's observed alleles into the set of PHG
 * haplotype ids that carry it.
 *
 * This is step 2 of the VCF-to-readMapping pipeline. Step 1 ([BuildVcfHaplotypeIndex]) does the
 * expensive HaplotypeGraph work once, up front; this step needs only the small flat index file it
 * produces - no HaplotypeGraph is built or needed here.
 *
 * For every genotype call in --impute-vcf, each called allele (one per gamete - e.g. 2 for a
 * diploid site) is treated as one independent piece of evidence, exactly like a single sequencing
 * read: it contributes a +1 (or, with --use-allele-depth, its AD value) to the count for whichever
 * set of haplotypes the index says carries that allele at that position. Both alleles of a
 * heterozygous call land in the same per-sample readMapping map, since - like a real read - a
 * called allele doesn't reveal which parental chromosome it came from.
 *
 * Output files are named <sampleName>_readMapping.txt (the suffix FindPaths requires), one per
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

    val imputeVcf by option(help = "Full path to the VCF file to impute readMapping evidence from. Required parameter.")
        .required()
        .validate { require(File(it).isFile) { "$it is not a valid file" } }

    val outputDir by option("-o", "--output-dir", help = "Directory to write one <sampleName>_readMapping.txt file per VCF sample. Created if it does not exist. Required parameter.")
        .required()

    val useAlleleDepth by option("--use-allele-depth", help = "If set, weight each called allele's count by its AD (allele depth) FORMAT value instead of a flat 1. Falls back to 1 when AD is absent or malformed for a given genotype.")
        .flag(default = false)

    override fun run() {
        logCommand(this)

        //trim a trailing slash so "$outputDirPath/..." below never produces a double slash,
        //regardless of whether the user passed one
        val outputDirPath = outputDir.trimEnd('/')
        File(outputDirPath).mkdirs()

        val (perSampleReadMappings, stats) = buildReadMappings(indexFile, imputeVcf, useAlleleDepth)

        perSampleReadMappings.forEach { (sampleName, hapIdMapping) ->
            val outputFile = "$outputDirPath/${sampleName}_readMapping.txt"
            AlignmentUtils.exportReadMapping(outputFile, hapIdMapping, sampleName, Pair(imputeVcf, ""))
        }
        writePathKeyFile(outputDirPath, perSampleReadMappings.keys)

        logSummary(stats, perSampleReadMappings)
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
     * Converts one impute-VCF [record] into per-sample count deltas. Returns an empty map when
     * the record's position isn't covered by [index]. [stats] is mutated in place with a counter
     * for every edge case encountered.
     */
    fun processRecord(
        record: VariantContext,
        index: Map<Position, VcfHaplotypeIndexPosition>,
        useAlleleDepth: Boolean,
        stats: ReadMappingStats
    ): Map<String, Map<List<String>, Int>> {

        stats.totalImputeRecords++

        val indexPosition = index[Position(record.contig, record.start)]
        if (indexPosition == null) {
            stats.positionsNotInIndex++
            return emptyMap()
        }
        stats.positionsMatched++

        val deltas = mutableMapOf<String, MutableMap<List<String>, Int>>()

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
                val sampleMap = deltas.getOrPut(genotype.sampleName) { mutableMapOf() }
                sampleMap[hapIds] = (sampleMap[hapIds] ?: 0) + weight
            }
        }

        return deltas
    }

    /**
     * Loads [indexFile], streams [imputeVcfFile] once, and accumulates a
     * Map<List<hapId>, count> per sample. Memory is bounded by the number of distinct hapId sets
     * observed across all samples, not by the number of VCF records - the same bound the existing
     * kmer/ropebwt readMapping builders already accept.
     *
     * Every sample named in the VCF header gets an entry in the returned map, even if it never
     * contributes any evidence, so a downstream keyfile-driven step can assume every sample's
     * readMapping file exists.
     */
    fun buildReadMappings(
        indexFile: String,
        imputeVcfFile: String,
        useAlleleDepth: Boolean
    ): Pair<Map<String, Map<List<String>, Int>>, ReadMappingStats> {

        val index = VcfHaplotypeIndexUtils.readIndex(indexFile)
        val stats = ReadMappingStats()
        val perSample = mutableMapOf<String, MutableMap<List<String>, Int>>()

        val allSampleNames: List<String>
        VCFFileReader(File(imputeVcfFile), false).use { reader ->
            allSampleNames = reader.header.sampleNamesInOrder
            reader.forEach { record ->
                processRecord(record, index, useAlleleDepth, stats).forEach { (sampleName, hapIdCounts) ->
                    val sampleMap = perSample.getOrPut(sampleName) { mutableMapOf() }
                    hapIdCounts.forEach { (hapIds, count) -> sampleMap[hapIds] = (sampleMap[hapIds] ?: 0) + count }
                }
            }
        }

        //every VCF sample gets a file, even with zero evidence
        allSampleNames.forEach { perSample.getOrPut(it) { mutableMapOf() } }

        return perSample to stats
    }

    /**
     * Writes a pathKeyFile.txt (sampleName\tfilename per row) to [outputDir], so the directory
     * this command wrote to can be passed directly as `find-paths --path-keyfile`. [outputDir]
     * must not have a trailing slash.
     */
    fun writePathKeyFile(outputDir: String, sampleNames: Collection<String>) {
        getBufferedWriter("$outputDir/pathKeyFile.txt").use { writer ->
            writer.write("sampleName\tfilename\n")
            sampleNames.forEach { writer.write("$it\t$outputDir/${it}_readMapping.txt\n") }
        }
    }

    fun logSummary(stats: ReadMappingStats, perSampleReadMappings: Map<String, Map<List<String>, Int>>) {
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
        myLogger.info("  samples with readMapping files written: ${perSampleReadMappings.size}")

        val samplesWithNoEvidence = perSampleReadMappings.filterValues { it.isEmpty() }.keys
        if (samplesWithNoEvidence.isNotEmpty()) {
            myLogger.warn("${samplesWithNoEvidence.size} sample(s) had zero resolvable evidence: $samplesWithNoEvidence")
        }

        if (stats.totalImputeRecords > 0 && stats.positionsMatched == 0) {
            myLogger.error(
                "No impute-VCF position matched the index. This usually means the index and the " +
                    "impute VCF use different contig names."
            )
        }
    }
}
