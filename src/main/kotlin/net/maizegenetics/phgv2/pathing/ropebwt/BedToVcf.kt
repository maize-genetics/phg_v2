package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.UsageError
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import com.google.common.collect.TreeRangeMap
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.utils.createGenericHeader
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * A founder path for one sample: for each contig, the reference intervals over which the sample's
 * two haplotypes are assigned to a pair of founders.
 *
 * Keyed by contig with plain integer positions rather than by [net.maizegenetics.phgv2.utils.Position].
 * `Position.compareTo` ignores the position entirely when the contigs differ, and can report 0 for
 * two positions that are not equal -- `Position("chr1", 5)` against `Position("1", 900)` compares
 * equal, because it strips "chr" before comparing numerically. A [TreeRangeMap] depends on a total
 * order consistent with equality, so splitting by contig first keeps the comparison to integers
 * within one contig, where it is unambiguous.
 */
typealias FounderPath = Map<String, RangeMap<Int, Pair<String, String>>>

/**
 * Composes imputed founder paths into a VCF by reading the founder alleles out of a reference panel.
 *
 * The path says which founders a sample descends from over each interval; the panel says what allele
 * each founder carries at each site. Together they give the sample's genotype at every site in the
 * panel, which is what makes this the density-raising step: feed a path inferred from sparse data
 * and a dense panel, and the output carries the dense panel's sites.
 *
 * Ported from the grits `BedToVcf`, with three changes:
 *
 *  - the duplicated VCF header builders are dropped in favour of
 *    [net.maizegenetics.phgv2.utils.createGenericHeader]
 *  - paths are keyed by contig and integer position rather than by `Position` (see [FounderPath])
 *  - the recognised BED file names include PHG's own `<sample>_imputed_path.bed`, which the grits
 *    naming did not cover (see [sampleNameOf])
 *
 * ## The panel is read as haploid founders
 *
 * Each panel sample contributes **one** allele per site, its first. A founder panel is normally
 * inbred, so this is usually exact; where a founder is heterozygous at a site its second allele is
 * ignored and the first is used for both haplotypes that descend from it. That is a deliberate
 * simplification of the output step, not an oversight.
 *
 * ## Sites with no path
 *
 * A sample whose path does not cover a site gets no genotype in that record, which the writer emits
 * as missing. Paths from `impute-path-from-ps4g` are gapless within a contig, so in practice this
 * arises for contigs the path does not mention at all.
 *
 * ## INFO fields
 *
 * Records are built from the panel's own [VariantContext], so the panel's INFO fields carry through
 * unchanged. Some of them -- `AF`, `DP`, `NS` -- describe the panel rather than the imputed samples
 * and are misleading here. Preserved for now because the affinity-model comparison was scored
 * against output produced this way; changing it is a decision about output content, not a fix.
 */
class BedToVcf : CliktCommand(help = "Compose imputed founder paths (BED) into a VCF using a reference panel") {

    private val myLogger = LogManager.getLogger(BedToVcf::class.java)

    val bedDir by option(help = "Directory of founder-path BED files, one or more per sample. " +
            "Recognised names are <sample>_imputed_path.bed as written by impute-path-from-ps4g, " +
            "<sample>_chr<contig>_imputed.bed, and <sample>.bed. Several files for one sample are " +
            "merged. Required parameter.")
        .required()
        .validate { require(File(it).isDirectory) { "$it is not a directory" } }

    val referencePanelVcf by option(help = "Reference panel VCF supplying each founder's allele at " +
            "each site. Its sample names must match the founder names in the BED files, and its " +
            "sites determine the sites in the output. Only the first allele of each panel genotype " +
            "is used. Required parameter.")
        .required()
        .validate { require(File(it).exists()) { "$it is not a valid file" } }

    val outputFile by option(help = "The VCF to write. Required parameter.")
        .required()

    override fun run() {
        logCommand(this)
        val paths = readPaths(File(bedDir))
        // A UsageError rather than a require: an empty directory is a mistake in the invocation, so
        // it should exit with a readable message rather than an uncaught exception, matching how
        // ImputePathFromPs4g reports an empty --contigs-to-use file.
        if (paths.isEmpty()) throw UsageError(
            "No usable BED files found in $bedDir. Expected names are <sample>_imputed_path.bed, " +
                    "<sample>_chr<contig>_imputed.bed, or <sample>.bed."
        )
        writeVcf(File(referencePanelVcf), paths, File(outputFile))
    }

    /**
     * Reads every BED file in [dir], grouped into one [FounderPath] per sample. Several files for
     * one sample -- the per-contig form -- are merged into that sample's single path.
     */
    fun readPaths(dir: File): Map<String, FounderPath> {
        val paths = mutableMapOf<String, MutableMap<String, RangeMap<Int, Pair<String, String>>>>()
        val bedFiles = dir.listFiles { file: File -> file.isFile && file.extension == "bed" }
            ?.sortedBy { it.name } ?: emptyList()
        for (bedFile in bedFiles) {
            val sample = sampleNameOf(bedFile)
            val forSample = paths.getOrPut(sample) { mutableMapOf() }
            var records = 0
            readPath(bedFile) { contig, start, end, parent1, parent2 ->
                forSample.getOrPut(contig) { TreeRangeMap.create() }
                    .put(Range.closed(start, end), Pair(parent1, parent2))
                records++
            }
            myLogger.info("${bedFile.name}: $records intervals for sample $sample")
        }
        return paths
    }

    /**
     * Strips the known suffixes to recover the sample name.
     *
     * `impute-path-from-ps4g` writes one file per sample named `<sample>_imputed_path.bed`, while
     * the per-contig convention used elsewhere is `<sample>_chr<contig>_imputed.bed`. Both have to
     * be recognised, or a path written by this project's own imputer would be read back under a
     * sample name with the suffix still attached.
     */
    fun sampleNameOf(bedFile: File): String {
        val name = bedFile.nameWithoutExtension
        if (name.endsWith("_imputed_path")) return name.removeSuffix("_imputed_path")
        val perContig = Regex("""^(.+)_chr[A-Za-z0-9._-]+_imputed$""").matchEntire(name)
        if (perContig != null) return perContig.groupValues[1]
        return name
    }

    /**
     * Parses one BED file, calling [record] for each interval with **1-based inclusive** bounds.
     *
     * The file is 0-based half-open per the BED specification; the conversion to the 1-based
     * inclusive bounds a VCF position is compared against is `start + 1` and `end`. A four-column
     * file -- what a haploid path writes -- is read as homozygous, with the single parent used for
     * both haplotypes.
     */
    fun readPath(bedFile: File, record: (contig: String, start: Int, end: Int, parent1: String, parent2: String) -> Unit) {
        bedFile.forEachLine { line ->
            if (line.isBlank() || line.startsWith("chrom") || line.startsWith("#")) return@forEachLine
            val fields = line.split('\t')
            if (fields.size < 4) return@forEachLine
            val start = fields[1].trim().toInt()
            val end = fields[2].trim().toInt()
            // A zero-length interval covers no bases and has no 1-based inclusive form.
            if (end <= start) return@forEachLine
            val parent1 = fields[3].trim()
            val parent2 = if (fields.size >= 5) fields[4].trim() else parent1
            record(fields[0], start + 1, end, parent1, parent2)
        }
    }

    /**
     * The same [FounderPath] structure [readPaths] builds, but from paths held in memory rather than
     * read from BED files, so a command that has just inferred them can compose without writing any.
     *
     * Intervals are 0-based half-open, as a BED record is, and become 1-based inclusive here by the
     * same `(start + 1, end)` conversion [readPath] applies -- the two routes have to agree, or a path
     * would land differently depending on whether it went through a file.
     */
    fun pathsFromIntervals(
        paths: Map<String, List<PathInterval<Pair<String, String>>>>
    ): Map<String, FounderPath> = paths.mapValues { (_, intervals) ->
        val byContig = mutableMapOf<String, RangeMap<Int, Pair<String, String>>>()
        for (interval in intervals) {
            if (interval.end <= interval.start) continue
            byContig.getOrPut(interval.contig) { TreeRangeMap.create() }
                .put(Range.closed(interval.start + 1, interval.end), interval.call)
        }
        byContig
    }

    /**
     * Streams [panelVcf] and writes one output record per panel site, carrying a genotype for every
     * sample in [paths] whose path covers that site.
     *
     * The panel is read one record at a time and never held, so memory is bounded by the paths.
     */
    fun writeVcf(panelVcf: File, paths: Map<String, FounderPath>, outputFile: File) {
        val samples = paths.keys.sorted()
        var sites = 0L
        var genotypesWritten = 0L
        VariantContextWriterBuilder()
            .unsetOption(Options.INDEX_ON_THE_FLY)
            .setOutputFile(outputFile)
            .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
            .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
            .build().use { writer ->
                writer.writeHeader(createGenericHeader(samples, emptySet()))
                VCFFileReader(panelVcf, false).use { reader ->
                    for (variant in reader) {
                        val genotypes = composeGenotypes(variant, paths, samples)
                        sites++
                        genotypesWritten += genotypes.size
                        writer.add(VariantContextBuilder(variant).genotypes(genotypes).make())
                    }
                }
            }
        myLogger.info("Wrote $outputFile: $sites sites, ${samples.size} samples, " +
                "$genotypesWritten genotypes called")
    }

    /**
     * Builds the genotypes for one panel site: look the site up in each sample's path to get its
     * founder pair, then read each founder's allele out of the panel record.
     *
     * A founder the panel does not carry at this site yields a no-call for that haplotype rather
     * than dropping the sample, so a partially unknown genotype is still reported.
     */
    fun composeGenotypes(
        variant: VariantContext,
        paths: Map<String, FounderPath>,
        samples: List<String>
    ): List<Genotype> {
        val founderAlleles = founderToAllele(variant)
        val genotypes = ArrayList<Genotype>(samples.size)
        for (sample in samples) {
            val founders = paths[sample]?.get(variant.contig)?.get(variant.start) ?: continue
            val allele1 = founderAlleles[founders.first] ?: Allele.NO_CALL
            val allele2 = founderAlleles[founders.second] ?: Allele.NO_CALL
            genotypes.add(GenotypeBuilder(sample).alleles(listOf(allele1, allele2)).make())
        }
        return genotypes
    }

    /**
     * Each panel sample's allele at this site: its first called allele, per the haploid-founder
     * reading described on the class. Samples with no called allele are omitted, so a founder
     * missing here becomes a no-call rather than a wrong call.
     */
    fun founderToAllele(variant: VariantContext): Map<String, Allele> {
        val alleles = HashMap<String, Allele>(variant.nSamples * 2)
        for (genotype in variant.genotypes) {
            val first = genotype.alleles.firstOrNull() ?: continue
            if (first.isNoCall) continue
            alleles[genotype.sampleName] = first
        }
        return alleles
    }
}
