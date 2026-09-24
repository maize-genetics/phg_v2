package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.UsageError
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.choice
import com.github.ajalt.clikt.parameters.types.double
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.getBufferedWriter
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * Infers, for every sample in a VCF, which reference-panel founders it descends from along each contig.
 *
 * The observations are **genotypes**, not reads. A candidate founder pair implies a genotype by Mendelian
 * sampling, the sample's genotype is known, and [VcfGenotypeEmissionProbability] scores one against the
 * other; the Viterbi recursion then finds the founder path that best explains the whole contig. Nothing
 * reconstructs reads, which is the point -- a ps4g intermediate would have to invent a read count for a
 * call made at an unknown depth.
 *
 * Pair the resulting path with a **denser** panel through `bed-to-vcf` and the output carries that
 * panel's sites: that is how a sparse assay becomes a dense genotype set.
 *
 * ## How it differs from `impute-path-from-ps4g`
 *
 * `--prob-correct` and the choice of haploid or diploid path mean the same as they do there, and haploid
 * is the same mechanism: a diploid path at an inbreeding coefficient of one, where only the homozygous
 * states are reachable.
 *
 * Two things are different.
 *
 * **The transition is charged per base, not per observation.** There is no bin grid here -- the HMM's
 * positions are the VCF's own sites, at whatever spacing the assay happens to have -- so a fixed
 * per-step penalty would make the effective recombination rate depend on marker density. `--prob-switch`
 * is therefore a rate: the probability of a switch across `--prob-switch-distance`, scaled to each
 * step's actual distance by [transitionLogsByDistance]. It is also expressed as the probability of
 * switching rather than of staying, which is the way round users think about recombination.
 *
 * **There is no presence/absence correction.** That correction exists because a read cannot support a
 * founder that has no sequence at a locus. A genotype carries no such asymmetry: if the panel has no
 * call for a founder the emission simply treats that founder's pairs as uninformative there.
 *
 * ## What the path covers
 *
 * A contig's path spans the first to the last site shared with the panel, and claims nothing beyond
 * either. `--extend-to-contig-ends` carries the terminal intervals to the contig's ends instead. See
 * [applyContigBounds].
 *
 * ## Assumptions
 *
 * The reference panel is haploid, or diploid and **mostly homozygous** -- each founder is read as
 * contributing gametes from its own called alleles, so a heterozygous panel site is handled correctly
 * but a panel that is heterozygous everywhere is not what this models. The sample VCF is **unphased**;
 * genotypes are compared as unordered pairs.
 */
class ImputePathFromVcf : CliktCommand(help = "Impute founder paths for the samples in a VCF, using a reference panel VCF") {

    private val myLogger = LogManager.getLogger(ImputePathFromVcf::class.java)

    val toImputeVcf by option(help = "VCF holding the samples to impute. Coordinate sorted, unphased, " +
            "and using the same contig names as the panel. Required parameter.")
        .required()
        .validate { require(File(it).exists()) { "$it is not a valid file" } }

    val panelVcf by option(help = "Reference panel VCF whose samples are the candidate founders. " +
            "Coordinate sorted, and assumed haploid or mostly homozygous. Its ##contig header lines " +
            "define the contig order both files must follow. Required parameter.")
        .required()
        .validate { require(File(it).exists()) { "$it is not a valid file" } }

    val outPathDir by option(help = "Directory for the imputed paths, written as " +
            "<sampleName>_imputed_path.bed. Coordinates are 0-based half-open. Required parameter.")
        .required()

    val pathType by option(help = "The type of path to find. 'haploid' infers a single founder per " +
            "position, which is a diploid path with the inbreeding coefficient set to 1. 'diploid' " +
            "infers a pair.")
        .choice("haploid", "diploid")
        .default("haploid")

    val probCorrect by option(help = "The probability that a genotype call is correct. A mismatch " +
            "between the observed genotype and the one a founder pair predicts costs ln(1 - this), " +
            "about -3.9 at the default. Default = 0.98")
        .double()
        .default(0.98)
        .validate { require(it > 0.5 && it < 1.0) { "prob-correct must be between 0.5 and 1.0 exclusive" } }

    val probSwitch by option(help = "The probability of a path switch (a recombination) across " +
            "--prob-switch-distance, scaled to each step's actual distance. Conceptually a " +
            "recombination rate, but in practice a smoothing parameter: the read-based pipeline's " +
            "tuned equivalent is about 1300 times less recombination than real maize, because it has " +
            "to resist noisy evidence. At the default about 3 to 4 consecutive mismatching sites " +
            "justify a switch. Default = 1e-4")
        .double()
        .default(1e-4)
        .validate { require(it > 0.0 && it < 1.0) { "prob-switch must be between 0 and 1 exclusive" } }

    val probSwitchDistance by option(help = "The distance --prob-switch is quoted over, in base " +
            "pairs. Default = 1000000")
        .double()
        .default(1_000_000.0)
        .validate { require(it > 0.0) { "prob-switch-distance must be positive" } }

    val inbreedCoef by option(help = "The inbreeding coefficient, used for diploid paths. Only 0.0 " +
            "and 1.0 are supported: an intermediate value needs the general Viterbi scan, whose " +
            "transition matrix would have to be rebuilt at every position once the transition varies " +
            "with distance. --path-type haploid sets this to 1.0 regardless. Default = 0.0")
        .double()
        .default(0.0)
        .validate {
            require(it == 0.0 || it == 1.0) {
                "inbreed-coef must be 0.0 or 1.0 for this command; $it would need the general scan, " +
                        "which cannot take a distance-scaled transition"
            }
        }

    val extendToContigEnds by option(help = "Carry the first and last interval of each contig out to " +
            "the contig's ends. Off by default, so a path spans only the first to the last site shared " +
            "with the panel and a denser panel's sites outside that span get no call. There is no " +
            "evidence of ancestry beyond the terminal markers; set this to assume it continues.")
        .flag()

    val contigsToUse by option(help = "Comma-separated contigs to impute, or a file with one per " +
            "line. All contigs shared by the two VCFs if omitted.")
        .default("")

    override fun run() {
        logCommand(this)
        File(outPathDir).mkdirs()
        val result = imputeFounderPaths(File(toImputeVcf), File(panelVcf), parameters())
        logPairingCounts(result.counts, myLogger)
        writeFounderPathBeds(result.paths, File(outPathDir), pathType == "haploid", myLogger)
    }

    /** This command's options as the parameter bundle the shared path finder takes. */
    fun parameters() = VcfPathParameters(
        pathType = pathType,
        inbreedCoef = inbreedCoef,
        probCorrect = probCorrect,
        probSwitch = probSwitch,
        probSwitchDistance = probSwitchDistance,
        contigsToUse = buildContigSet(contigsToUse),
        extendToContigEnds = extendToContigEnds
    )

    /** Comma-separated list, or a file with one contig per line, or empty for all. */
    fun buildContigSet(value: String): Set<String> {
        if (value.isBlank()) return emptySet()
        val asFile = File(value)
        val contigs = if (asFile.exists()) {
            asFile.readLines().map { it.trim() }.filter { it.isNotEmpty() }
        } else {
            value.split(",").map { it.trim() }.filter { it.isNotEmpty() }
        }
        // A UsageError rather than a require, so a mistake in the invocation exits with a readable
        // message instead of a stack trace, as ImputePathFromPs4g does for the same option.
        if (contigs.isEmpty()) throw UsageError(
            "--contigs-to-use was given as '$value' but yielded no contigs. An empty set would mean " +
                    "something different from omitting the option, so it is reported instead."
        )
        return contigs.toSet()
    }
}

/**
 * Everything the VCF path finder needs beyond the two files, so the chained VCF-to-VCF command can
 * drive it without reconstructing a Clikt command.
 */
data class VcfPathParameters(
    val pathType: String = "haploid",
    val inbreedCoef: Double = 0.0,
    val probCorrect: Double = 0.98,
    val probSwitch: Double = 1e-4,
    val probSwitchDistance: Double = 1_000_000.0,
    val contigsToUse: Set<String> = emptySet(),
    val extendToContigEnds: Boolean = false
) {
    /** The coefficient actually used: a haploid path is a diploid path at 1. */
    val coefficient: Double get() = if (pathType == "haploid") 1.0 else inbreedCoef
}

/** A founder path per sample, with what the site pairing saw while producing it. */
data class VcfPathResult(
    val paths: Map<String, List<PathInterval<Pair<String, String>>>>,
    val counts: PairingCounts
)

/**
 * Infers a founder path for every sample in [toImputeVcf] against the founders in [panelVcf].
 *
 * Returns the intervals rather than writing them, so the chained command can compose them straight
 * into a VCF with no BED ever reaching disk. Each sample is its own Viterbi run over the same
 * [ContigSites]: the panel's founder alleles are shared, each sample's genotypes are not.
 */
fun imputeFounderPaths(
    toImputeVcf: File,
    panelVcf: File,
    parameters: VcfPathParameters,
    logger: org.apache.logging.log4j.Logger = LogManager.getLogger("ImputeFounderPaths")
): VcfPathResult {
    val hmm = ViterbiHMM(parameters.coefficient, 1.0 - parameters.probSwitch, parameters.probCorrect)
    val intervals = mutableMapOf<String, MutableList<PathInterval<Pair<String, String>>>>()

    // Contig lengths, where the panel declares them, needed only by --extend-to-contig-ends. They are
    // optional in the VCF specification and real panels omit them.
    val contigLengths: Map<String, Int> = VCFFileReader(panelVcf, false).use { reader ->
        reader.fileHeader.contigLines
            .filter { it.genericFields["ID"] != null && it.genericFields["length"] != null }
            .associate { it.genericFields["ID"]!! to it.genericFields["length"]!!.toInt() }
    }

    val counts = pairVcfSites(toImputeVcf, panelVcf, parameters.contigsToUse) { sites ->
        val distances = stepDistances(sites.positions)
        // A haploid recursion switches among founders; the diploid fast path switches one haplotype at
        // a time, so its alternatives are counted in founders too. The divisor is the same.
        val (lnNoSwitch, lnSwitch) = transitionLogsByDistance(
            distances, parameters.probSwitch, parameters.probSwitchDistance, sites.nFounders)

        for (sampleIndex in sites.sampleNames.indices) {
            val emission = VcfGenotypeEmissionProbability(sites, sampleIndex, parameters.probCorrect)
            val path = hmm.findDiploidStatePath(
                sites.nFounders, sites.nSites,
                emission::getDiploidEmissionProbabilityArray,
                null, lnNoSwitch, lnSwitch
            )
            val calls = path.first.mapIndexed { index, state ->
                Pair(
                    Position(sites.contig, sites.positions[index]),
                    Pair(
                        sites.founderNames[state / sites.nFounders],
                        sites.founderNames[state % sites.nFounders]
                    )
                )
            }
            // binSize of 1 because the positions are base pairs, not bin indices: the midpoint cut
            // rule then places each boundary halfway between two real sites.
            intervals.getOrPut(sites.sampleNames[sampleIndex]) { mutableListOf() }
                .addAll(applyContigBounds(
                    pathToIntervals(calls, 1), sites.positions.first(),
                    contigLengths[sites.contig], parameters.extendToContigEnds))
        }
        logger.info("${sites.contig}: ${sites.nSites} sites, ${sites.nSamples} samples, " +
                "${sites.nFounders} founders")
    }
    return VcfPathResult(intervals, counts)
}

/**
 * Logs what the pairing saw. The absent-from-panel count is the one to read: those sites are silently
 * dropped, so a sample and panel that do not correspond -- a different reference build, different
 * contig naming, a panel that does not cover the assay -- would otherwise yield a confident path drawn
 * from a fraction of the evidence.
 */
fun logPairingCounts(counts: PairingCounts, logger: org.apache.logging.log4j.Logger) {
    val offered = counts.sitesUsed + counts.sampleSitesNotInPanel
    val percentUsed = if (offered > 0) 100.0 * counts.sitesUsed / offered else 0.0
    logger.info("Sites used from the sample VCF: ${counts.sitesUsed}")
    logger.info("Sites in the sample VCF not present in the panel: ${counts.sampleSitesNotInPanel} " +
            "(${"%.1f".format(100.0 - percentUsed)}% of $offered offered)")
    logger.info("Sites in the panel not present in the sample VCF: ${counts.panelSitesNotInSample}")
    if (counts.duplicateKeysSkipped > 0) {
        logger.info("Records skipped for sharing a contig, position and REF with one already used: " +
                "${counts.duplicateKeysSkipped}")
    }
    if (counts.sitesUsed == 0L) {
        logger.warn("No sites were shared by the two VCFs. Check that they use the same contig names " +
                "and the same reference coordinates.")
    } else if (percentUsed < 50.0) {
        logger.warn("Only ${"%.1f".format(percentUsed)}% of the sample VCF's sites are in the panel. " +
                "The path rests on that fraction alone.")
    }
}

/** One BED per sample, named as `impute-path-from-ps4g` names its own so `bed-to-vcf` reads either. */
fun writeFounderPathBeds(
    paths: Map<String, List<PathInterval<Pair<String, String>>>>,
    outputDir: File,
    haploid: Boolean,
    logger: org.apache.logging.log4j.Logger
) {
    outputDir.mkdirs()
    for ((sample, sampleIntervals) in paths) {
        val file = File(outputDir, "${sample}_imputed_path.bed")
        getBufferedWriter(file).use { writer ->
            writer.write(if (haploid) "chrom\tstart\tend\tparent1\n"
                         else "chrom\tstart\tend\tparent1\tparent2\n")
            for (interval in sampleIntervals) {
                val call = if (haploid) {
                    check(interval.call.first == interval.call.second) {
                        "haploid path produced a heterozygous call ${interval.call}"
                    }
                    interval.call.first
                } else "${interval.call.first}\t${interval.call.second}"
                writer.write("${interval.contig}\t${interval.start}\t${interval.end}\t$call\n")
            }
        }
        logger.info("Wrote $file: ${sampleIntervals.size} intervals")
    }
}

/**
 * Trims a contig's path to the sites it was actually inferred from, or, when [extendToContigEnds] is
 * set, carries it out to both ends of the contig.
 *
 * [pathToIntervals] is asymmetric: its first interval starts at 0 rather than at the first observation,
 * while its last ends at the last observation. On a uniform bin grid that is invisible, because the
 * first and last occupied bins sit near the ends of the contig anyway. On VCF sites it is not -- a
 * marker set covering the middle 60% of a chromosome would have the leading 20% filled in and the
 * trailing 20% left uncalled, for no reason a user could infer.
 *
 * The default is to claim neither: the path starts at the first shared site and ends at the last, and a
 * denser panel's sites outside that span get no call. Beyond the terminal markers there is no evidence
 * of ancestry, and a telomere-proximal recombination is exactly where it would be missed.
 *
 * [extendToContigEnds] fills both ends instead, for a caller who would rather assume the terminal
 * ancestry continues -- reasonable when the markers nearly reach the ends, or when a complete call set
 * matters more than the risk at the edges. The trailing end uses [contigLength] where the panel declares
 * one; `##contig` lines are optional and real panels omit them, so where it is absent the interval is
 * carried to [Int.MAX_VALUE], which covers any position a denser panel could hold.
 *
 * @param firstPosition the first position the path was inferred at, 1-based as a VCF carries it.
 */
fun applyContigBounds(
    intervals: List<PathInterval<Pair<String, String>>>,
    firstPosition: Int,
    contigLength: Int?,
    extendToContigEnds: Boolean
): List<PathInterval<Pair<String, String>>> {
    if (intervals.isEmpty()) return intervals
    val bounded = intervals.toMutableList()

    // 0-based half-open, so a start of firstPosition - 1 makes firstPosition the first base covered.
    val start = if (extendToContigEnds) 0 else firstPosition - 1
    if (start != bounded.first().start) bounded[0] = bounded.first().copy(start = start)

    if (extendToContigEnds) {
        val end = contigLength ?: Int.MAX_VALUE
        if (end > bounded.last().end) bounded[bounded.lastIndex] = bounded.last().copy(end = end)
    }

    // Defensive: an interval that trimming empties would be zero-length, which has no 1-based
    // inclusive form and BedToVcf rejects. pathToIntervals cannot produce one here, since its first cut
    // is the midpoint of the first two sites and so never falls below the first site.
    return bounded.filter { it.end > it.start }
}
