package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.UsageError
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.choice
import com.github.ajalt.clikt.parameters.types.double
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

    val contigsToUse by option(help = "Comma-separated contigs to impute, or a file with one per " +
            "line. All contigs shared by the two VCFs if omitted.")
        .default("")

    override fun run() {
        logCommand(this)
        File(outPathDir).mkdirs()
        val paths = imputePaths()
        writeBedFiles(paths, File(outPathDir))
    }

    /**
     * Infers a path per sample, returning the intervals rather than a file, so the chained VCF-to-VCF
     * command can compose them directly without a BED ever being written.
     *
     * Each sample is a separate Viterbi run over the same [ContigSites], because the panel's founder
     * alleles are shared but each sample's genotypes are its own.
     */
    fun imputePaths(): Map<String, List<PathInterval<Pair<String, String>>>> {
        val coefficient = if (pathType == "haploid") 1.0 else inbreedCoef
        val hmm = ViterbiHMM(coefficient, 1.0 - probSwitch, probCorrect)
        val intervals = mutableMapOf<String, MutableList<PathInterval<Pair<String, String>>>>()

        val counts = pairVcfSites(File(toImputeVcf), File(panelVcf), buildContigSet(contigsToUse)) { sites ->
            val distances = stepDistances(sites.positions)
            // A haploid recursion switches among founders; the diploid fast path switches one haplotype
            // at a time, so its alternatives are also counted in founders. The divisor is the same.
            val (lnNoSwitch, lnSwitch) =
                transitionLogsByDistance(distances, probSwitch, probSwitchDistance, sites.nFounders)

            for (sampleIndex in sites.sampleNames.indices) {
                val emission = VcfGenotypeEmissionProbability(sites, sampleIndex, probCorrect)
                val path = hmm.findDiploidStatePath(
                    sites.nFounders, sites.nSites,
                    emission::getDiploidEmissionProbabilityArray,
                    null,
                    lnNoSwitch, lnSwitch
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
                // binSize of 1 because the positions are already base pairs, not bin indices: the
                // midpoint cut rule then places each boundary halfway between two real sites.
                intervals.getOrPut(sites.sampleNames[sampleIndex]) { mutableListOf() }
                    .addAll(pathToIntervals(calls, 1))
            }
            myLogger.info("${sites.contig}: ${sites.nSites} sites, ${sites.nSamples} samples, " +
                    "${sites.nFounders} founders")
        }

        reportCounts(counts)
        return intervals
    }

    /**
     * Logs what the pairing saw. The absent-from-panel count is the one to read: those sites are
     * silently dropped, so a sample and panel that do not correspond -- a different reference build,
     * different contig naming, a panel that does not cover the assay -- would otherwise produce a
     * confident path from almost no evidence.
     */
    private fun reportCounts(counts: PairingCounts) {
        val offered = counts.sitesUsed + counts.sampleSitesNotInPanel
        val percentUsed = if (offered > 0) 100.0 * counts.sitesUsed / offered else 0.0
        myLogger.info("Sites used from the sample VCF: ${counts.sitesUsed}")
        myLogger.info("Sites in the sample VCF not present in the panel: " +
                "${counts.sampleSitesNotInPanel} (${"%.1f".format(100.0 - percentUsed)}% of " +
                "$offered offered)")
        myLogger.info("Sites in the panel not present in the sample VCF: ${counts.panelSitesNotInSample}")
        if (counts.duplicateKeysSkipped > 0) {
            myLogger.info("Records skipped for sharing a contig, position and REF with one already " +
                    "used: ${counts.duplicateKeysSkipped}")
        }
        if (counts.sitesUsed == 0L) {
            myLogger.warn("No sites were shared by the two VCFs. Check that they use the same contig " +
                    "names and the same reference coordinates.")
        } else if (percentUsed < 50.0) {
            myLogger.warn("Only ${"%.1f".format(percentUsed)}% of the sample VCF's sites are in the " +
                    "panel. The path rests on that fraction alone.")
        }
    }

    /** One BED per sample, matching what `impute-path-from-ps4g` writes so `bed-to-vcf` reads either. */
    private fun writeBedFiles(
        paths: Map<String, List<PathInterval<Pair<String, String>>>>,
        outputDir: File
    ) {
        val haploid = pathType == "haploid"
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
            myLogger.info("Wrote $file: ${sampleIntervals.size} intervals")
        }
    }

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
