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
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * Raises the density of a sample VCF: infers which founders each sample descends from, then reads its
 * genotypes back out of a denser panel.
 *
 * Two panels, doing two different jobs. The first has to overlap the sample's own sites, because that
 * is the evidence the path is inferred from; the second supplies the output's sites, so the denser it
 * is the more the output gains. They may be the same file, in which case nothing is gained but nothing
 * breaks.
 *
 *     sample VCF + panel VCF  ->  founder path  ->  high-density panel  ->  imputed VCF
 *
 * Both halves already exist as commands -- `impute-path-from-vcf` and `bed-to-vcf` -- and this calls the
 * same code they do. What it adds is that **nothing intermediate reaches disk**: the paths are handed
 * across in memory. `--bed-dir` writes them anyway, for looking at which founders were chosen, which is
 * expected to be the exception rather than the rule.
 *
 * ## What is checked before any work starts
 *
 * Founder names have to line up: the path names founders from the first panel, and the second panel is
 * asked for those founders' alleles. Where a name is missing the composition yields a no-call and says
 * nothing, so a mismatch would quietly produce a VCF of no-calls. The first panel's samples must
 * therefore be a **subset** of the second's, checked up front and reported by name. A superset is fine
 * -- extra founders in the dense panel are simply never asked for.
 *
 * ## Assumptions, inherited
 *
 * Both panels are haploid or mostly homozygous, and the sample VCF is unphased. Note the two halves
 * treat a heterozygous founder differently, each in the way that suits its job: path inference uses
 * both of its alleles under Mendelian sampling, while composition, which has to name one allele and
 * cannot tell which was inherited, emits a no-call. A site where one founder of a pair is heterozygous
 * therefore comes out half called, as `1/.`. See [BedToVcf].
 *
 * ## What the output covers
 *
 * A path spans the first to the last site the sample shares with `--panel-vcf`, so dense-panel sites
 * outside that span get no call. `--extend-to-contig-ends` carries the terminal ancestry out to both
 * ends of each contig instead. See [applyContigBounds].
 */
class ImputeVcfFromVcf : CliktCommand(help = "Impute a higher-density VCF for the samples in a VCF, via an inferred founder path") {

    private val myLogger = LogManager.getLogger(ImputeVcfFromVcf::class.java)

    val toImputeVcf by option(help = "VCF holding the samples to impute. Coordinate sorted and " +
            "unphased. Required parameter.")
        .required()
        .validate { require(File(it).exists()) { "$it is not a valid file" } }

    val panelVcf by option(help = "Reference panel VCF used to infer the path. Its sites must overlap " +
            "the sample VCF's, since that overlap is the evidence, and its samples must be a subset of " +
            "the high-density panel's. Required parameter.")
        .required()
        .validate { require(File(it).exists()) { "$it is not a valid file" } }

    val highDensityPanelVcf by option(help = "Reference panel VCF supplying the output's sites and " +
            "alleles. Denser than --panel-vcf is the point, though the same file is allowed. Required " +
            "parameter.")
        .required()
        .validate { require(File(it).exists()) { "$it is not a valid file" } }

    val outputFile by option(help = "The imputed VCF to write. Required parameter.")
        .required()

    val bedDir by option(help = "Optional directory for the intermediate founder paths, as " +
            "<sampleName>_imputed_path.bed. Nothing intermediate is written unless this is given; " +
            "supply it to see which founders were chosen.")
        .default("")

    val pathType by option(help = "The type of path to find. 'haploid' infers a single founder per " +
            "position, which is a diploid path with the inbreeding coefficient set to 1. 'diploid' " +
            "infers a pair.")
        .choice("haploid", "diploid")
        .default("haploid")

    val probCorrect by option(help = "The probability that a genotype call is correct. Default = 0.98")
        .double()
        .default(0.98)
        .validate { require(it > 0.5 && it < 1.0) { "prob-correct must be between 0.5 and 1.0 exclusive" } }

    val probSwitch by option(help = "The probability of a path switch (a recombination) across " +
            "--prob-switch-distance, scaled to each step's actual distance. Default = 1e-4")
        .double()
        .default(1e-4)
        .validate { require(it > 0.0 && it < 1.0) { "prob-switch must be between 0 and 1 exclusive" } }

    val probSwitchDistance by option(help = "The distance --prob-switch is quoted over, in base " +
            "pairs. Default = 1000000")
        .double()
        .default(1_000_000.0)
        .validate { require(it > 0.0) { "prob-switch-distance must be positive" } }

    val inbreedCoef by option(help = "The inbreeding coefficient, used for diploid paths. Only 0.0 " +
            "and 1.0 are supported; --path-type haploid sets it to 1.0 regardless. Default = 0.0")
        .double()
        .default(0.0)
        .validate {
            require(it == 0.0 || it == 1.0) {
                "inbreed-coef must be 0.0 or 1.0 for this command; $it would need the general Viterbi " +
                        "scan, which cannot take a distance-scaled transition"
            }
        }

    val extendToContigEnds by option(help = "Carry the first and last interval of each contig out to " +
            "the contig's ends. Off by default, so a path spans only the first to the last site shared " +
            "with the panel and a denser panel's sites outside that span get no call. There is no " +
            "evidence of ancestry beyond the terminal markers; set this to assume it continues.")
        .flag()

    val contigsToUse by option(help = "Comma-separated contigs to impute, or a file with one per " +
            "line. All contigs shared by the sample VCF and --panel-vcf if omitted.")
        .default("")

    override fun run() {
        logCommand(this)
        requireFoundersPresentInDensePanel()

        val pathFinder = ImputePathFromVcf()
        val parameters = VcfPathParameters(
            pathType = pathType,
            inbreedCoef = inbreedCoef,
            probCorrect = probCorrect,
            probSwitch = probSwitch,
            probSwitchDistance = probSwitchDistance,
            contigsToUse = pathFinder.buildContigSet(contigsToUse),
            extendToContigEnds = extendToContigEnds
        )

        val result = imputeFounderPaths(File(toImputeVcf), File(panelVcf), parameters, myLogger)
        logPairingCounts(result.counts, myLogger)
        if (result.paths.isEmpty()) throw UsageError(
            "No founder path could be inferred, so there is nothing to compose. The two VCFs shared " +
                    "no sites; check that they use the same contig names and reference coordinates."
        )

        if (bedDir.isNotBlank()) {
            writeFounderPathBeds(result.paths, File(bedDir), pathType == "haploid", myLogger)
        } else {
            myLogger.info("Holding ${result.paths.size} founder paths in memory; pass --bed-dir to " +
                    "write them out")
        }

        val composer = BedToVcf()
        composer.writeVcf(File(highDensityPanelVcf), composer.pathsFromIntervals(result.paths),
            File(outputFile))
    }

    /**
     * Fails, before any work, if the path panel names founders the dense panel does not carry.
     *
     * This is the mismatch worth catching early: composition looks each founder up by name in the dense
     * panel and emits a no-call when it is absent, so a naming difference between the two panels
     * produces a complete VCF of no-calls with nothing in the log to say why.
     */
    private fun requireFoundersPresentInDensePanel() {
        val pathFounders = VCFFileReader(File(panelVcf), false).use { it.fileHeader.sampleNamesInOrder.toList() }
        val denseFounders = VCFFileReader(File(highDensityPanelVcf), false)
            .use { it.fileHeader.sampleNamesInOrder.toSet() }
        val missing = pathFounders.filterNot { denseFounders.contains(it) }
        if (missing.isNotEmpty()) throw UsageError(
            "${missing.size} of ${pathFounders.size} founders in --panel-vcf are absent from " +
                    "--high-density-panel-vcf, so a path naming them could not be composed: " +
                    "${missing.take(10).joinToString(", ")}" +
                    (if (missing.size > 10) ", and ${missing.size - 10} more" else "") +
                    ". The path panel's samples must be a subset of the high-density panel's."
        )
        myLogger.info("${pathFounders.size} founders, all present in the high-density panel " +
                "(${denseFounders.size} samples)")
    }
}
