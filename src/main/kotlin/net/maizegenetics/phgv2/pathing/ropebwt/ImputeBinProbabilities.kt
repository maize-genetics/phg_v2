package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.UsageError
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.check
import com.github.ajalt.clikt.parameters.options.convert
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.choice
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.pathing.MostLikelyPs4gParents
import net.maizegenetics.phgv2.pathing.PathInputFile
import net.maizegenetics.phgv2.utils.getBufferedWriter
import org.apache.logging.log4j.LogManager
import java.io.File
import java.nio.file.Paths
import kotlin.io.path.createDirectories

/**
 * CLI command that computes posterior haplotype (bin) probabilities from PS4G read-mapping data
 * using the forward-backward algorithm.
 *
 * Unlike [ImputePathFromPs4g], which uses Viterbi to infer a single most-likely path, this command
 * runs a positional forward-backward HMM to produce per-position posterior probabilities over the
 * candidate haplotype states. It supports both `haploid` imputation (states are individual parents)
 * and `diploid` imputation (states are ordered parent pairs).
 *
 * Input is one or more PS4G files produced by `align-reads`, supplied either through a tab-delimited
 * key file (`--path-keyfile`) or directly (`--read-file`). For each sample, a per-contig HMM is
 * built from the same model as [ImputePathFromPs4g], so the posteriors describe the model whose
 * most likely path that command reports:
 *  - a uniform initial state distribution,
 *  - transitions from [SwitchTransitionStep], charged for the base pairs between consecutive bins
 *    with reads at [probSwitch] per [probSwitchDistance]; a diploid path's two haplotypes switch
 *    independently,
 *  - or, for a diploid path with [selfingGenerations], the transitions and initial distribution of a
 *    line selfed from a biparental F1 ([SelfingTransitions], applied by [SelfingTransitionStep]),
 *    at the same distance-scaled rate,
 *  - emission probabilities from [GameteSetEmissionProbability] driven by [probCorrect], with the
 *    presence/absence correction when [presenceFile] is supplied (diploid only).
 *
 * Results are written per sample to `<outputDir>/<sampleName>_imputed_probabilities.txt`.
 */
class ImputeBinProbabilities: CliktCommand(help = "Calculate posterior haplotype probabilities from a Ps4g file using an HMM.")  {

    val readInputFiles: PathInputFile by mutuallyExclusiveOptions<PathInputFile>(
        option(
            "--path-keyfile", help = "Name of tab-delimited key file.  Columns for samplename and filename" +
                    " are required. Files must be ps4g files created by align-reads. A value must be entered for " +
                    "either --path-keyfile or --read-files."
        )
            .convert { PathInputFile.KeyFile(it) },
        option("--read-file", help = "The name of a ps4g file created by align-reads.")
            .convert { PathInputFile.ReadFiles(it) }
    ).single().required()

    val outputDir by option(
        help = "The directory where the haplotype probabilities will be written for each sample. " +
                "File names will be <sampleName>_imputed_probabilities.txt."
    )
        .required()

    val imputeType by option(
        help = "The type of genotype probabilities to calculate. Must be lower case 'haploid' or 'diploid' (without quotes)."
    )
        .choice("haploid", "diploid")
        .default("haploid")

    val probCorrect by option(help = "The probability that a read maps to correct haplotype. Default = 0.98")
        .double()
        .default(0.98)
        .check { it in 0.5..1.0 }

    val probSwitch by option(
        help = "The probability of a path switch (a recombination) across --prob-switch-distance, " +
                "scaled to each step's actual distance: the base pairs between consecutive bins with " +
                "reads. A smoothing parameter more than a recombination rate, so it does not change " +
                "meaning with --bin-size or with read coverage. Default = 1e-16"
    )
        .double()
        .default(1e-16)
        .validate { require(it > 0.0 && it < 1.0) { "prob-switch must be between 0 and 1 exclusive" } }

    val probSwitchDistance by option(help = "The distance --prob-switch is quoted over, in base " +
            "pairs. Default = 1000000")
        .double()
        .default(1_000_000.0)
        .validate { require(it > 0.0) { "prob-switch-distance must be positive" } }

    val selfingGenerations by option(help = "For diploid imputation, the number of generations of " +
            "selfing after a biparental F1 (S5 = F6 is 5). Replaces the model in which the two haplotypes " +
            "switch independently with the exact transitions of a selfed line, as in impute-path-from-ps4g, " +
            "so the posteriors give the probability of each homozygous and heterozygous state of a " +
            "recombinant inbred line. --prob-switch still sets the overall smoothing. Takes at most " +
            "${SelfingTransitions.MAX_PARENTS} parents, so use it with --n-parents (normally 2). " +
            "Default: not set, the independent-haplotype model.")
        .int()
        .validate { require(it >= 1) { "selfing-generations must be at least 1" } }

    val nParents by option(help = "Restrict the number of parents used for imputation to this number. " +
            "Default = 0 will use all parents.")
        .int()
        .default(0)

    val presenceFile by option(help = "Optional per-founder anchor presence table built from a " +
            "ropebwt3 lift file by `phg build-presence-table`. Used by diploid imputation only. Where a " +
            "founder has no alignable sequence in a window, every read matches the other founder " +
            "alone and the model would otherwise read that as evidence for homozygosity; supplying " +
            "this table lets the heterozygous state tie with the homozygous one instead, leaving " +
            "the transition prior to decide.")
        .default("")

    val pavThreshold by option(help = "A founder whose anchor presence in a window is at or below " +
            "this fraction is treated as absent there, for --presence-file. Default = 0.02.")
        .double()
        .default(0.02)

    val pavDamping by option(help = "Strength of the presence/absence correction, 0 to 1. A value of 0 gives " +
            "no correction. A value of 1 is full correction: when one founder is absent, the probability of a homozygote " +
            "equals the probability of a heterozygote. Default = 1.0.")
        .double()
        .default(1.0)
        .validate { require(it in 0.0..1.0) { "pav-damping must be between 0 and 1" } }

    val minProbRatio by option(help = "At each position only probabilities > min-prob-ratio * max probability will be reported. Default = 0.05")
        .double()
        .default(0.05)

    val binSize by option(help = "The size of the bin used for imputation.")
        .int()
        .default(256)

    val contigsToUse by option(help = "A list of contigs to be imputed. If no list is supplied, all contigs in the ps4g" +
            " file will be imputed. The value can be a comma-separated list or a file containing the list," +
            " with a single contig per line.")
        .default("")

    val myLogger = LogManager.getLogger(ImputeBinProbabilities::class.java)

    private var loadedPresenceTable: PresenceTable? = null

    /**
     * Clikt entry point. Logs the invocation, reports available JVM memory, ensures the
     * [outputDir] exists, and then delegates the imputation work to [imputeProbabilities].
     */
    override fun run() {
        logCommand(this)

        val maxMemory = Runtime.getRuntime().maxMemory()
        println("Max memory is ${maxMemory/1024/1024}mb")

        //create the outParentsDir, if it does not already exist
        if (outputDir.isNotBlank()) File(outputDir).mkdirs()

        if (selfingGenerations != null) {
            if (imputeType == "haploid") throw UsageError("--selfing-generations applies to --impute-type diploid only.")
            if (nParents == 0) throw UsageError("--selfing-generations takes at most " +
                    "${SelfingTransitions.MAX_PARENTS} parents but --n-parents = 0, which uses all parents. " +
                    "Set --n-parents to a value > 0 (normally 2).")
            if (nParents > SelfingTransitions.MAX_PARENTS) throw UsageError("--selfing-generations takes at most " +
                    "${SelfingTransitions.MAX_PARENTS} parents; set --n-parents (normally 2).")
        }

        //Get or create the output directory
        val pathToOutputDir = Paths.get(outputDir)
        pathToOutputDir.createDirectories()

        if (presenceFile.isNotBlank()) {
            require(File(presenceFile).exists()) { "--presence-file $presenceFile does not exist." }
            if (imputeType == "haploid") myLogger.warn("--presence-file is ignored for haploid imputation")
            myLogger.info("Loading founder presence table $presenceFile")
            loadedPresenceTable = PresenceTable(File(presenceFile))
            myLogger.info("Presence table window size = ${loadedPresenceTable!!.windowSize}, " +
                    "${loadedPresenceTable!!.nTaxa} taxa, PAV threshold = $pavThreshold")
        }

        imputeProbabilities()

    }

    /**
     * Runs the forward-backward HMM for every input PS4G file and writes posterior probabilities
     * to disk.
     *
     * For each sample this:
     *  1. Reads the PS4G file and collects its contigs, restricted to [contigsToUse] when that
     *     option is supplied.
     *  2. Determines the parent (state) set — restricted to the [nParents] most likely parents via
     *     [MostLikelyPs4gParents] when `nParents > 0`, otherwise every gamete present in the file.
     *  3. Builds the initial state probabilities -- uniform, or the selfed line's genotype
     *     frequencies with [selfingGenerations] -- and, per contig, the distance-scaled transitions
     *     ([SwitchTransitionStep], or [SelfingTransitionStep]) and the emission probabilities
     *     ([GameteSetEmissionProbability]: the homozygous diagonal for a haploid path, every ordered
     *     pair for a diploid one, with the presence/absence correction when a presence table is loaded).
     *  4. Runs [PositionalForwardBackward] over each contig's positions and writes the result to
     *     `<outputDir>/<sampleName>_imputed_probabilities.txt`.
     *
     * @throws IllegalArgumentException if no input read files are provided.
     */
    fun imputeProbabilities() {
        val keyFileLines = readInputFiles.getReadFiles()
        require(keyFileLines.isNotEmpty()) { "Must provide either --path-keyfile or --read-files." }
        val isHaploid = imputeType == "haploid"

        for (fileData in keyFileLines) {
            myLogger.info("Finding $imputeType probabilities for ${fileData.sampleName}")
            val ps4gReader = Ps4gFileReader(fileData.file1, ImputePathFromPs4g.buildContigSet(contigsToUse))

            //the reader has already dropped any contigs not in --contigs-to-use
            val contigs = ps4gReader.contigSet()
            myLogger.info("Contigs: $contigs")

            //if nParents > 0 find most likely parents, otherwise use the full parent set
            val parentSet = if (nParents > 0) {
                MostLikelyPs4gParents(ps4gReader, contigs).bestParents(nParents)
            } else {
                ps4gReader.gameteIndexMap().keys
            }

            //create a sorted parent list from the set
            val parentList = parentSet.sorted()

            //get the states names for the parent list, because that will be needed to write the output
            val gameteIndexMap = ps4gReader.gameteIndexMap()
            val stateNames = if (isHaploid) {
                parentList.map { gameteIndexMap[it] }
            } else {
                parentList.flatMap { parent1 -> parentList.map { parent2 -> "(${gameteIndexMap[parent1]},${gameteIndexMap[parent2]})" } }

            }

            val numberOfParents = parentSet.size
            myLogger.info("Parent set: $parentSet")

            val numberOfStates = if (isHaploid) numberOfParents else numberOfParents * numberOfParents
            // run() has already refused --selfing-generations for haploid imputation.
            val selfing = selfingGenerations?.let { SelfingTransitions(it, numberOfParents) }
            // Uniform, as in ViterbiHMM: no state is favored before the reads are seen. A selfed line
            // starts from its genotype frequencies, as impute-path-from-ps4g does.
            val initialStateProbs = selfing?.initialProbabilities()
                ?: DoubleArray(numberOfStates) { 1.0 / numberOfStates }
            val outputFilepath = Paths.get(outputDir).resolve("${fileData.sampleName}_imputed_probabilities.txt")
            val header = "contig\tposition\tstate\tprobability\n"

            getBufferedWriter(outputFilepath.toFile()).use { writer ->
                writer.write(header)

                for (contig in contigs) {

                    //Generate list of (index, gamete name) in order by index
                    val readMapForContig = ps4gReader.readMapForContig(contig)
                    check(readMapForContig != null) { "read data for contig $contig was null for ${fileData.sampleName}" }

                    val numberOfPositions = readMapForContig.size

                    // the presence/absence correction only rewrites heterozygous cells, so a haploid path never uses it
                    val emissionProbability = GameteSetEmissionProbability(readMapForContig, parentSet, probCorrect,
                        if (isHaploid) null else loadedPresenceTable, contig, gameteIndexMap, binSize,
                        pavThreshold, pavDamping)
                    val lnEmissionArray = if (isHaploid) {
                        {ndx: Int -> emissionProbability.getHomozygousEmissionProbabilityArray(ndx)}
                    } else {
                        {ndx: Int -> emissionProbability.getDiploidEmissionProbabilityArray(ndx)}
                    }

                    val distances = binDistances(readMapForContig.keys)
                    val transition = if (selfing != null) {
                        SelfingTransitionStep(selfing, distances, probSwitch, probSwitchDistance)
                    } else {
                        val (lnNoSwitch, lnSwitch) = transitionLogsByDistance(distances,
                            probSwitch, probSwitchDistance, numberOfParents)
                        SwitchTransitionStep(numberOfParents, !isHaploid, lnNoSwitch, lnSwitch)
                    }

                    val imputedResult = PositionalForwardBackward(numberOfStates, numberOfPositions,
                        initialStateProbs, transition, lnEmissionArray).run()

                    val positionProbabilities = imputedResult.posteriorStateProbabilities
                    val positionList = readMapForContig.keys.sorted()
                    positionList.mapIndexed { index, positionIndex ->
                        val probabilityArray = positionProbabilities[index]
                        val maxProb = probabilityArray.max()
                        val minProb = maxProb * minProbRatio
                        for (stateIndex in 0 until numberOfStates) {
                            if (probabilityArray[stateIndex] > minProb) {
                                val start = binSize * positionIndex
                                writer.write("$contig\t$start\t${stateNames[stateIndex]}\t${probabilityArray[stateIndex]}")
                                writer.newLine()

                            }
                        }
                    }

                }
            }
        }
    }

    /**
     * Base pairs between consecutive bins with reads, in position order, as [ImputePathFromPs4g]
     * charges them; see [stepDistances].
     */
    private fun binDistances(binPositions: Collection<Int>): IntArray {
        val basePositions = binPositions.sorted().map { bin ->
            val position = bin.toLong() * binSize
            check(position <= Int.MAX_VALUE) { "bin $bin at --bin-size $binSize is beyond ${Int.MAX_VALUE} bp" }
            position.toInt()
        }.toIntArray()
        return stepDistances(basePositions)
    }

}
