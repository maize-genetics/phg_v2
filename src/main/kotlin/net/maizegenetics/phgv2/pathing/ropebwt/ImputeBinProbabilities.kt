package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.convert
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
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
 * built from:
 *  - an initial state distribution derived from [inbreedCoef] and the number of parents,
 *  - a transition matrix parameterized by [probSame] (haploid) or [DiploidTransitionProbability] (diploid),
 *  - emission probabilities from [EmissionProbabilityForForwardBackward] driven by [probCorrect].
 *
 * Results are written per sample to `<outputDir>/<sampleName>_imputed_probabilities.txt`.
 */
class ImputeBinProbabilities: CliktCommand(help = "Impute best haplotypes from a Ps4g file.")  {

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
        help = "The directory where the imputed assembly haplotypes will be written for each sample. " +
                "File names will be <sampleName>_imputed_path.bed. Coordinates are 0-based half-open."
    )
        .required()

    val imputeType by option(
        help = "The type of path to find. Must be lower case 'haploid' or 'diploid' (without quotes). " +
                "'haploid' infers a single path through the graph. 'diploid' infers a pair of paths. Required parameter."
    )
        .choice("haploid", "diploid")
        .default("haploid")

    val probCorrect by option(help = "The probability that a read maps to correct haplotype. Default = 0.98")
        .double()
        .default(0.98)

    val probSame by option(
        help = "The probability that a path stays on the same gamete when transitioning between " +
                "two adjacent positions. (1 - probability of a recombination). Default = 0.9999"
    )
        .double()
        .default(0.9999)

    val inbreedCoef by option(
        help = "The inbreeding coefficient (between 0.0 and 1.0). " +
                "This parameter is used only for diploid paths. Default = 0.0"
    )
        .double()
        .default(0.0)

    val nParents by option(help = "Restrict the number of parents used for imputation to this number. " +
            "Default = 0 will use all parents.")
        .int()
        .default(0)

    val minProbRatio by option(help = "At each position only probabilities > min-prob-ratio * max probability will be reported. Default = 0.05")
        .double()
        .default(0.05)

    val binSize by option(help = "The size of the bin used for imputation.")
        .int()
        .default(256)

    val myLogger = LogManager.getLogger(ImputeBinProbabilities::class.java)
    val pSwitch = 1.0 - probSame

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

        //Get or create the output directory
        val pathToOutputDir = Paths.get(outputDir)
        pathToOutputDir.createDirectories()

        imputeProbabilities()

    }

    /**
     * Runs the forward-backward HMM for every input PS4G file and writes posterior probabilities
     * to disk.
     *
     * For each sample this:
     *  1. Reads the PS4G file and collects its contigs, excluding any whose name starts with `scaf`.
     *  2. Determines the parent (state) set — restricted to the [nParents] most likely parents via
     *     [MostLikelyPs4gParents] when `nParents > 0`, otherwise every gamete present in the file.
     *  3. Builds the [initialStateProbs], transition matrix (haploid: [probSame]/[pSwitch];
     *     diploid: log transition probabilities from [DiploidTransitionProbability]), and, per contig,
     *     the emission probabilities from [EmissionProbabilityForForwardBackward].
     *  4. Runs [PositionalForwardBackward] over each contig's positions and writes the result to
     *     `<outputDir>/<sampleName>_imputed_probabilities.txt`.
     *
     * @throws IllegalArgumentException if no input read files are provided.
     */
    fun imputeProbabilities() {
        val keyFileLines = readInputFiles.getReadFiles()
        require(keyFileLines.isNotEmpty()) { "Must provide either --path-keyfile or --read-files." }
        val isHaploid = imputeType == "haploid"

        val initialStateProbs = initialStateProbabilityArray()

        for (fileData in keyFileLines) {
            myLogger.info("Finding $imputeType probabilities for ${fileData.sampleName}")
            val ps4gReader = Ps4gFileReader(fileData.file1)

            //do not use contigs starting with scaf
            val contigs = ps4gReader.contigSet().filter { !it.startsWith("scaf") }
            myLogger.info("Contigs: $contigs")

            //if nParents > 0 find most likely parents, otherwise use the full parent set
            val parentSet = if (nParents > 0) {
                MostLikelyPs4gParents(ps4gReader, contigs.toSet()).bestParents(nParents)
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

            val transitionMatrix = transitionMatrix(numberOfParents)

            val numberOfStates = if (isHaploid) numberOfParents else numberOfParents * numberOfParents
            val outputFilepath = Paths.get(outputDir).resolve("${fileData.sampleName}_imputed_probabilities.txt")
            val header = "contig\tposition\tstate\tprobability"

            getBufferedWriter(outputFilepath.toFile()).use { writer ->
                writer.write(header)

                for (contig in contigs) {

                    //Generate list of (index, gamete name) in order by index
                    val readMapForContig = ps4gReader.readMapForContig(contig)
                    check(readMapForContig != null) { "read data for contig $contig was null for ${fileData.sampleName}" }

                    val numberOfPositions = readMapForContig.size

                    val emissionProbability = EmissionProbabilityForForwardBackward(readMapForContig, parentSet, probCorrect)
                    val emissionArray = if (isHaploid) {
                        {ndx: Int -> emissionProbability.getHaploidEmissionProbabilityArray(ndx)}
                    } else {
                        {ndx: Int -> emissionProbability.getDiploidEmissionProbabilityArray(ndx)}
                    }

                    val imputedResult = PositionalForwardBackward(numberOfStates, numberOfPositions,
                        initialStateProbs, transitionMatrix, emissionArray).run()

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

    private fun initialStateProbabilityArray(): DoubleArray {
        return if (imputeType == "haploid") {
            DoubleArray(nParents) {1.0/nParents}
        } else {
            val homozygoteProbabillity = inbreedCoef / nParents
            val heterozygoteProbability = (1.0 - inbreedCoef) /(nParents * nParents - nParents)
            DoubleArray(nParents * nParents) {ndx -> if (ndx % nParents == ndx / nParents) homozygoteProbabillity
            else heterozygoteProbability}
        }
    }

    private fun transitionMatrix(numberOfParents: Int): DoubleArray {
        return if (imputeType == "haploid") {
            DoubleArray(numberOfParents * numberOfParents) {ndx -> if (ndx % numberOfParents == ndx / numberOfParents) probSame else pSwitch }
        } else {
            val transitionMatrix = DoubleArray(numberOfParents * numberOfParents)
            var ptr = 0
            for (index1 in 0 until numberOfParents) {
                for (index2 in 0 until numberOfParents) {
                    for (index3 in 0 until numberOfParents) {
                        for (index4 in 0 until numberOfParents) {
                            transitionMatrix[ptr++] = DiploidTransitionProbability(probSame, inbreedCoef, numberOfParents)
                                .calculateLn(Pair(index1, index2), Pair(index3, index4))
                        }
                    }
                }
            }
            transitionMatrix
        }
    }

}