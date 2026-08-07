package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.UsageError
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.convert
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.choice
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.pathing.MostLikelyPs4gParents
import net.maizegenetics.phgv2.pathing.PathInputFile
import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.getBufferedReader
import net.maizegenetics.phgv2.utils.getBufferedWriter
import org.apache.logging.log4j.LogManager
import java.io.File
import java.nio.file.InvalidPathException
import java.nio.file.Path
import java.nio.file.Paths
import kotlin.io.path.createDirectories
import kotlin.io.path.exists

class ImputePathFromPs4g: CliktCommand(help = "Impute best haplotypes from a Ps4g file.") {

    companion object {
        /**
         * Parses the value of --contigs-to-use into the set of contig names to be imputed.
         *
         * If [contigString] names an existing file, the contigs are read from it, one per line.
         * Otherwise the value is treated as a comma-separated list of contig names. A blank value
         * returns an empty set, which [Ps4gFileReader] takes to mean "use every contig in the file".
         */
        fun buildContigSet(contigString: String): Set<String> {
            return if (contigString.isNotBlank()) {
                try {
                    val filePath = Path.of(contigString)
                    if(filePath.exists()) {
                        val tempContigSet = getBufferedReader(filePath.toFile())
                            .use { it.readLines()}.map {it.trim()}.filter{it.isNotBlank()}.toSet()
                        if (tempContigSet.isEmpty()) throw UsageError(
                            "--contigs-to-use names the file $contigString, which contains no contig names."
                        )
                        tempContigSet
                    }
                    else splitContigList(contigString)
                } catch (e: InvalidPathException) {
                    //the value cannot be a file name, so it must be a list of contigs
                    splitContigList(contigString)
                }
            } else emptySet()
        }

        /**
         * Splits a comma-separated list of contig names, trimming each name and dropping any that
         * are blank.
         */
        private fun splitContigList(contigString: String): Set<String> =
            contigString.split(",").map { it.trim() }.filter { it.isNotBlank() }.toSet()
    }

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

    val outPathDir by option(
        help = "The directory where the imputed assembly haplotypes will be written for each sample. " +
                "File names will be <sampleName>_imputed_path.bed. Coordinates are 0-based half-open."
    )
        .required()

    val pathType by option(
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

    val nParents by option(help = "Restrict the number of parents used for diploid imputation to this number. " +
            "Default = 0 will use all parents.")
        .int()
        .default(0)

    val binSize by option(help = "The bin size used to create the ps4g file. Default = 256.")
        .int()
        .default(256)

    val expandBins by option(help = "Write one output record per bin instead of merging adjacent bins " +
            "that have identical parents into a single record. Default = false (adjacent bins are merged).")
        .flag()

    val contigsToUse by option(help = "A list of contigs to be imputed. If no list is supplied, all contigs in the ps4g" +
            " file will be imputed. The value can be a comma-separated list or a file containing the list," +
            " with a single contig per line.")
        .default("")

    val myLogger = LogManager.getLogger(ImputePathFromPs4g::class.java)

    /**
     * Entry point for the command. Creates the output directory and dispatches to either
     * [imputeHaploidPath] or [imputeDiploidPath] depending on the value of [pathType].
     */
    override fun run() {

        logCommand(this)

        val maxMemory = Runtime.getRuntime().maxMemory()
        println("Max memory is ${maxMemory/1024/1024}mb")

        //create the outParentsDir, if it does not already exist
        if (outPathDir.isNotBlank()) File(outPathDir).mkdirs()

        val isHaploid = pathType == "haploid"

        //Get or create the output directory
        val pathToOutputDir = Paths.get(outPathDir)
        pathToOutputDir.createDirectories()

        if (isHaploid) imputeHaploidPath(pathToOutputDir)
        else imputeDiploidPath(pathToOutputDir)

    }

    /**
     * Imputes a single (haploid) haplotype path for each ps4g file supplied via --path-keyfile
     * or --read-file. For every sample, each contig is run through [ViterbiHMM.findHaploidPath]
     * and the resulting path is written to <sampleName>_imputed_path.bed in [outputDir] as
     * chrom/start/end/parent1 records.
     */
    fun imputeHaploidPath(outputDir: Path) = imputePaths<String>(
        outputDir = outputDir,
        header = "chrom\tstart\tend\tparent1\n",
        formatCall = { it },
        parentSelector = { ps4gReader, _ -> ps4gReader.gameteIndexMap().keys },
        pathFinder = { hmm, contig, gameteIndexMap, readMap, parentSet ->
            hmm.findHaploidPath(contig, gameteIndexMap, readMap, parentSet)
        }
    )

    /**
     * Imputes a pair of (diploid) haplotype paths for each ps4g file supplied via --path-keyfile
     * or --read-file. For every sample, the parent set is optionally reduced to the [nParents] most
     * likely parents via [MostLikelyPs4gParents], then each contig is run through
     * [ViterbiHMM.findDiploidPath]. The resulting paths are written to <sampleName>_imputed_path.bed
     * in [outputDir] as chrom/start/end/parent1/parent2 records.
     *
     * The parent pair is ordered, so adjacent bins are only merged when parent1 and parent2 both
     * match; a phase switch between (A, B) and (B, A) starts a new record.
     */
    fun imputeDiploidPath(outputDir: Path) = imputePaths<Pair<String, String>>(
        outputDir = outputDir,
        header = "chrom\tstart\tend\tparent1\tparent2\n",
        formatCall = { "${it.first}\t${it.second}" },
        parentSelector = { ps4gReader, contigs ->
            //if the number of likely parents is > 0 and < number of genomes, find the likely parents
            val numberOfGenomes = ps4gReader.gameteIndexMap().size
            myLogger.info("Getting parent set for $nParents parents.")
            if (nParents in 1..<numberOfGenomes) {
                MostLikelyPs4gParents(ps4gReader, contigs.toSet()).bestParents(nParents)
            } else ps4gReader.gameteIndexMap().keys
        },
        pathFinder = { hmm, contig, gameteIndexMap, readMap, parentSet ->
            hmm.findDiploidPath(contig, gameteIndexMap, readMap, parentSet)
                .map { Pair(it.first, Pair(it.second, it.third)) }
        }
    )

    /**
     * Runs the imputation pipeline shared by the haploid and diploid paths: for every input ps4g
     * file, finds a path through each contig — restricted to [contigsToUse] when that option is
     * supplied — and writes the result to <sampleName>_imputed_path.bed in [outputDir].
     *
     * The per-bin path returned by [pathFinder] is converted to reference coordinates by
     * [pathToIntervals], which merges adjacent bins with equal calls unless --expand-bins is set.
     *
     * @param T the parent call for one bin: the parent name for haploid paths, an ordered
     *   (parent1, parent2) pair for diploid paths. Adjacent bins merge when their calls are `==`.
     * @param header the output header line, including its trailing newline.
     * @param formatCall renders a call as the tab-delimited parent column(s) of a record.
     * @param parentSelector chooses the candidate parent gamete indices given the reader and the
     *   contigs being imputed.
     * @param pathFinder runs the Viterbi algorithm for one contig, returning (position, call) pairs
     *   ordered by bin position.
     */
    private fun <T> imputePaths(
        outputDir: Path,
        header: String,
        formatCall: (T) -> String,
        parentSelector: (Ps4gFileReader, List<String>) -> Set<Int>,
        pathFinder: (ViterbiHMM, String, Map<Int, String>, Map<Int, MutableList<Ps4gGameteSet>>, Set<Int>) -> List<Pair<Position, T>>
    ) {
        val keyFileLines = readInputFiles.getReadFiles()
        require(keyFileLines.isNotEmpty()) { "Must provide either --path-keyfile or --read-files." }

        for (fileData in keyFileLines) {
            myLogger.info("Finding $pathType path for ${fileData.sampleName}")
            val ps4gReader = Ps4gFileReader(fileData.file1, buildContigSet(contigsToUse))

            //the reader has already dropped any contigs not in --contigs-to-use
            val contigs = ps4gReader.contigSet()
            myLogger.info("Contigs: $contigs")

            val parentSet = parentSelector(ps4gReader, contigs.sorted())
            myLogger.info("Parent set: $parentSet")

            val outputFilepath = outputDir.resolve("${fileData.sampleName}_imputed_path.bed")
            getBufferedWriter(outputFilepath.toFile()).use { writer ->

                writer.write(header)

                for (contig in contigs) {

                    //Generate list of (index, gamete name) in order by index
                    val readMapForContig = ps4gReader.readMapForContig(contig)
                    check(readMapForContig != null) { "read data for contig $contig was null for ${fileData.sampleName}" }

                    val startTime = System.nanoTime()
                    val contigPath = pathFinder(
                        ViterbiHMM(inbreedCoef, probSame, probCorrect),
                        contig, ps4gReader.gameteIndexMap(), readMapForContig, parentSet
                    )
                    myLogger.info("elapsed time for $contig was ${(System.nanoTime() - startTime) / 1_000_000_000.0} sec")

                    pathToIntervals(contigPath, binSize, mergeAdjacent = !expandBins).forEach { interval ->
                        writer.write("${interval.contig}\t${interval.start}\t${interval.end}\t${formatCall(interval.call)}\n")
                    }
                }
            }
        }
    }
}