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
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.choice
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.pathing.MostLikelyPs4gParents
import net.maizegenetics.phgv2.pathing.PathInputFile
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

    val nParents by option(help = "Restrict the number of parents used for imputation to this number. " +
            "Default = 0 will use all parents.")
        .int()
        .default(0)

    val presenceFile by option(help = "Optional per-founder anchor presence table built from a " +
            "ropebwt3 lift file by `phg build-presence-table`. Used by diploid paths only. Where a " +
            "founder has no alignable sequence in a window, every read matches the other founder " +
            "alone and the model would otherwise read that as evidence for homozygosity; supplying " +
            "this table lets the heterozygous state tie with the homozygous one instead, leaving " +
            "the transition prior to decide.")
        .default("")

    val pavThreshold by option(help = "A founder whose anchor presence in a window is at or below " +
            "this fraction is treated as absent there, for --presence-file. Validated against " +
            "gVCF deletion calls: 0.02 flags windows that are genuinely deleted about 90% of the " +
            "time, and 0.05 about 88%. The looser value scores marginally better overall but results in " +
            "more homozygous calls, so 0.02 is the default. Default = 0.02.")
        .double()
        .default(0.02)

    val pavDamping by option(help = "Strength of the presence/absence correction, 0 to 1. A value of 0 gives " +
            "no correction. A value of 1 is full correction: when one founder is absent, the probability of a homozygote " +
            "equals the probability of a heterozygote . Default = 1.0.")
        .double()
        .default(1.0)

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

    private var loadedPresenceTable: PresenceTable? = null

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

        if (presenceFile.isNotBlank()) {
            require(File(presenceFile).exists()) { "--presence-file $presenceFile does not exist." }
            myLogger.info("Loading founder presence table $presenceFile")
            loadedPresenceTable = PresenceTable(File(presenceFile))
            myLogger.info("Presence table window size = ${loadedPresenceTable!!.windowSize}, " +
                    "${loadedPresenceTable!!.nTaxa} taxa, PAV threshold = $pavThreshold")
        }

        if (isHaploid) imputeHaploidPath(pathToOutputDir)
        else imputeDiploidPath(pathToOutputDir)

    }

    /**
     * Chooses the candidate parent gametes, honoring --n-parents. Shared by both path types so a
     * single-path imputation restricts its candidates exactly as a diploid one does; the former
     * standalone haploid path ignored the option and always used every gamete.
     */
    private fun selectParents(ps4gReader: Ps4gFileReader, contigs: List<String>): Set<Int> {
        val numberOfGenomes = ps4gReader.gameteIndexMap().size
        myLogger.info("Getting parent set for $nParents parents.")
        return if (nParents in 1..<numberOfGenomes) {
            MostLikelyPs4gParents(ps4gReader, contigs.toSet()).bestParents(nParents)
        } else ps4gReader.gameteIndexMap().keys
    }

    /**
     * Imputes a single (haploid) haplotype path for each ps4g file supplied via --path-keyfile
     * or --read-file, writing it to <sampleName>_imputed_path.bed in [outputDir] as
     * chrom/start/end/parent1 records.
     *
     * This is [ViterbiHMM.findDiploidPath] at an inbreeding coefficient of 1, not a separate
     * algorithm. At that coefficient every transition into a heterozygous state has probability
     * zero, so only the homozygous states are reachable, the recursion is the haploid one over
     * the candidate parents, and the emission reduces to the diagonal of the gameteset model.
     * Measured over fifty inbred read sets -- twenty-five in panel, twenty-five held out -- this
     * agrees with the former standalone haploid path to a mean of 0.005 percentage points, and is
     * identical on twenty-one of them; the residual is the binomial coefficient the old haploid
     * emission carried, which does not matter in practice.
     *
     * Running one model rather than two is the point: a single-path imputation can no longer
     * silently diverge from a diploid imputation of homozygous material.
     */
    fun imputeHaploidPath(outputDir: Path) = imputePaths(
        outputDir = outputDir,
        header = "chrom\tstart\tend\tparent1\n",
        // Both founders of a reachable state are the same, so one column is written; the check
        // guards against that invariant being broken upstream rather than papering over it.
        formatCall = {
            check(it.first == it.second) {
                "haploid path produced a heterozygous call (${it.first}, ${it.second})"
            }
            it.first
        },
        parentSelector = ::selectParents,
        inbreedingCoefficient = 1.0
    )

    /**
     * Imputes a pair of (diploid) haplotype paths for each ps4g file supplied via --path-keyfile
     * or --read-file. For every sample, the parent set is optionally reduced to the [nParents] most
     * likely parents via [MostLikelyPs4gParents], then each contig is run through
     * [ViterbiHMM.findDiploidPath] at an inbreeding coefficient of 0, so the two haplotypes switch
     * independently. The resulting paths are written to <sampleName>_imputed_path.bed
     * in [outputDir] as chrom/start/end/parent1/parent2 records.
     *
     * The parent pair is ordered, so adjacent bins are only merged when parent1 and parent2 both
     * match; a phase switch between (A, B) and (B, A) starts a new record.
     */
    fun imputeDiploidPath(outputDir: Path) = imputePaths(
        outputDir = outputDir,
        header = "chrom\tstart\tend\tparent1\tparent2\n",
        formatCall = { "${it.first}\t${it.second}" },
        parentSelector = ::selectParents,
        inbreedingCoefficient = 0.0
    )

    /**
     * Runs the imputation pipeline shared by the haploid and diploid paths: for every input ps4g
     * file, finds a path through each contig — restricted to [contigsToUse] when that option is
     * supplied — and writes the result to <sampleName>_imputed_path.bed in [outputDir].
     *
     * Each contig is run through [ViterbiHMM.findDiploidPath] with transitions from [binTransitions],
     * and the per-bin path is converted to reference coordinates by [pathToIntervals], which merges
     * adjacent bins with equal calls unless --expand-bins is set.
     *
     * @param header the output header line, including its trailing newline.
     * @param formatCall renders a bin's ordered (parent1, parent2) call as the tab-delimited parent
     *   column(s) of a record. Adjacent bins merge when their calls are `==`.
     * @param parentSelector chooses the candidate parent gamete indices given the reader and the
     *   contigs being imputed.
     * @param inbreedingCoefficient 1 for a haploid path, 0 for a diploid one.
     */
    /**
     * Per-haplotype log transition probabilities for the bins of a contig that have reads, in
     * position order: each step is charged for the base pairs between consecutive bins
     * ((bin index difference) x --bin-size) at --prob-switch per --prob-switch-distance, so the
     * switch penalty does not change with bin size or read coverage.
     */
    private fun binTransitions(binPositions: Collection<Int>, nParents: Int): Pair<DoubleArray, DoubleArray> {
        val basePositions = binPositions.sorted().map { bin ->
            val position = bin.toLong() * binSize
            check(position <= Int.MAX_VALUE) { "bin $bin at --bin-size $binSize is beyond ${Int.MAX_VALUE} bp" }
            position.toInt()
        }.toIntArray()
        return transitionLogsByDistance(stepDistances(basePositions), probSwitch, probSwitchDistance, nParents)
    }

    private fun imputePaths(
        outputDir: Path,
        header: String,
        formatCall: (Pair<String, String>) -> String,
        parentSelector: (Ps4gFileReader, List<String>) -> Set<Int>,
        inbreedingCoefficient: Double
    ) {
        val hmm = ViterbiHMM(inbreedingCoefficient, 1.0 - probSwitch, probCorrect, binSize,
            loadedPresenceTable, pavThreshold, pavDamping)
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
                    val (lnNoSwitch, lnSwitch) = binTransitions(readMapForContig.keys, parentSet.size)
                    val contigPath = hmm.findDiploidPath(contig, ps4gReader.gameteIndexMap(),
                        readMapForContig, parentSet, lnNoSwitch, lnSwitch)
                        .map { Pair(it.first, Pair(it.second, it.third)) }
                    myLogger.info("elapsed time for $contig was ${(System.nanoTime() - startTime) / 1_000_000_000.0} sec")

                    pathToIntervals(contigPath, binSize, mergeAdjacent = !expandBins).forEach { interval ->
                        writer.write("${interval.contig}\t${interval.start}\t${interval.end}\t${formatCall(interval.call)}\n")
                    }
                }
            }
        }
    }
}