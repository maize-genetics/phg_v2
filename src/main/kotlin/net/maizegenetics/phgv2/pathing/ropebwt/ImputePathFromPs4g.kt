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
import org.apache.logging.log4j.LogManager
import java.io.File
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.StandardOpenOption
import kotlin.io.path.bufferedWriter
import kotlin.io.path.createDirectories

class ImputePathFromPs4g: CliktCommand(help = "Impute best haplotypes from a Ps4g file.") {

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
                "File names will be <sampleName>_imputed_path.txt. If no directory name is supplied, the files will not be written."
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
     * or --read-file. For every sample, each non-scaffold contig is run through [ViterbiHMM.findHaploidPath]
     * and the resulting path is written to <sampleName>_imputed_path.bed in [outputDir] as
     * chrom/start/end/parent1 records.
     */
    fun imputeHaploidPath(outputDir: Path) {
        val keyFileLines = readInputFiles.getReadFiles()
        require(keyFileLines.isNotEmpty()) { "Must provide either --path-keyfile or --read-files." }

        for (fileData in keyFileLines) {
            myLogger.info("Finding $pathType path for ${fileData.sampleName}")
            val ps4gReader = Ps4gFileReader(fileData.file1)
            val contigs = ps4gReader.contigSet().filter{!it.startsWith("scaf")}
            val parentSet = ps4gReader.gameteIndexMap().keys

            val outputFilepath = outputDir.resolve("${fileData.sampleName}_imputed_path.bed")
            outputFilepath.bufferedWriter().use { writer ->
                writer.write("chrom\tstart\tend\tparent1\n")
            }


            for (contig in contigs) {

                //Generate list of (index, gamete name) in order by index
                val readMapForContig = ps4gReader.readMapForContig(contig)
                check(readMapForContig != null) { "read data for contig $contig was null for ${fileData.sampleName}" }
                val startTime = System.nanoTime()
                val contigPath = ViterbiHMM(inbreedCoef, probSame, probCorrect)
                    .findHaploidPath(contig, ps4gReader.gameteIndexMap(), readMapForContig, parentSet)
                myLogger.info("elapsed time for $contig was ${(System.nanoTime() - startTime) / 1_000_000_000.0} sec")

                //write to the output file
                outputFilepath.bufferedWriter(
                    options = arrayOf(StandardOpenOption.APPEND)
                ).use { writer ->
                    if (contigPath.size > 1) {
                        var startPos = 1
                        var endPos = (contigPath[0].first.position + contigPath[1].first.position) / 2 * binSize
                        writer.write("${contigPath[0].first.contig}\t$startPos\t$endPos\t${contigPath[0].second}\n")

                        (1..(contigPath.size - 2)).forEach { ndx ->

                            startPos =
                                ((contigPath[ndx - 1].first.position + contigPath[ndx].first.position) / 2) * binSize + 1
                            endPos = (contigPath[ndx].first.position + contigPath[ndx + 1].first.position) / 2 * binSize
                            writer.write("${contigPath[ndx].first.contig}\t$startPos\t$endPos\t${contigPath[ndx].second}\n")
                        }

                    }

                    //write the last record
                    val ndx = contigPath.size - 1
                    val startPos = ((contigPath[ndx - 1].first.position + contigPath[ndx].first.position) / 2) * binSize + 1
                    val endPos = contigPath[ndx].first.position * binSize
                    writer.write("${contigPath[ndx].first.contig}\t$startPos\t$endPos\t${contigPath[ndx].second}\n")
                }

            }
        }


    }

    /**
     * Imputes a pair of (diploid) haplotype paths for each ps4g file supplied via --path-keyfile
     * or --read-file. For every sample, the parent set is optionally reduced to the [nParents] most
     * likely parents via [MostLikelyPs4gParents], then each non-scaffold contig is run through
     * [ViterbiHMM.findDiploidPath]. The resulting paths are written to <sampleName>_imputed_path.txt
     * in [outputDir] as chrom/start/end/parent1/parent2 records.
     */
    fun imputeDiploidPath(outputDir: Path) {
        val keyFileLines = readInputFiles.getReadFiles()
        require(keyFileLines.isNotEmpty()) { "Must provide either --path-keyfile or --read-files." }

        for (fileData in keyFileLines) {
            myLogger.info("Finding $pathType path for ${fileData.sampleName}")
            val ps4gReader = Ps4gFileReader(fileData.file1)
            val tmpContigs = ps4gReader.contigSet()

            //do not use contigs starting with scaf
            val contigs = tmpContigs.filter {!it.startsWith("scaf")}
            myLogger.info("Contigs: $contigs")

            val outputFilepath = outputDir.resolve("${fileData.sampleName}_imputed_path.txt")
            outputFilepath.bufferedWriter().use { writer ->
                writer.write("chrom\tstart\tend\tparent1\tparent2\n")
            }

            //if the number of likely parents is > 0 and < number of genomes, find the likely parents
            val numberOfGenomes = ps4gReader.gameteIndexMap().size

            myLogger.info("Getting parent set for $nParents parents.")
            val parentSet = if (nParents in 1..<numberOfGenomes) {
                MostLikelyPs4gParents(ps4gReader, contigs.toSet()).bestParents(nParents)
            } else ps4gReader.gameteIndexMap().keys

            myLogger.info("Parent set: $parentSet")
            for (contig in contigs) {

                //Generate list of (index, gamete name) in order by index
                val readMapForContig = ps4gReader.readMapForContig(contig)
                check(readMapForContig != null) { "read data for contig $contig was null for ${fileData.sampleName}" }

                val startTime = System.nanoTime()
                val contigPath = ViterbiHMM(inbreedCoef, probSame, probCorrect)
                    .findDiploidPath(contig, ps4gReader.gameteIndexMap(), readMapForContig, parentSet)
                myLogger.info("elapsed time for $contig was ${(System.nanoTime() - startTime) / 1_000_000_000.0} sec")

                //write to the output file
                outputFilepath.bufferedWriter(
                    options = arrayOf(StandardOpenOption.APPEND)
                ).use { writer ->
                    var startPos = 1
                    var endPos = (contigPath[0].first.position + contigPath[1].first.position) / 2 * binSize
                    writer.write("${contigPath[0].first.contig}\t$startPos\t$endPos\t${contigPath[0].second}\t${contigPath[0].third}\n")
                    (1..(contigPath.size - 2)).forEach { ndx ->

                        startPos =
                            ((contigPath[ndx - 1].first.position + contigPath[ndx].first.position) / 2) * binSize + 1
                        endPos = (contigPath[ndx].first.position + contigPath[ndx + 1].first.position) / 2 * binSize
                        writer.write("${contigPath[ndx].first.contig}\t$startPos\t$endPos\t${contigPath[ndx].second}\t${contigPath[ndx].third}\n")
                    }
                    //write the last record
                    val ndx = contigPath.size - 1
                    startPos = ((contigPath[ndx - 1].first.position + contigPath[ndx].first.position) / 2) * binSize + 1
                    endPos = contigPath[ndx].first.position * binSize
                    writer.write("${contigPath[ndx].first.contig}\t$startPos\t$endPos\t${contigPath[ndx].second}\t${contigPath[ndx].third}\n")
                }

            }

        }
    }
}