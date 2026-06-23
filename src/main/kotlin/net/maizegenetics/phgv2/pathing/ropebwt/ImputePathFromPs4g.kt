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
import net.maizegenetics.phgv2.pathing.KeyFileData
import net.maizegenetics.phgv2.pathing.PathFinderHMMPS4G
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
        .default(0.95)

    val probSame by option(
        help = "The probability that a path stays on the same gamete when transitioning between " +
                "two adjacent positions. (1 - probability of a recombination)"
    )
        .double()
        .default(0.001)

    val inbreedCoef by option(
        help = "The inbreeding coefficient (between 0.0 and 1.0). " +
                "This parameter is used only for diploid paths."
    )
        .double()
        .default(0.0)

    val likelyParents by option(help = "Restrict the number of parents used for diploid imputation to this number. " +
            "Default = 0 will use all parents.")
        .int()
        .default(0)

    val myLogger = LogManager.getLogger(ImputePathFromPs4g::class.java)

    override fun run() {

        logCommand(this)

        //create the outParentsDir, if it does not already exist
        if (outPathDir.isNotBlank()) File(outPathDir).mkdirs()
        val outputDir = Paths.get(outPathDir)

        val isHaploid = pathType == "haploid"

        //Get or create the output directory
        val pathToOutputDir = Paths.get(outPathDir)
        pathToOutputDir.createDirectories()

        if (isHaploid) imputeHaploidPath(pathToOutputDir)
        else imputeDiploidPath(pathToOutputDir)

    }

    fun imputeHaploidPath(outputDir: Path) {
        val keyFileLines = readInputFiles.getReadFiles()
        require(keyFileLines.isNotEmpty()) { "Must provide either --path-keyfile or --read-files." }


        val pathFinder = PathFinderHMMPS4G(probCorrect, probSame, inbreedCoef)


        for (fileData in keyFileLines) {
            myLogger.info("Finding $pathType path for ${fileData.sampleName}")
            val ps4gReader = Ps4gFileReader(fileData.file1)
            val contigs = ps4gReader.contigSet()

            val outputFilepath = outputDir.resolve("${fileData.sampleName}_imputed_path.txt")
            outputFilepath.bufferedWriter().use { writer ->
                writer.write("contig\tposition\tgenome\n")
            }

            for (contig in contigs) {

                //Generate list of (index, gamete name) in order by index
                val readMapForContig = ps4gReader.readMapForContig(contig)
                check(readMapForContig != null) { "read data for contig $contig was null for ${fileData.sampleName}" }
                val contigPath = pathFinder.findHaploidPath(contig, ps4gReader.gameteIndexMap(), readMapForContig)

                //write to the output file
                outputFilepath.bufferedWriter(
                    options = arrayOf(StandardOpenOption.APPEND)
                ).use { writer ->
                    contigPath.forEach {
                        writer.write("${it.first.contig}\t${it.first.position}\t${it.second}\n")
                    }
                }

            }
        }


    }

    fun imputeDiploidPath(outputDir: Path) {
        val keyFileLines = readInputFiles.getReadFiles()
        require(keyFileLines.isNotEmpty()) { "Must provide either --path-keyfile or --read-files." }


        val pathFinder = PathFinderHMMPS4G(probCorrect, probSame, inbreedCoef)

        for (fileData in keyFileLines) {
            myLogger.info("Finding $pathType path for ${fileData.sampleName}")
            val ps4gReader = Ps4gFileReader(fileData.file1)
            val contigs = ps4gReader.contigSet()

            val outputFilepath = outputDir.resolve("${fileData.sampleName}_imputed_path.txt")
            outputFilepath.bufferedWriter().use { writer ->
                writer.write("contig\tposition\tgenome1\tgenome2\n")
            }

            for (contig in contigs) {

                //Generate list of (index, gamete name) in order by index
                val readMapForContig = ps4gReader.readMapForContig(contig)
                check(readMapForContig != null) { "read data for contig $contig was null for ${fileData.sampleName}" }
                val contigPath = pathFinder.findDiploidPath(contig, ps4gReader.gameteIndexMap(), readMapForContig)

                //write to the output file
                outputFilepath.bufferedWriter(
                    options = arrayOf(StandardOpenOption.APPEND)
                ).use { writer ->
                    contigPath.forEach {

                        writer.write("${it.first.contig}\t${it.first.position}\t${it.second}\t${it.third}\n")
                    }
                }

            }

        }
    }
}