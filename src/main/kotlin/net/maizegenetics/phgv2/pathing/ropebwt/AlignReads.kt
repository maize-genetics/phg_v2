package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.cli.logCommand
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * This command serves as a high-level wrapper around the ropebwt3 mem algorithm,
 * providing a streamlined interface for aligning short sequencing reads to a
 * pre-built FM-index reference.
 *
 * This command identifies Super-Maximal Exact Matches (SMEMs) between query reads
 * and the indexed reference, efficiently mapping reads to genomic coordinates.
 * Results are exported in BED format for easy visualization and downstream PS4G
 * construction.
 */
class AlignReads : CliktCommand(help = "Align reads to ropebwt3 index using maximal exact matches") {

    private val myLogger = LogManager.getLogger(AlignReads::class.java)

    val index by option(
        help = "Path to the ropebwt3 FM-index (.fmd) file used as the reference for read alignment."
    ).required()

    val query by option(
        help = "Input FASTQ file containing reads to align. Supports both uncompressed and .gz-compressed files."
    ).required()

    val minSmemLen by option(
        "--min-smem-len",
        help = "Minimum SMEM (Super-Maximal Exact Match) length used as a seed for alignment. Larger values produce fewer but more specific matches."
    ).int()
        .default(31)

    val threads by option(
        help = "Number of threads to use for parallel processing during alignment. Improves performance on multicore systems."
    ).int()
        .default(1)

    val output by option(
        "-o", "--output",
        help = "Output file in BED format containing aligned read intervals, mapping quality, and strand information."
    ).required()

    val condaEnvPrefix by option(
        help = "Prefix for the conda environment to use. If provided, this should be the full path to the conda environment."
    ).default("")

    override fun run() {
        logCommand(this)

        myLogger.info("Validating inputs...")
        validateInputs()

        myLogger.info("Running ropebwt3 mem alignment...")
        myLogger.info("  Index: $index")
        myLogger.info("  Query: $query")
        myLogger.info("  Min SMEM length: $minSmemLen")
        myLogger.info("  Threads: $threads")
        myLogger.info("  Output: $output")

        runRopebwt3Mem()

        myLogger.info("Alignment complete. Results written to: $output")
    }

    /**
     * Validate that required input files exist and are accessible.
     */
    private fun validateInputs() {
        if (!File(index).exists()) {
            throw IllegalArgumentException("Index file does not exist: $index")
        }
        if (!File(query).exists()) {
            throw IllegalArgumentException("Query file does not exist: $query")
        }

        // Ensure output directory exists
        val outputFile = File(output)
        outputFile.parentFile?.mkdirs()
    }

    /**
     * Execute the ropebwt3 mem command using ProcessBuilder.
     *
     * Command structure: ropebwt3 mem -t<threads> -l<minSmemLen> <index> <query> > <output>
     */
    private fun runRopebwt3Mem() {
        val prefixArg = if (condaEnvPrefix.isNotBlank()) {
            Pair("-p", condaEnvPrefix)
        } else {
            Pair("-n", "phgv2-conda")
        }

        val ropebwt3Command = mutableListOf(
            "conda", "run", prefixArg.first, prefixArg.second,
            "ropebwt3", "mem",
            "-t$threads",
            "-l$minSmemLen",
            index,
            query
        )

        myLogger.info("Executing command: ${ropebwt3Command.joinToString(" ")} > $output")

        val processBuilder = ProcessBuilder(ropebwt3Command)
            .redirectOutput(File(output))
            .redirectError(ProcessBuilder.Redirect.INHERIT)

        val process = processBuilder.start()
        val exitCode = process.waitFor()

        if (exitCode != 0) {
            throw IllegalStateException("ropebwt3 mem command failed with exit code: $exitCode")
        }
    }
}
