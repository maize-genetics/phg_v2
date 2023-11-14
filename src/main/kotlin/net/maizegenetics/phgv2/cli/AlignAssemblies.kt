package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.int
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.launch
import kotlinx.coroutines.runBlocking
import kotlinx.coroutines.withContext
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.stream.Collectors

/**
 * This will align assemblies to a reference genome.
 * Uses anchorwave's proali to align, which handles
 * genome alignment with relocation variation,
 * chromosome fusion or whole genome duplication.
 *
 * This function allows for multiple assembly alignments to be run
 * in parallel.  Users may specify number of alignments to run in parallel
 * and the total number of threads available to be split between the alignments.
 * Anchorwave takes 10G of RAM for the dynamic program aspects, and then another
 * 10G+ per thread, depending on the processor.
 *
 * The table below shows the memory usage for a single assembly alignment
 * based on processor type: (From Baoxing Song, 2023-11-10)
 *
 *  ---------------------------------------------
 *  | Processor | peak memory (Gb) | wall time  |
 *  ---------------------------------------------
 *  |  SSE2    |    20.1          | 26:47:17    |
 *  ---------------------------------------------
 *  |  SSE4.1  |    20.6          | 24:05:07    |
 *  ---------------------------------------------
 *  |  AVX2    |    20.1          | 21:40:00    |
 *  ---------------------------------------------
 *  |  AVX512  |    20.1          | 18:31:39    |
 *  ---------------------------------------------
 *  | ARM      |    34.2          | 18:08:57    |
 *  ---------------------------------------------
 *
 *  For a machine with 128G RAM, you would use 10 threads
 *  divided between your parallel runs.  For example, if you
 *  have 10 assemblies to align, you could run 2 parallel
 *  alignments, each using 5 threads.
 *
 */
class AlignAssemblies : CliktCommand(help="Align assemblies using anchorwave") {

    private val myLogger = LogManager.getLogger(AlignAssemblies::class.java)

    val gff by option(help = "Full path to the reference gff file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--gff must not be blank"
            }
        }

    val reference by option(help = "Full path to reference fasta file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--reference must not be blank"
            }
        }

    val assemblies by option(
        "-a",
        "--assemblies",
        help = "File containing list of assemblies to align, 1 per line, full path to file"
    )
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--assemblies must not be blank"
            }
        }

    val outputDir by option("-o", "--output-dir", help = "Directory where temporary and final files will be written")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--output-dir must not be blank"
            }
        }

    val totalThreads by option(help = "Number of threads available.  These will be split among the alginments that are run in parallel")
        .int()
        .default(1)

    val inParallel by option(
        help = "Number of assemblies to simultaneously process. " +
                "If you have 10 threads and the in-parallel value is 2, then 2 assemblies will be aligned at a time, each using 5 threads." +
                "The anchorwave application can take up to 30G per thread for each assembly processed, plus some overhead." +
                "Consider this memory factor when providing values for the total-threads and in-parallel."
    )
        .int()
        .default(1)

    val refMaxAlignCov by option(help = "Anchorwave proali parameter R, indicating reference genome maximum alignment coverage.")
        .int()
        .default(1)

    val queryMaxAlignCov by option(help = "Anchorwave proali parameter Q, indicating query genome maximum alignment coverage.")
        .int()
        .default(1)

    data class InputChannelData(
        val refFasta: String,
        val asmFasta: String,
        val outputDir: String,
        val gffFile: String,
        val refSamOutFile: String
    )

    override fun run() {

        // create CDS fasta from reference and gff3 file
        val cdsFasta = "$outputDir/ref.cds.fasta"

        createCDSfromRefData(reference, gff, cdsFasta, outputDir)

        // create list of assemblies to align from the assemblies file
        val assembliesList = File(assemblies).readLines().filter { it.isNotBlank() }

        // run minimap2 for ref to refcds
        val justNameRef = File(reference).nameWithoutExtension
        val samOutFile = "${justNameRef}.sam"
        val refSamOutFile = "${outputDir}/${samOutFile}"

        val builder = ProcessBuilder(
            "conda", "run", "-n", "phgv2-conda", "minimap2", "-x", "splice", "-t", totalThreads.toString(), "-k", "12",
            "-a", "-p", "0.4", "-N20", reference, cdsFasta, "-o", refSamOutFile
        )
        val redirectError = "$outputDir/minimap2Ref_error.log"
        val redirectOutput = "$outputDir/minimap2Ref_output.log"
        builder.redirectOutput(File(redirectOutput))
        builder.redirectError(File(redirectError))

        myLogger.info("Ref minimap Command: " + builder.command().joinToString(" "));

        val process = builder.start()
        val error = process.waitFor()
        if (error != 0) {
            myLogger.error("minimap2 for $reference run via ProcessBuilder returned error code $error")
            throw IllegalStateException("Error running minimap2 for reference: $error")
        }

        runAnchorWaveMultiThread(reference, assembliesList, cdsFasta, gff, refSamOutFile)

    }

    private fun createCDSfromRefData(refFasta: String, gffFile: String, cdsFasta: String, outputDir: String): Boolean {

        // val command = "anchorwave gff2seq -r ${refFasta} -i ${gffFile} -o ${cdsFasta} "
        // Need to set the conda environment here to access anchorwave
        val builder = ProcessBuilder(
            "conda",
            "run",
            "-n",
            "phgv2-conda",
            "anchorwave",
            "gff2seq",
            "-r",
            refFasta,
            "-i",
            gffFile,
            "-o",
            cdsFasta
        )

        val redirectOutput = "$outputDir/anchorwave_gff2seq_output.log"
        val redirectError = "$outputDir/anchorwave_gff2seq_error.log"
        builder.redirectOutput(File(redirectOutput))
        builder.redirectError(File(redirectError))

        myLogger.info("createCDSfromRefData command:" + builder.command().joinToString(" "))
        try {
            val process = builder.start()
            val error = process.waitFor()
            if (error != 0) {
                myLogger.error("createCDSfromRefData run via ProcessBuilder returned error code $error")
                return false
            }
            return true
        } catch (e: Exception) {
            myLogger.error("Error: could not execute anchorwave command. Run anchorwave manually and retry.")
            return false
        }

    }

    /**
     * This will run anchorwave for each assembly in the assemblies list.
     * Each anchorwave process will be run in a separate thread.
     */
    private fun runAnchorWaveMultiThread(
        refFasta: String,
        assemblies: List<String>,
        cdsFasta: String,
        gffFile: String,
        refSamOutFile: String
    ) {
        runBlocking {
            // Setup
            val inputChannel = Channel<InputChannelData>()

            // The input channel gets data needed to run minimap2 and align with anchorwave
            launch {
                myLogger.info("Adding entries to the inputChannel:")
                assemblies.forEach { asmFile ->

                    // Column names were checked for validity above

                    myLogger.info("Adding: $asmFile for processing")
                    inputChannel.send(InputChannelData(refFasta, asmFile, outputDir, gffFile, refSamOutFile))
                }
                myLogger.info("Done Adding data to the inputChannel:")
                inputChannel.close() // Need to close this here to show the workers that it is done adding more data
            }

            // Do not need a coroutine that "joins" the threads as they will all
            // terminate above when there is no more data on the input channel

            // This calls anchorwave's proali, and minimap2 scripts to process the alignments
            val workerThreads = (1..inParallel).map { run ->
                launch { alignAssembly(inputChannel, cdsFasta, gffFile) }
            }

        }
    }

    private suspend fun alignAssembly(inputChannel: Channel<InputChannelData>, cdsFasta: String, gffFile: String) =
        withContext(
            Dispatchers.Default
        ) {
            val threadsPerRun = totalThreads / inParallel
            println("alignAssembly: totalThreads: $totalThreads, inParallel: $inParallel, threadsPerRun: $threadsPerRun")
            for (assemblyEntry in inputChannel) {
                // Column names were checked for validity above
                val justName = File(assemblyEntry.asmFasta).nameWithoutExtension
                val samShort = "${justName}.sam"
                val asmSamFile = "${assemblyEntry.outputDir}/${samShort}"

                myLogger.info("alignAssembly: asmFileFull: ${assemblyEntry.asmFasta}, outputFile: $asmSamFile")


                val builder = ProcessBuilder(
                    "conda",
                    "run",
                    "-n",
                    "phgv2-conda",
                    "minimap2",
                    "-x",
                    "splice",
                    "-t",
                    threadsPerRun.toString(),
                    "-k",
                    "12",
                    "-a",
                    "-p",
                    "0.4",
                    "-N20",
                    assemblyEntry.asmFasta,
                    cdsFasta,
                    "-o",
                    asmSamFile
                )

                val redirectError = "${assemblyEntry.outputDir}/minimap2_${justName}_error.log"
                val redirectOutput = "${assemblyEntry.outputDir}/minimap2_${justName}_output.log"
                myLogger.info("redirectError: $redirectError")
                builder.redirectOutput(File(redirectOutput))
                builder.redirectError(File(redirectError))
                myLogger.info(" begin minimap assembly Command: " + builder.command().joinToString(" "))

                val process = builder.start()
                val error = process.waitFor()
                if (error != 0) {
                    myLogger.error("minimap2 for assembly ${assemblyEntry.asmFasta} run via ProcessBuilder returned error code $error")
                    throw IllegalStateException("alignAssembly: error running minimap2 for ${justName}: $error")
                }
                // We have the SAM File, call proali to align with anchorwave
                runAnchorwaveProali(
                    gffFile,
                    assemblyEntry.refFasta,
                    assemblyEntry.asmFasta,
                    cdsFasta,
                    assemblyEntry.refSamOutFile,
                    asmSamFile,
                    threadsPerRun.toString()
                )
            }
        }

    private fun runAnchorwaveProali(
        gffFile: String,
        refFasta: String,
        asmFasta: String,
        cdsFasta: String,
        refSam: String,
        asmSam: String,
        threadsPerRun: String
    ) {

        val justNameAsm = File(asmFasta).nameWithoutExtension
        val justNameRef = File(refFasta).nameWithoutExtension

        val anchorsproFile = "${outputDir}/${justNameAsm}_${justNameRef}.anchorspro"
        // Using just the assembly name.  This facilitates programmatically creating
        // a GVCF keyfile from the maf keyfiles.  It will be understood that the maf
        // file name is <assemblyFastaNoExtension>.maf
        val outputFile = "${outputDir}/${justNameAsm}.maf"
        val builder = ProcessBuilder(
            "conda",
            "run",
            "-n",
            "phgv2-conda",
            "anchorwave",
            "proali",
            "-i",
            gffFile,
            "-r",
            refFasta,
            "-as",
            cdsFasta,
            "-a",
            asmSam,
            "-ar",
            refSam,
            "-s",
            asmFasta,
            "-n",
            anchorsproFile,
            "-R",
            refMaxAlignCov.toString(),
            "-Q",
            queryMaxAlignCov.toString(),
            "-t",
            threadsPerRun,
            "-o",
            outputFile
        )

        val redirectError = "${outputDir}/proali_${justNameAsm}_outputAndError.log"
        myLogger.info("redirectError: $redirectError")

        // NOTE: anchowave proali has an output file parameter, unlike minimap2, which uses
        // command line redirection (ie >) to write the .SAM file
        //
        // Do not use the anchorwave output file parameter as the message redirect output
        // for ProcessBuilder  or you'll end up with the text "AnchorWave done!" in the middle of
        // a reference sequence, causing problems in parsing, and causing the sequence itself to be
        // on a separate line from the data describing it.  Instead, send all program output messages
        // to the error file.
        // Here, we  write all program message output (ie non MAF file output) to a single file
        // named as above.

        builder.redirectOutput(File(redirectError))
        builder.redirectError(File(redirectError))
        myLogger.info(
            "runAnchorwaveProali proali Command for ${justNameAsm}: " + builder.command().joinToString (" ")
        )

        try {
            val process = builder.start()
            val error = process.waitFor()
            if (error != 0) {
                myLogger.error("proali for assembly $asmFasta run via ProcessBuilder returned error code $error")
                throw IllegalStateException("runAnchorwaveProali: error running proali for $justNameAsm")
            }
        } catch (e: Exception) {
            myLogger.error("Error: could not execute anchorwave command. Run anchorwave manually and retry.")
        }
    }

}