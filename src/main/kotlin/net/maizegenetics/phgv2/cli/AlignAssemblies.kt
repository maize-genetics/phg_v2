package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
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
 */
class AlignAssemblies : CliktCommand() {

    private val myLogger = LogManager.getLogger(AlignAssemblies::class.java)

    val gff by option(help = "Full path to the reference gff file")
        .required()

    val ref by option(help = "Full path to reference fasta file")
        .required()

    val assemblies by option(
        "-a",
        "--assemblies",
        help = "File containing list of assemblies to align, 1 per line, full path to file"
    )
        .required()

    val outputDir by option("-o", "--outputDir", help = "Directory where temporary and final files will be written")
        .required()

    val threads by option(help = "Number of threads to use for each assembly processed")
        .int()
        .default(1)

    val runs by option(
        help = "Number of assemblies to simultaneously process. " +
                "If you have 4 threads and 2 runs, then 2 assemblies will be processed at a time, each using 4 threads." +
                "The anchorwave application can take up to 50G per thread for each assembly processed, plus some overhead." +
                "Consider this memory factor when providing values for threadsPerRun and numRuns."
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

        createCDSfromRefData(ref, gff, cdsFasta, outputDir)

        // create list of assemblies to align from the assemblies file
        val assembliesList = File(assemblies).readLines().filter { it.isNotBlank() }

        // run minimap2 for ref to refcds
        val justNameRef = File(ref).nameWithoutExtension
        val samOutFile = "${justNameRef}.sam"
        val refSamOutFile = "${outputDir}/${samOutFile}"

        val builder = ProcessBuilder(
            "conda", "run", "-n", "phgv2-conda", "minimap2", "-x", "splice", "-t", threads.toString(), "-k", "12",
            "-a", "-p", "0.4", "-N20", ref, cdsFasta, "-o", refSamOutFile
        )
        val redirectError = "$outputDir/minimap2Ref_error.log"
        val redirectOutput = "$outputDir/minimap2Ref_output.log"
        builder.redirectOutput(File(redirectOutput))
        builder.redirectError(File(redirectError))

        myLogger.info("Ref minimap Command: " + builder.command().stream().collect(Collectors.joining(" ")));
        val process = builder.start()
        val error = process.waitFor()
        if (error != 0) {
            myLogger.error("minimap2 for $ref run via ProcessBuilder returned error code $error")
            throw IllegalStateException("Error running minimap2 for reference: $error")
        }

        runAnchorWaveMultiThread(ref, assembliesList, cdsFasta, gff, refSamOutFile)

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

        myLogger.info("createCDSfromRefData command:" + builder.command().stream().collect(Collectors.joining(" ")))
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

            // This calls mummer4 scripts to process the alignments
            val workerThreads = (1..runs).map { run ->
                launch { alignAssembly(inputChannel, cdsFasta, gffFile) }
            }

        }
    }

    private suspend fun alignAssembly(inputChannel: Channel<InputChannelData>, cdsFasta: String, gffFile: String) =
        withContext(
            Dispatchers.Default
        ) {
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
                    threads.toString(),
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
                myLogger.info(
                    " begin minimap assembly Command: " + builder.command().stream().collect(Collectors.joining(" "))
                );
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
                    asmSamFile
                )
            }
        }

    private fun runAnchorwaveProali(
        gffFile: String,
        refFasta: String,
        asmFasta: String,
        cdsFasta: String,
        refSam: String,
        asmSam: String
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
            threads.toString(),
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
            "runAnchorwaveProali proali Command for ${justNameAsm}: " + builder.command().stream()
                .collect(Collectors.joining(" "))
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