package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import net.maizegenetics.phgv2.utils.verifyURI
import org.apache.logging.log4j.LogManager
import java.io.File

class Hvcf2GvcfHandleMultipleSamples :
    CliktCommand(help = "Create g.vcf file for a PHG pathing h.vcf using data from existing PHG created g.vcf files") {

    private val myLogger = LogManager.getLogger(Hvcf2Gvcf::class.java)

    val hvcfDir by option(help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        .required()

    val condaEnvPrefix by option(help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
        .default("")

    val outputDir by option(help = "Output directory for the gVCF files.  If not provided, the current working directory is used.")
        .default("")

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        .default("")

    val referenceFile by option(help = "Path to local Reference FASTA file needed for sequence dictionary")
        .required()

    override fun run() {

        logCommand(this)

        val dbPath = dbPath.ifBlank {
            System.getProperty("user.dir")
        }

        // Verify the tiledbURI - verifyURI will throw an exception if the URI is not valid
        verifyURI(dbPath, "gvcf_dataset", condaEnvPrefix)

        runBlocking {
            exportAllGvcfFiles(outputDir, dbPath, condaEnvPrefix)
        }

    }

    private suspend fun exportAllGvcfFiles(
        outputDir: String,
        dbPath: String,
        condaEnvPrefix: String
    ): Boolean {

        val uri = "$dbPath/gvcf_dataset"
        val allSamples = LoadVcf().getTileDBSampleLists(uri)

        val processingChannel = Channel<Deferred<Boolean>>(25)

        CoroutineScope(Dispatchers.IO).launch {

            allSamples.chunked(60).forEach { sampleList ->
                processingChannel.send(async {
                    exportGvcfFiles(
                        outputDir,
                        dbPath,
                        sampleList,
                        condaEnvPrefix
                    )
                })
            }

            processingChannel.close()

        }

        for (deferred in processingChannel) {
            try {
                deferred.await()
            } catch (e: Exception) {
                myLogger.error("Error exporting gvcf files: ${e.message}")
                return false
            }
        }

        return true

    }

    // tiledbvcf export \
    // --uri my_vcf_dataset \
    // --output-format v \
    // --output-dir exported-vcfs
    private fun exportGvcfFiles(
        outputDir: String,
        dbPath: String,
        samples: List<String>,
        condaEnvPrefix: String
    ): Boolean {

        // Set up the conda environment portion of the command
        val condaCommandPrefix = if (condaEnvPrefix.isNotBlank()) {
            mutableListOf("conda", "run", "-p", condaEnvPrefix)
        } else {
            mutableListOf("conda", "run", "-n", "phgv2-tiledb")
        }

        val dataCommand = mutableListOf(
            "tiledbvcf",
            "export",
            "--uri",
            "$dbPath/gvcf_dataset",
            "--output-format",
            "v",
            "--sample-names",
            samples.joinToString(","),
            "--output-dir",
            outputDir
        )

        // Join the conda and data portions of the command
        val command = condaCommandPrefix + dataCommand
        val builder = ProcessBuilder(command)

        val redirectError = "$outputDir/export_gvcf_error.log"
        val redirectOutput = "$outputDir/export_gvcf_output.log"
        builder.redirectOutput(File(redirectOutput))
        builder.redirectError(File(redirectError))

        // Log the command
        myLogger.info("Command: " + builder.command().joinToString(" "))

        // Start the process
        val process = builder.start()
        val error = process.waitFor()

        // Handle error if the process fails
        if (error != 0) {
            myLogger.error("tiledbvcf export returned error code $error")
            throw IllegalStateException(
                "Error running tiledbvcf export of dataset $dbPath/gvcf_dataset: $error"
            )
        }

        return true

    }

}