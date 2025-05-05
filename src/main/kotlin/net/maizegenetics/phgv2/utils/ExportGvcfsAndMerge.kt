package net.maizegenetics.phgv2.utils

import biokotlin.util.MergeGVCFUtils.mergeGVCFs
import htsjdk.variant.vcf.VCFFileReader
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import net.maizegenetics.phgv2.cli.LoadVcf
import org.apache.logging.log4j.LogManager
import java.io.File

private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.utils.ExportGvcfsAndMerge")

const val GVCF_EXPORT_DIR = "gvcf_phg_export_dir"
const val MASTER_GVCF_NAME = "master_gvcf.g.vcf"

/**
 * This function exports gvcf files from a tiledb database and merges them into a master gvcf file.
 * It checks if the master gvcf file exists and if it contains all samples.
 * If it does not, it exports all gvcf files from the tiledb database
 * and merges them into the master gvcf file.
 */
fun exportGvcfFilesAndMerge(
    dbPath: String,
    condaEnvPrefix: String = ""
) {

    // Create the gvcf export directory if it does not exist
    val gvcfDir = "$dbPath/$GVCF_EXPORT_DIR"
    if (!File(gvcfDir).exists()) {
        val result = File(gvcfDir).mkdirs()
        if (!result) {
            myLogger.error("Failed to create gvcf export directory: $gvcfDir")
            throw IllegalStateException("Failed to create gvcf export directory: $gvcfDir")
        }
    }

    // This is the master GVCF file that will be created by merging all gvcf files
    // if it does not exist and if it does not contain all samples
    val masterGVCF = "$dbPath/$MASTER_GVCF_NAME"

    // Get a list of all samples in the existing master GVCF file
    val samplesFromMaster = if (File(masterGVCF).exists()) {
        samplesFromVCFFile(masterGVCF)
    } else {
        emptyList()
    }

    // Get list of all samples in the gvcf dataset of the tiledb database
    val gvcfDataset = "$dbPath/gvcf_dataset"
    val allSamples = LoadVcf().getTileDBSampleLists(gvcfDataset)

    // Check if the master GVCF file exists and if it contains all samples
    // If it does not, export missing gvcf files from the tiledb database
    // and merge them into the master GVCF file
    if (!samplesFromMaster.containsAll(allSamples)) {

        // Get list of samples that are missing from the GVCF directory
        val missingSamples = allSamples
            .filter { !File("$gvcfDir/$it.vcf").exists() }

        runBlocking {
            exportGvcfFilesForSamples(missingSamples, gvcfDir, dbPath, condaEnvPrefix)
        }

        val filesToMerge = allSamples
            .map { "$gvcfDir/$it.vcf" }

        File(masterGVCF).delete()
        mergeGVCFs(filesToMerge, masterGVCF)

    }

}

/**
 * This function reads a VCF file and returns a list of sample names.
 */
fun samplesFromVCFFile(filename: String): List<String> {

    val vcfFile = File(filename)
    if (!vcfFile.exists()) {
        throw IllegalArgumentException("samplesFromVCFFile: VCF file does not exist: $filename")
    }

    VCFFileReader(vcfFile, false).use {
        val header = it.header
        val samples = header.sampleNamesInOrder
        return samples
    }

}

/**
 * This function exports all gvcf files from the tiledb database.
 * It uses coroutines to process the export in parallel.
 *
 * @param outputDir The directory where the gvcf files will be exported.
 * @param dbPath The path to the tiledb database.
 * @param condaEnvPrefix The prefix for the conda environment to use. Default is an empty string.
 *
 * @return A list of all samples processed.
 */
private suspend fun exportGvcfFilesForSamples(
    sampleNames: List<String>,
    outputDir: String,
    dbPath: String,
    condaEnvPrefix: String = ""
) {

    // Create the log directory if it does not exist
    val logDir = "$outputDir/logs"
    if (!File(logDir).exists()) {
        File(logDir).mkdirs()
    }

    val processingChannel = Channel<Deferred<Boolean>>(25)

    CoroutineScope(Dispatchers.IO).launch {

        sampleNames.chunked(1).forEach { sampleList ->
            processingChannel.send(async {
                exportGvcfFiles(
                    sampleList,
                    outputDir,
                    dbPath,
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
            throw IllegalStateException("Error exporting gvcf files: ${e.message}")
        }
    }

}

// function to export from tiledb the gvcf files if they don't exist
// If we get the tiledb-java API working for MAC, we can hold these in memory
// while processing and skip the export.
private fun exportGvcfFiles(
    sampleNames: List<String>,
    outputDir: String,
    dbPath: String,
    condaEnvPrefix: String = ""
): Boolean {

    val logDir = "$outputDir/logs"

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
        sampleNames.joinToString(","),
        "--output-dir",
        outputDir
    )

    // Join the conda and data portions of the command
    val command = condaCommandPrefix + dataCommand
    val builder = ProcessBuilder(command)

    val redirectError = "$logDir/export_${sampleNames[0]}_error.log"
    val redirectOutput = "$logDir/export_${sampleNames[0]}_output.log"
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