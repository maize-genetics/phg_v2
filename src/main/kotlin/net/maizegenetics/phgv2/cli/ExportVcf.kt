package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * This will export the give samples dataset to a hvcf file.
 */
class ExportVcf : CliktCommand(help = "Export given samples to an h.vcf file") {

    private val myLogger = LogManager.getLogger(ExportVcf::class.java)

    val dbPath by option(help = "Folder name where TileDB datasets are stored")
        .required()

    val datasetType by option(help = "Type of dataset to export: choices are gvcf or hvcf, defaults to hvcf")
        .default("hvcf")

    val sampleNames by option(help = "Comma separated list of Sample names to export")
        .required()

    val outputDir by option("-o", "--outputDir", help = "Directory where temporary and final files will be written")
        .required()

    override fun run() {

        // This is the tiledbvcf command we want to run:
        // Doing this with a ProcessBuilder and using the phg_v2 conda environment
        // tiledbvcf export --uri tiledb/hvcf_dataset -O v --sample-names Ref --output-dir exported-vcfs

        val dtype =  if (datasetType == "gvcf") "gvcf_dataset" else "hvcf_dataset"
        val builder = ProcessBuilder(
            "conda",
            "run",
            "-n",
            "phgv2-conda",
            "tiledbvcf",
            "export",
            "--uri",
            "$dbPath/$dtype",
            "-O",
            "v",
            "--sample-names",
            sampleNames,
            "--output-dir",
            outputDir
        )
        val redirectError = "$outputDir/export_${dtype}_error.log"
        val redirectOutput = "$outputDir/export_$dtype}_output.log"
        builder.redirectOutput(File(redirectOutput))
        builder.redirectError(File(redirectError))

        myLogger.info("ExportVcf Command: " + builder.command().joinToString(" "))
        val process = builder.start()
        val error = process.waitFor()
        if (error != 0) {
            myLogger.error("tiledbvcf export for: $sampleNames run via ProcessBuilder returned error code $error")
            throw IllegalStateException("Error running tiledbvcf export of dataset $dbPath/$dtype for: $sampleNames. error: $error")
        }

    }

}