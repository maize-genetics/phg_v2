package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * This will export the give samples dataset to a hvcf file.
 */
class ExportHvcf : CliktCommand() {

    private val myLogger = LogManager.getLogger(ExportHvcf::class.java)

    val dbpath by option(help = "Folder name where TileDB datasets are stored")
        .required()

    val sampleNames by option(help = "Sample names to export")
        .required()

    val outputDir by option("-o", "--outputDir", help = "Directory where temporary and final files will be written")
        .required()

    override fun run() {

        // This is the tiledbvcf command we want to run:
        // Doing this with a ProcessBuilder and using the phg_v2 conda environment
        // tiledbvcf export --uri tiledb/hvcf_dataset -O v --sample-names Ref --output-dir exported-vcfs

        val builder = ProcessBuilder(
            "conda",
            "run",
            "-n",
            "phgv2-conda",
            "tiledbvcf",
            "export",
            "--uri",
            "$dbpath/hvcf_dataset",
            "-O",
            "v",
            "--sample-names",
            sampleNames,
            "--output-dir",
            outputDir
        )
        val redirectError = "$outputDir/export_hvcf_error.log"
        val redirectOutput = "$outputDir/export_hvcf_output.log"
        builder.redirectOutput(File(redirectOutput))
        builder.redirectError(File(redirectError))

        myLogger.info("ExportHvcf Command: " + builder.command().joinToString(" "))
        val process = builder.start()
        val error = process.waitFor()
        if (error != 0) {
            myLogger.error("tiledbvcf export for: $sampleNames run via ProcessBuilder returned error code $error")
            throw IllegalStateException("Error running tiledbvcf export for: $sampleNames. error: $error")
        }

    }

}