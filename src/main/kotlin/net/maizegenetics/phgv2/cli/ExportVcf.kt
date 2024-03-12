package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.options.convert
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import net.maizegenetics.phgv2.utils.verifyURI
import org.apache.logging.log4j.LogManager
import java.io.File



// Enum indicating type of input for sample names
enum class SampleFormatEnum(private val value: String) {
    FILE("File"), LIST("List");
}

/**
 * Sealed classes to handle either a list of sample names or a file containing sample names, 1 per line.
 */
sealed class SampleProcessing {
   // abstract fun getSampleOption(): Pair<SampleFormatEnum,String>
    // Would be nice if this could return the command itself, but
    // the command requires the value of the other input parameters,
    // which are not known here.
    abstract fun getExportCommand(): List<String>
    data class SampleFile(val sampleFile: String): SampleProcessing() {
        @Override
        override fun getExportCommand(): List<String> {
            check(File(sampleFile).exists()) { "Samples file $sampleFile does not exist." }
            return listOf(SampleFormatEnum.FILE.toString(), sampleFile)
        }
    }

    data class SampleList(val sampleList: String): SampleProcessing() {
        @Override
        override fun getExportCommand(): List<String> {
            return listOf(SampleFormatEnum.LIST.toString(), sampleList)
        }
    }
}

/**
 * This will export the given samples dataset to an hvcf file.
 */
class ExportVcf : CliktCommand(help = "Export given samples to an h.vcf file") {

    private val myLogger = LogManager.getLogger(ExportVcf::class.java)

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        .default("")

    val datasetType by option(help = "Type of dataset to export: choices are gvcf or hvcf, defaults to hvcf")
        .default("hvcf")


    // Get the sample names from either a list or a file
    // the dataclass SampleProcessing.getExportCommand() a list with 2 strings: the type of input (file or list) and the value
    val samples: SampleProcessing by mutuallyExclusiveOptions(
        option("--sample-names", help = "Comma separated list of Sample names to export. Either sampleNames or samplesFile must be provided").convert{ SampleProcessing.SampleList(it) },
        option("--sample-file", help = "Text file with list of sample names to export, one per line. Either sampleNames or samplesFile must be provided").convert{ SampleProcessing.SampleFile(it)}
    ).single().required()


    val outputDir by option("-o", "--outputDir", help = "Directory where temporary and final files will be written")
        .required()

    override fun run() {

        val dbPath = if (dbPath.isBlank()) {
            System.getProperty("user.dir")
        } else {
            dbPath
        }

        // Verify the tiledbURI - an exception is thrown from verifyURI if the URI is not valid
        val validDB = verifyURI(dbPath,"hvcf_dataset")

        // This is the tiledbvcf command we want to run:
        // Doing this with a ProcessBuilder and using the phg_v2 conda environment
        // tiledbvcf export --uri tiledb/hvcf_dataset -O v --sample-names Ref --output-dir exported-vcfs
        // or ... if the user provides a file with sample names
        // tiledbvcf export --uri tiledb/hvcf_dataset -O v --samples-file sampleNames.txt --output-dir exported-vcfs

        val dtype =  if (datasetType == "gvcf") "gvcf_dataset" else "hvcf_dataset"

        // Tiledbvcf can take either a file with samplenames, or a comma-separated list of sample names
        // setup the command based on user input type.
        var command = if (samples.getExportCommand()[0] == SampleFormatEnum.FILE.toString()) arrayOf(
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
            "--samples-file",
            samples.getExportCommand()[1],
            "--output-dir",
            outputDir
        ) else arrayOf(
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
            samples.getExportCommand()[1],
            "--output-dir",
            outputDir
        )


        val builder = ProcessBuilder(*command)

        val redirectError = "$outputDir/export_${dtype}_error.log"
        val redirectOutput = "$outputDir/export_$dtype}_output.log"
        builder.redirectOutput(File(redirectOutput))
        builder.redirectError(File(redirectError))

        myLogger.info("ExportVcf Command: " + builder.command().joinToString(" "))
        val process = builder.start()
        val error = process.waitFor()
        if (error != 0) {
            myLogger.error("tiledbvcf export for: $samples run via ProcessBuilder returned error code $error")
            throw IllegalStateException("Error running tiledbvcf export of dataset $dbPath/$dtype for: $samples. error: $error")
        }

    }

}