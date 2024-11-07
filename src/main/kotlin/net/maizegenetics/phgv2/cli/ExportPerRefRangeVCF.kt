package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.convert
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.utils.verifyURI
import org.apache.logging.log4j.LogManager
import java.io.File

class ExportPerRefRangeVCF: CliktCommand(help = "Export given samples to merged vcf file on a per-range basis") {

    private val myLogger = LogManager.getLogger(ExportVcf::class.java)

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        .default("")

    // Get the sample names from either a list or a file
    // the dataclass SampleProcessing.getExportCommand() a list with 2 strings: the type of input (file or list) and the value
    val samples: SampleProcessing by mutuallyExclusiveOptions(
        option(
            "--sample-names",
            help = "Comma separated list of Sample names to export. Either --sample-names or --samples-file must be provided"
        ).convert { SampleProcessing.SampleList(it) },
        option(
            "--sample-file",
            help = "Text file with list of sample names to export, one per line. Either --sample-names or --samples-file must be provided"
        ).convert { SampleProcessing.SampleFile(it) }
    ).single().required()


    val outputDir by option("-o", "--output-dir", help = "Directory where temporary and final files will be written")
        .required()

    val regionsFile by option(help = "A bedfile or vcf file containing the regions to be exported. Regions can be single base pair positions. File extension must be either .bed or .vcf.")
        .default("")

    // not currently used - will figure out if this makes sense later.
    val batchSize by option(help= "Number of samples to export at a time.  Default is 5")
        .int()
        .default(5)


    val condaEnvPrefix by option (help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
        .default("")

    override fun run() {
        //TODO("Not yet implemented")

        val dbPath = if (dbPath.isBlank()) {
            System.getProperty("user.dir")
        } else {
            dbPath
        }

        // Verify the tiledbURI - an exception is thrown from verifyURI if the URI is not valid
        val validDB = verifyURI(dbPath, "gvcf_dataset",condaEnvPrefix)
        runExportByRegion(dbPath, outputDir, regionsFile, samples, condaEnvPrefix, batchSize)

    }

    // This function will run the export command for individual regions based on user input and the regions file
    // It should be batched at some point, but for now, it runs each region individually
    fun runExportByRegion(dbPath: String, outputDir: String, regionsFile: String, samples: SampleProcessing,
                          condaEnvPrefix: String, batchSize: Int) {
        // Setup the conda environment portion of the command
        var commandPrefix = if (condaEnvPrefix.isNotBlank()) mutableListOf("conda","run","-p",condaEnvPrefix) else mutableListOf("conda","run","-n","phgv2-conda")

        val regions = convertRegionFileToRegionList(regionsFile)
        val batches = regions.chunked(batchSize) // not sending in batches right now - will figure that out later

        // THis isn't correct - we want to send a process builder command for each entry in
        // the batch, running batch.size number of samples at a time

        regions.forEach { region ->
            // Create the data portion of the command for this batch

            val outputFile = "$outputDir/export_gvcf_${region.replace(":","_")}.vcf"
            var dataCommand = if (samples.getExportCommand()[0] == SampleFormatEnum.FILE.toString()) mutableListOf(
                "tiledbvcf",
                "export",
                "--uri",
                "$dbPath/gvcf_dataset",
                "-m",
                "-b",
                "65536",
                "-O",
                "v",
                "--samples-file",
                samples.getExportCommand()[1],
                "-o",
                outputFile,
                "--regions",
                region,
                "--tiledb-config",
                "sm.memory_budget=10737418240,sm.memory_budget_var=21474836480,sm.skip_unary_partitioning_budget_check=true"
            ) else mutableListOf(
                "tiledbvcf",
                "export",
                "--uri",
                "$dbPath/gvcf_dataset",
                "-m",
                "-b",
                "65536",
                "-O",
                "v",
                "--sample-names",
                samples.getExportCommand()[1],
                "-o",
                outputFile,
                "--regions",
                region,
                "--tiledb-config",
                "sm.memory_budget=10737418240,sm.memory_budget_var=21474836480,sm.skip_unary_partitioning_budget_check=true"
            )
            // Join the conda and data portions of the command
            val command = commandPrefix + dataCommand
            val builder = ProcessBuilder(command)

            val redirectError = "$outputDir/export_gvcf_error_${region.replace(":","_")}.log"  // Log file per command
            val redirectOutput = "$outputDir/export_gvcf_output_${region.replace(":","_")}.log"
            builder.redirectOutput(File(redirectOutput))
            builder.redirectError(File(redirectError))

            // Log the command
            //myLogger.info("ExportVcf Command for region: $region")
            myLogger.info("Command: " + builder.command().joinToString(" "))

            // Start the process
            val process = builder.start()
            val error = process.waitFor()

            // Handle error if the process fails
            if (error != 0) {
                myLogger.error("tiledbvcf export for region: ${region}} returned error code $error")
                throw IllegalStateException("Error running tiledbvcf export of dataset $dbPath/gvcf_dataset for batch: ${region}. error: $error")
            }
        }
    }

    // Take the regions file, which should be bed file format, and create a list
    // of regions to send to tiledbvcf export.  Each entry should be of the form
    // "chr:start-end"
    fun convertRegionFileToRegionList(regionsFile: String): List<String> {

        val regions = mutableListOf<String>()
        if (regionsFile.isNotBlank()) {
            val file = File(regionsFile)
            if (!file.exists()) {
                throw IllegalArgumentException("Regions file $regionsFile does not exist")
            }

            // Assume bed file format with CHROM, START, END as the first 3 columns
            // additional columns will be ignored
            file.forEachLine { line ->
                // if the line is empty or is a comment line, skip it
                if (line.isBlank() || line.startsWith("#")) {
                    return@forEachLine
                }

                val fields = line.split("\t")
                if (fields.size < 3) {
                    throw IllegalArgumentException("Format error on line $line : Bed file $regionsFile must have at least 3 fields per line indicating CHROM, START, END")
                }
                // Add 1 to the start position to make it 1-based
                regions.add("${fields[0]}:${fields[1].toInt() + 1}-${fields[2]}")
            }

        }
        return regions
    }
}