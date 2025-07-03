package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.convert
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.utils.verifyURI
import net.maizegenetics.phgv2.utils.writeBedfileFromVcf
import org.apache.logging.log4j.LogManager
import java.io.File
import java.nio.file.Files
import kotlin.io.path.absolutePathString


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
    data class SampleFile(val sampleFile: String) : SampleProcessing() {
        @Override
        override fun getExportCommand(): List<String> {
            check(File(sampleFile).exists()) { "Samples file $sampleFile does not exist." }
            return listOf(SampleFormatEnum.FILE.toString(), sampleFile)
        }
    }

    data class SampleList(val sampleList: String) : SampleProcessing() {
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

    val datasetType by option(help = "Type of dataset to export: choices are gvcf or hvcf")
        .default("hvcf")


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

    val condaEnvPrefix by option (help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
        .default("")

    override fun run() {

        logCommand(this)

        // if using a regions file, the output vcfs can contain duplicate sequential reference blocks which need to be deleted
        // in that case write the vcf files to a temp directory, then write the de-duped vcfs to the output.
        // If a regions-file is specified, check for its existence
        val workingOutputDirectory = if (regionsFile.isNotBlank()) {
            require(File(regionsFile).exists()) {"$regionsFile does not exist."}
            val tmpDir = Files.createTempDirectory("vcfOut").toFile()
            tmpDir.deleteOnExit()
            tmpDir
        } else File(outputDir)


        val dbPath = if (dbPath.isBlank()) {
            System.getProperty("user.dir")
        } else {
            dbPath
        }

        // Verify the tiledbURI - an exception is thrown from verifyURI if the URI is not valid
        verifyURI(dbPath, "hvcf_dataset",condaEnvPrefix)

        // This is the tiledbvcf command we want to run:
        // Doing this with a ProcessBuilder and using the phg_v2 conda environment
        // tiledbvcf export --uri tiledb/hvcf_dataset -O v --sample-names Ref --output-dir exported-vcfs
        // or ... if the user provides a file with sample names
        // tiledbvcf export --uri tiledb/hvcf_dataset -O v --samples-file sampleNames.txt --output-dir exported-vcfs

        val dtype = if (datasetType == "gvcf") "gvcf_dataset" else "hvcf_dataset"

        // Setup the conda environment portion of the command
        var command = if (condaEnvPrefix.isNotBlank()) mutableListOf("conda","run","-p",condaEnvPrefix) else mutableListOf("conda","run","-n","phgv2-tiledb")

        // Tiledbvcf can take either a file with samplenames, or a comma-separated list of sample names
        // setup the command based on user input type.
        var dataCommand = if (samples.getExportCommand()[0] == SampleFormatEnum.FILE.toString()) mutableListOf(
            //"conda",
            //"run",
            //"-n",
            //"phgv2-tiledb",
            "tiledbvcf",
            "export",
            "--uri",
            "$dbPath/$dtype",
            "-O",
            "v",
            "--samples-file",
            samples.getExportCommand()[1],
            "--output-dir",
            workingOutputDirectory.absolutePath
        ) else mutableListOf(
            //"conda",
            //"run",
            //"-n",
            //"phgv2-tiledb",
            "tiledbvcf",
            "export",
            "--uri",
            "$dbPath/$dtype",
            "-O",
            "v",
            "--sample-names",
            samples.getExportCommand()[1],
            "--output-dir",
            workingOutputDirectory.absolutePath
        )

        if (regionsFile.isNotBlank()) {
            when (File(regionsFile).extension) {
                "vcf" -> {
                    val tmpFilename = Files.createTempFile("tempBedfile", ".bed").absolutePathString()
                    writeBedfileFromVcf(regionsFile, tmpFilename)
                    dataCommand.add("--regions-file")
                    dataCommand.add(tmpFilename)
                }

                "bed" -> {
                    dataCommand.add("--regions-file")
                    dataCommand.add(regionsFile)
                }

                else -> throw IllegalArgumentException("Regions file $regionsFile must end in .bed or .vcf.")
            }
        }

        // join the 2 conda and data portion of the commands
        command.addAll(dataCommand)
        val builder = ProcessBuilder(command)

        //val redirectError = "$outputDir/export_${dtype}_error.log"
        //val redirectOutput = "$outputDir/export_${dtype}_output.log"

        val listDirCmd = listOf("conda", "run", "-n", "phgv2-tiledb", "ls", "-alh", outputDir)
        ProcessBuilder(listDirCmd).start().waitFor()


        try {
            //builder.redirectOutput(File(redirectOutput))
            //builder.redirectError(File(redirectError))
        } catch (e: Exception) {
            myLogger.error("Error setting up ProcessBuilder for tiledbvcf export command: ${e.message}")
            throw IllegalStateException("Error setting up ProcessBuilder redirect: ${e.message}", e)
        }

        val exportCommand = builder.command().joinToString(" ")
        myLogger.info("ExportVcf Command: " + builder.command().joinToString(" "))
        val process = try {
            // workingOutputDirectory.mkdirs()
            builder.start()
        } catch (e: Exception) {
            myLogger.error("Error running tiledbvcf export command: ${e.message}")
            throw IllegalStateException("Error running tiledbvcf export command: ${e.message}", e)
        }
        val stderr = process.errorStream.bufferedReader().readText()
        val error = process.waitFor()
        if (error != 0) {
            myLogger.error("tiledbvcf export for: $samples run via ProcessBuilder returned error code $error")
            throw IllegalStateException("Error running tiledbvcf export of dataset $dbPath/$dtype for: $samples. exportCommand: $exportCommand error: $error stderr: $stderr")
        }

        if (regionsFile.isNotBlank()) {
            //get rid of duplicate reference blocks and write the resulting files to the output directory
            val fileList = workingOutputDirectory.listFiles().filter { it.name.endsWith(".vcf") }
            val finalOutputDirectory = File(outputDir)
            for (tmpFile in fileList) {
                val outputFile = finalOutputDirectory.resolve(tmpFile.name)
                deleteDuplicateSequentialRefBlocks(tmpFile, outputFile)
            }

            //finished with the temporary workingOutputDirectory so delete it and contents
            workingOutputDirectory.deleteRecursively()
        }

        val typeSamples = samples.getExportCommand()
        val type = typeSamples[0]
        val sampleNames = typeSamples[1]

        when (type) {
            SampleFormatEnum.FILE.toString() -> {
                File(sampleNames).readLines().forEach { sample ->
                    when (datasetType) {
                        "gvcf" -> File("$outputDir/$sample.vcf").renameTo(File("$outputDir/${sample}.g.vcf"))
                        "hvcf" -> File("$outputDir/$sample.vcf").renameTo(File("$outputDir/${sample}.h.vcf"))
                    }
                }
            }

            SampleFormatEnum.LIST.toString() -> {
                sampleNames.split(",").forEach { sample ->
                    when (datasetType) {
                        "gvcf" -> File("$outputDir/$sample.vcf").renameTo(File("$outputDir/${sample}.g.vcf"))
                        "hvcf" -> File("$outputDir/$sample.vcf").renameTo(File("$outputDir/${sample}.h.vcf"))
                    }
                }
            }
        }

    }

    private fun deleteDuplicateSequentialRefBlocks(inputVcf: File, outputVcf: File) {
        VCFFileReader(inputVcf, false).use { vcfReader ->
            val writer = VariantContextWriterBuilder()
                .unsetOption(Options.INDEX_ON_THE_FLY)
                .setOutputFile(outputVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .build()

            writer.writeHeader(vcfReader.header)

            var previousContext: VariantContext? = null
            for (varctxt in vcfReader) {
                if (areVariantContextsDifferent(previousContext, varctxt)) {
                    writer.add(varctxt)
                    previousContext = varctxt
                }
            }

            writer.close()

        }
    }

    private fun areVariantContextsDifferent(vc1: VariantContext?, vc2: VariantContext?): Boolean {
        return when {
            vc1 == null -> true
            vc2 == null -> true
            vc1.contig != vc2.contig -> true
            vc1.start != vc2.start -> true
            vc2.end != vc2.end -> true
            else -> false
        }
    }
}