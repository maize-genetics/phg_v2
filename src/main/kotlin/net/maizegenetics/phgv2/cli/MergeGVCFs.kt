package net.maizegenetics.phgv2.cli

import biokotlin.util.MergeGVCFUtils.mergeGVCFs
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import org.apache.logging.log4j.LogManager
import java.io.File

class MergeGVCFs : CliktCommand(help = "Merge GVCF files into Single VCF file") {

    private val myLogger = LogManager.getLogger(MergeGVCFs::class.java)

    val inputDir by option(help = "Full path to input GVCF file directory")
        .required()

    val bedFile by option(help = "Path to the bed file used to subset the output VCF files. ")

    val outputFile by option(help = "Full path to output VCF file")
        .required()

    override fun run() {

        logCommand(this)

        // Checks to ensure that the input directory exists
        require(File(inputDir).isDirectory) { "Input GVCF directory does not exist: $inputDir" }

        // Checks to ensure that the output file does not exist
        require(!File(outputFile).isFile) { "Output file already exists: $outputFile" }

        mergeGVCFs(inputDir, outputFile, bedFile)

    }

}