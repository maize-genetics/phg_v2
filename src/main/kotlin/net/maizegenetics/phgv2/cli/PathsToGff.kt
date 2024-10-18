package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.convert
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import net.maizegenetics.phgv2.utils.loadGFFsToGff3Feature
import net.maizegenetics.phgv2.utils.makeGffFromHvcf
import org.apache.logging.log4j.LogManager

/**
 * This class is used to create a GFF file from a PHG imputation Path H.vcf file.
 * it takes as input
 *   1. a keyfile that links sample names to GFF files, (the GFF files
 *      must exist at the path specified in the keyfile)
 *   2. PHG paths file:  an h.vcf file created by the PHG imputation pipeline
 *   3. An optional output file where the pseudo GFF file will be written.
 */

class PathsToGff: CliktCommand( help = "Create a GFF file from a PHG imputation Path H.vcf file") {

    private val myLogger = LogManager.getLogger(PathsToGff::class.java)

    val keyFile by option(help = "Tab-delimited file containing 2 column: SampleName and GFF.  GFF column contains full path name for the GFF file for that sample.  SampleName contains the sample name for that assembly, e.g. B73 or CML247, and must match the name as stored in TileDB and AGC. ")
        .required()

    val hvcfFile by option(help = "Full path to the H.vcf file")
        .required()

    val outputFile by option("-o", "--output-dir", help = "Directory where pseudo GFF file will be written")
        .default("")

    // Pre-compile the Regex pattern - used when creating the output fasta file names
    val HVCF_PATTERN = Regex("""(\.hvcf|\.h\.vcf|\.hvcf\.gz|\.h\.vcf\.gz)$""")

    override fun run() {
        val resultsTreeCenter = loadGFFsToGff3Feature(keyFile)
        myLogger.info("Calling makeGffFromPath for file ${hvcfFile}")
        val time = System.nanoTime()
        // makeGffFromHvcf will return a Set<Gff3Feature> , but will also write a file if an
        // output file is specified.
        // This flexibility is to allow for user to store the results in memory for further
        // analysis from a jupyter notebook or other environment.
        val taxonPathGFF = makeGffFromHvcf(hvcfFile, resultsTreeCenter,outputFile)
        val endingTime = (System.nanoTime() - time)/1e9
        if (outputFile == "") {
            myLogger.info("makeGffFromHvcf for file ${hvcfFile} took ${endingTime} seconds, no output file written")
        } else {
            myLogger.info("makeGffFromHvcf for file ${hvcfFile} took ${endingTime} seconds, file written to ${taxonPathGFF}")
        }

    }
}