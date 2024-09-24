package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import net.maizegenetics.phgv2.utils.getBufferedWriter
import net.maizegenetics.phgv2.utils.retrieveAgcData
import net.maizegenetics.phgv2.utils.verifyURI
import org.apache.logging.log4j.LogManager

/**
 * This class is used to list the sample names in both the AGC compressed file and the tiledb datasets.
 * It is a command that allows users to see what data they have loaded, and make decisions on
 * which sample data to include when imputing, analysing or exporting data.
 *
 * It takes as input a path to the tiledb datasets and an output file.  Data on sample names
 * is written to a tab-delimited file with columns for the sample name, and whether the sample
 * is in the AGC file, the gvcf dataset and/or the hvcf dataset.
 */
class ListSamples: CliktCommand(help = "List the sample names in both the AGC compressed file and the tiledb datasets.") {
    private val myLogger = LogManager.getLogger(ListSamples::class.java)

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        .default("")

    val condaEnvPrefix by option (help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
        .default("")

    val outputFile by option(help = "File where output data will be written")
        .required()

    override fun run() {
        val tileDb = if (dbPath.isBlank()) System.getProperty("user.dir") else dbPath

        // Verify the tiledbURI - an exception is thrown from verifyURI if the URI is not valid
        verifyURI(tileDb, "hvcf_dataset",condaEnvPrefix)

        myLogger.info("URI verified: $tileDb , calling createSampleLists")
        // call command to get list of samples from the AGC file and both tiledb datasets
        createSampleLists(tileDb, condaEnvPrefix, outputFile)
    }

    // Function to create the output file with the list of samples in the AGC file and the TileDB datasets
    // There may be some samples in the hvcf that are not in the gvcf, in particular the reference genome
    // We would expect that otherwise the samples in the gvcf and hvcf datasets would be the same but it
    // is possible users have only loaded a portion of each type vcf.
    fun createSampleLists(dbPath: String, condaEnvPrefix: String, outputFile: String) {
        // call command to get list of samples from the AGC file
        myLogger.info("Getting agc sample list")
        val agcSamples = retrieveAgcData(dbPath,listOf("listset"),condaEnvPrefix).sorted()

        // call command to get list of samples from the TileDB datasets, first from hvcf_dataset
        var uri = dbPath + "/gvcf_dataset"
        myLogger.info("Getting gvcf sample list")
        val gvcfSamples = LoadVcf().getTileDBSampleLists(uri).sorted()

        uri = dbPath + "/hvcf_dataset"
        myLogger.info("Getting hvcf sample list")
        val hvcfSamples = LoadVcf().getTileDBSampleLists(uri).sorted()

        // take the agcSamples, gvcfSamples and hvcfSamples and create a union of the three lists
        val allSamples = agcSamples.union(gvcfSamples).union(hvcfSamples).sorted()
        val headerLine = "SampleName\tInAGC\tInGVCF\tInHVCF\n"

        // write the sample list to the output file.
        getBufferedWriter(outputFile).use { writer ->
            writer.write(headerLine)
            for (sample in allSamples) {
                val inAGC = if (agcSamples.contains(sample)) "Y" else "N"
                val inGVCF = if (gvcfSamples.contains(sample)) "Y" else "N"
                val inHVCF = if (hvcfSamples.contains(sample)) "Y" else "N"
                writer.write("$sample\t$inAGC\t$inGVCF\t$inHVCF\n")
            }
        }

        myLogger.info("Sample list data written to $outputFile")
    }
}