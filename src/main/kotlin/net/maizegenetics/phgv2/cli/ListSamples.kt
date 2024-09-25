package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.convert
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.enum
import net.maizegenetics.phgv2.utils.getBufferedWriter
import net.maizegenetics.phgv2.utils.retrieveAgcData
import net.maizegenetics.phgv2.utils.verifyURI
import org.apache.logging.log4j.LogManager
import kotlin.test.DefaultAsserter.fail

/**
 * This class is used to list the sample names in both the AGC compressed file and the tiledb datasets.
 * It is a command that allows users to see what data they have loaded, and make decisions on
 * which sample data to include when imputing, analysing or exporting data.
 *
 * It takes as input a path to the tiledb datasets, a dataset type for output,  and an output file.
 * By default, the output will be a file, one sample per line, showing
 * the sample names that exist in the hvcf dataset.
 *
 * The optional parameter --data-set, allows the user to specify the dataset from which to pull sample names.  The default is hvcf.
 * The options are agc, gvcf, hvcf, and all.  If all is selected, the samples from all datasets are listed in a tab-delimited file
 * with columns for the sample name, and whether the sample is in the AGC file, the gvcf dataset and/or the hvcf dataset.
 *
 * If "hvcf", "gvcf" or "agc" is selected, the output will be a file with a list of sample names in the selected dataset, one per line
 *
 * Another optional parameter is the conda environment prefix.  If provided, this should be the full path to the conda environment.
 *
 */

// Create an enum class to hold the dataset options
enum class DatasetOptions {
    agc, gvcf, hvcf, all
}
class ListSamples: CliktCommand(help = "List the sample names in both the AGC compressed file and the tiledb datasets.") {
    private val myLogger = LogManager.getLogger(ListSamples::class.java)

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        .default("")

    // Datset option must be one from the enum, default to hvcf
    val dataSet by option(help = "Name of the dataset from which to list the samples. Options are: all, agc, gvcf, hvcf.")
        .enum<DatasetOptions>()
        .default(DatasetOptions.hvcf)

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
        createSampleLists(tileDb, condaEnvPrefix, outputFile, dataSet)
    }

    // Function to create the output file with the list of samples in the AGC file and the TileDB datasets
    // There may be some samples in the hvcf that are not in the gvcf, in particular the reference genome
    // We would expect that otherwise the samples in the gvcf and hvcf datasets would be the same but it
    // is possible users have only loaded a portion of each type vcf.
    fun createSampleLists(dbPath: String, condaEnvPrefix: String, outputFile: String, dataset:DatasetOptions) {


        // call command to get list of samples from the AGC file
        myLogger.info("Getting agc sample list")
        val agcSamples = if (dataset == DatasetOptions.agc || dataset == DatasetOptions.all) retrieveAgcData(dbPath,listOf("listset"),condaEnvPrefix).sorted()
            else emptyList()

        // call command to get list of samples from the TileDB datasets, first from hvcf_dataset
        var uri = dbPath + "/gvcf_dataset"
        myLogger.info("Getting gvcf sample list")
        val gvcfSamples = if (dataset == DatasetOptions.gvcf || dataset == DatasetOptions.all) LoadVcf().getTileDBSampleLists(uri).sorted() else emptyList()

        uri = dbPath + "/hvcf_dataset"
        myLogger.info("Getting hvcf sample list")
        val hvcfSamples = if (dataset == DatasetOptions.hvcf || dataset == DatasetOptions.all) LoadVcf().getTileDBSampleLists(uri).sorted() else emptyList()

        // take the agcSamples, gvcfSamples and hvcfSamples and create a union of the three lists
        val allSamples = agcSamples.union(gvcfSamples).union(hvcfSamples).sorted()
        val headerLine = "SampleName\tInAGC\tInGVCF\tInHVCF\n"

        // Create outputFile based on the dataset option
        when (dataset) {
            DatasetOptions.agc -> {
                // write the values from the agcSamples list to the output file, one per line.
                getBufferedWriter(outputFile).use { writer ->
                    for (sample in agcSamples) {
                        writer.write("$sample\n")
                    }
                }
            }
            DatasetOptions.gvcf -> {
                // write the values from the gvcfSamples list to the output file, one per line.
                getBufferedWriter(outputFile).use { writer ->
                    for (sample in gvcfSamples) {
                        writer.write("$sample\n")
                    }
                }
            }
            DatasetOptions.hvcf ->  {
                // write the values from the hvcfSamples list to the output file, one per line.
                getBufferedWriter(outputFile).use { writer ->
                    for (sample in hvcfSamples) {
                        writer.write("$sample\n")
                    }
                }
            }
            DatasetOptions.all -> {
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
            }
        }
        myLogger.info("Sample list data for ${dataset} written to $outputFile")
    }
}