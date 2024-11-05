package net.maizegenetics.phgv2.cli

import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.pathing.BuildKmerIndex
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * This class will take imputed hvcf files and create tables that can be used for imputation analysis
 */
class ImputationResultsTables: CliktCommand(help = "Metrics for WeiYun's project") {
    private val myLogger = LogManager.getLogger(ImputationResultsTables::class.java)
    val hvcfDir by option("--hvcf-dir", help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        .required()

    val outputDir by option (help = "Output directory for the gVCF files.  If not provided, the current working directory is used.")
        .default("")

    // Not sure we need these
//    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
//        .default("")
//    val condaEnvPrefix by option (help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
//        .default("")

    // Do we need the reference?  I'm thinking not
//    val referenceFile by option(help = "Path to local Reference FASTA file needed for sequence dictionary")
//        .required()
    override fun run() {
        //TODO("Not yet implemented")
        createHvcfMetrics(hvcfDir, outputDir)
    }

    fun createHvcfMetrics(hvcfDir: String,  outputDir: String) {

        // First, build a graph from the hvcf files
        val hvcfFiles = File(hvcfDir).listFiles { file -> file.name.endsWith(".h.vcf") || file.name.endsWith(".h.vcf.gz") }.map { it.path }
        val startTime = System.nanoTime()
        val graph = HaplotypeGraph(hvcfFiles)
        val timeToBuildGraph = (System.nanoTime() - startTime) / 1e9
        myLogger.info("Time to build graph: $timeToBuildGraph")

        // We have the graph, we now want this data:
        // Table1: HapId to Gamete
        createHapidGameteTable(graph, outputDir)

        // Table2: SeedSample by RefRange: Value is HapID
        // Table3 (subset by RefRange i.e. 67K of these): RefPosition by Gamete: Value is GT


    }

    // Function to create a tab-delimited file with haplotype id and gametes
    // gametes (ie sample names) are a comma-separated list in the second column
    fun createHapidGameteTable(graph:HaplotypeGraph, outputDir:String) {

        val time = System.nanoTime()
        val hapidGameteFile = "${outputDir}/hapid_gamete.tsv"
        bufferedWriter(hapidGameteFile).use { writer ->
            graph.ranges().forEachIndexed { rangeIndex, range ->
                graph.hapIdToSampleGametes(range).forEach { (hapId, gametes) ->
                    // create a comma separated string of gametes from gametes
                    // The gametes are SampleGamete objects, which are a combination of a sample name and a gamete id
                    // For each gamete, substring it to remove everything from the last colon on
                    val gameteList = gametes.map { it.toString().substringBeforeLast(":") }
                        .joinToString(",")
                    writer.write("$hapId\t$gameteList\n")
                }
            }
        }
        myLogger.info("Time to write hapid_gamete.tsv: ${(System.nanoTime() - time) / 1e9}")
    }
}