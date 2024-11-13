package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import net.maizegenetics.phgv2.api.HaplotypeGraph
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * This class is a pipeline that runs a series of commands to create tables for the imputation results.
 *
 * The first table created uses h.vcf files from aligned assembly processing to create a table of haplotype IDs to SampleGamete.
 * The second table created uses h.vcf files from the imputation pipeline to create a table of seed samples by reference range.
 * The third function attempts to print a vcf of variants for each reference range:  This will create as many vcf files
 * as there are reference ranges.  The vcf files will be named based on the reference range contig:coordinates
 *
 */
class ImputationResultsTables: CliktCommand(help = "Metrics for WeiYun's project") {
    private val myLogger = LogManager.getLogger(ImputationResultsTables::class.java)
    // THis is the assembly hvcf dir - it should be renamed if we need as well an imputation hvcf dir
    val hvcfDir by option("--hvcf-dir", help = "Path to directory holding assembly alignment created hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        .required()

    val outputDir by option (help = "Output directory for the gVCF files.  If not provided, the current working directory is used.")
        .default("")

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
        // Table1: HapId to Gamete:
        // THis table expects the graph to be built from assembly alignment HVCF files
        HapidSampleTable().createHapidGameteTable(graph, outputDir)

        // Table2: SeedSample by RefRange: Value is HapID:  call SampleHapidByRange
        // This table expects the graph to be built from imputation pipeline HVCF files

        // Table3 (subset by RefRange i.e. 67K of these): RefPosition by Gamete: Value is GT


    }



}