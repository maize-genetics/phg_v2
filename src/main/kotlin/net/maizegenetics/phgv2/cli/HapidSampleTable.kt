package net.maizegenetics.phgv2.cli

import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import net.maizegenetics.phgv2.api.HaplotypeGraph
import org.apache.logging.log4j.LogManager
import java.io.File

class HapidSampleTable : CliktCommand(help = "Create a table of haplotype IDs to SampleGamete.  Can be multiple samples mapping to each hapid"){

    private val myLogger = LogManager.getLogger(HapidSampleTable::class.java)
    val hvcfDir by option("--hvcf-dir", help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        .required()

    val outputDir by option (help = "Output directory for the gVCF files.  If not provided, the current working directory is used.")
        .default("")

    override fun run() {
        // Get list of files, create graph, write table
        val hvcfFiles = File(hvcfDir).listFiles { file -> file.name.endsWith(".h.vcf") || file.name.endsWith(".h.vcf.gz") }.map { it.path }
        val startTime = System.nanoTime()
        val graph = HaplotypeGraph(hvcfFiles)
        val timeToBuildGraph = (System.nanoTime() - startTime) / 1e9
        myLogger.info("Time to build graph: $timeToBuildGraph")

        // Table1: HapId to Gamete
        createHapidGameteTable(graph, outputDir)
    }

    // Function to create a tab-delimited file with haplotype id and gametes.
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