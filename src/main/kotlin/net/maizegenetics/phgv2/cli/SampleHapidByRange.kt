package net.maizegenetics.phgv2.cli

import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.SampleGamete
import org.apache.logging.log4j.LogManager
import java.util.*

/**
 * Create a table of haplotype IDs by range.
 * This class takes a directory of HVCF files created from the imputation pipeline and creates a HaplotypeGraph.
 * It then writes a table of haplotype IDs (for each sample) by
 * reference range to a file.
 */
class SampleHapidByRange : CliktCommand(help = "Create a table of haplotype IDs by range") {

    private val myLogger = LogManager.getLogger(SampleHapidByRange::class.java)

    val inputDir by option(help = "Full path to input HVCF file directory")
        .required()

    val outputFile by option(help = "Full path to output table summary file")
        .required()

    override fun run() {

        logCommand(this)

        // Merge the HVCF files into a HaplotypeGraph
        val graph = HaplotypeGraph(inputDir)

        // Write the table of haplotype IDs by range to the output file
        writeTable(graph, outputFile)

    }

    // #CHROM POS                                B73                          SEEDGWAS1                         SEEDGWAS10
    // chr1   315223 <c9ecfe3967a71282f3ad7c41d48e0bbf> <b19364bc9a4c07a80986b1ee181446c2> <5c8e72b2e9f11ecc652d5b8e8d0e5bf3>
    // chr1   328917 <f162e742c4d30f151ae6276fbebe762c> <fdfdaa361c39cf5b6f13fad195d0e519> <283a8261c193212fd5cf43d208673322>
    // chr1   412242 <471d4abbf0545dede647e65915345648> <d6dd5ecea7fb4e6f77f9e630f601b7a8> <13e0ac1a8d12e1aedd6a5302d1e221fd>
    private fun writeTable(graph: HaplotypeGraph, outputFile: String) {

        bufferedWriter(outputFile).use { writer ->

            val samples = graph.samples()

            // Write header
            writer.write("#CHROM\tSTART\tEND")
            samples.forEach { sample ->
                writer.write("\t${sample}")
            }
            writer.write("\n")

            graph.ranges().forEach { range ->
                writer.write("${range.contig}\t${range.start}\t${range.end}")
                val sampleToHapid = graph.sampleGameteToHaplotypeId(range)
                val sampleToGametes = mutableMapOf<String, TreeSet<SampleGamete>>()
                sampleToHapid.keys.forEach { gamete ->
                    sampleToGametes.getOrPut(gamete.name) { TreeSet() }.add(gamete)
                }
                samples.forEach { sample ->
                    val hapids =
                        sampleToGametes[sample]?.joinToString(separator = "/") { "<${sampleToHapid[it]}>" } ?: "."
                    writer.write("\t$hapids")
                }
                writer.write("\n")
            }

        }

    }

}