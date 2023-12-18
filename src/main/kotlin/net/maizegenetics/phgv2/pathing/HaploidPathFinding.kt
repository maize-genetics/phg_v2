package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.int
import kotlinx.coroutines.channels.Channel
import net.maizegenetics.phgv2.api.HaplotypeGraph
import java.io.File

/**
 * Finds the most likely path through a HaplotypeGraph based on read mapping data. Input: a keyfile with taxa names
 * and read mapping file paths and a haplotype graph. Read mappings must be in the form of
 * haplotype list -> count of reads mapping to the haplotypes. The haplotype graph must be the one used to map reads or
 * a subset of that restricted to specific taxa and/or ReferenceRanges.
 *
 * Steps:
 * 1. Read the keyfile
 * 1. Build the HaplotypeGraph
 * 2. For each entry in the keyfile (multithreaded)
 *      a. Get the read mapping data
 *      b. Use Viterbi algorithm to find the best path
 *      c. Store the path (write hvcf)
 *
 * Outstanding issues:
 * How will haplotype graphs be reproduced? (A list of taxa and ReferenceRanges is sufficient to reproduce
 * a HaplotypeGraph, though additional parameters may be necessary. That probably needs to be stored in the
 * read mapping file or keyfile).
 *
 * Parameters:
 * Number of threads
 *
 *
 */
class HaploidPathFinding : CliktCommand(help = "Create gVCF and hVCF from Anchorwave MAF files") {
    val pathKeyfile by option(help = "key file with columns: SampleName, ReadMappingFiles")
        .required()
        .validate() { require(File(it).exists()) {"$it is not a valid file"} }

    val hvcfDir by option(help = "The directory containing the hvcf files used to build a HaplotypeGraph for path finding.")
        .required()
        .validate() { require(File(it).isDirectory)}

    val threads by option(help = "number of threads used to find paths.").int().default(3)

    //other parameters: probReadMappedCorrectly, useMostLikelyParents, maxParents, minCoverage, likelyParentFile

    private fun buildHaplotypeGraph(): HaplotypeGraph {
        val listOfHvcfFilenames = File(hvcfDir).listFiles().map { it.path }
        return HaplotypeGraph(listOfHvcfFilenames)
    }

    override fun run() {
        TODO("Not yet implemented")

    }

    private fun processKeyFile() {
        //build HaplotypeGraph
        val keyFile = KeyFile(File(pathKeyfile), listOf("SampleName","ReadMappingFiles"), listOf("SampleName"))
        val sampleNameColIndex = keyFile.columnNameMap["SampleName"]
        val readMappingColIndex = keyFile.columnNameMap["ReadMappingFiles"]

        //TODO migrate Viterbi
        //TODO loop through keyFile records to run Viterbi and store resulting paths

    }

    data class ReadMappingResult(val name: String) //placeholder
    data class Path(val name: String, val hapids: List<String>) //placeholder
    private fun processReadMappings() {
        //create a channel for read mappings
        val readMappingChannel = Channel<ReadMappingResult>(10)

        //create a channel for paths
        val pathChannel = Channel<Path>(10)

        //load read mappings for each sample into a channel
        //create worker threads to process entries from the read mapping channel

        //create a coroutine to store paths

    }

}