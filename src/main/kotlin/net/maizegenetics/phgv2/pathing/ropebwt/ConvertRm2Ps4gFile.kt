package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.headerCommand
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.utils.Position
import org.apache.logging.log4j.LogManager


/**
 * Command to convert the Read Mapping file to a Position Support Gamete File
 *
 * Read Mapping file:
 * hapIdSet\tcount
 *
 *
 * Position Support 4 Gamete File:
 * gameteSet\tpos\tcount
 *
 * The Position support gamete File will keep track of the number of reads/snps that support a given gamete at a binned position.
 */
class ConvertRm2Ps4gFile : CliktCommand(help = "Convert Read Mapping file to Position Support Gamete File") {

    private val myLogger = LogManager.getLogger(ConvertRm2Ps4gFile::class.java)

    //TODO Add option for keyFile
    val readMappingFile by option(help = "Read Mapping file")
        .required()

    val outputDir by option(help = "Output directory")
        .required()

    val hvcfDir by option(help = "Directory containing the hvcf files")
        .required()

    val sortPositions by option(help = "Sort positions in the resulting PS4G file.")
        .flag(default = true)


    override fun run() {
        logCommand(this)
        val cliCommand = headerCommand(this)

        myLogger.info("Building Graph")
        val graph = HaplotypeGraph(hvcfDir)
        myLogger.info("Converting read mapping file")
        convertReadMappingFile(readMappingFile, outputDir, graph, cliCommand, sortPositions)

    }

    /**
     * Function to convert the read mapping file into a PS4G file
     */
    fun convertReadMappingFile(readMappingFile: String, outputDir: String, graph: HaplotypeGraph, cliCommand: String, sortPositions: Boolean = true) {

        myLogger.info("Loading in readMapping File: $readMappingFile")
        val (header,readMappings) = readInReadMappingFile(readMappingFile)


        myLogger.info("Converting readMappings to PS4GData")
        val (ps4GData, sampleGameteCount, gameteToIdxMap) = convertReadMappingDataToPS4G(readMappings, graph, sortPositions)


        val outputFile = PS4GUtils.buildOutputFileName(readMappingFile, outputDir)
        myLogger.info("Writing out PS4G file to $outputFile")
        PS4GUtils.writeOutPS4GFile(ps4GData, sampleGameteCount, gameteToIdxMap, outputFile, header, cliCommand)

    }

    /**
     * Function to read in the readMapping file and retain the header.
     */
    fun readInReadMappingFile(readMappingFile: String) : Pair<List<String>, Map<List<String>,Int>> {
        val lines = bufferedReader(readMappingFile).readLines().filter { it.isNotBlank() }
        val headerLines = mutableListOf<String>()
        val readMappingCounts = mutableMapOf<List<String>, Int>()
        for(line in lines) {
            if(line.startsWith("#")) {
                headerLines.add(line)
            } else if(!line.startsWith("HapIds")) {
                val splitLine = line.split("\t")
                val hapIds = splitLine[0].split(",")
                readMappingCounts[hapIds] = splitLine[1].toInt()
            }
        }

        return Pair(headerLines, readMappingCounts)
    }

    /**
     * Function to create the PS4G data information for the readMappings
     */
    fun convertReadMappingDataToPS4G(readMappings: Map<List<String>,Int>,
                                     graph: HaplotypeGraph,
                                     sortPositions: Boolean = true) : Triple<List<PS4GData>, Map<SampleGamete,Int>,Map<SampleGamete,Int>> {

        val gameteToIdxMap = graph.sampleGametesInGraph().mapIndexed { index, s -> Pair(s,index) }.toMap()

        val hapIdToRanges = graph.hapIdToRefRangeMap()

        val gameteCountMap = mutableMapOf<SampleGamete,Int>()

        val ps4GData = buildPS4GData(readMappings, hapIdToRanges, graph, gameteCountMap, gameteToIdxMap, sortPositions)

        return Triple(ps4GData, gameteCountMap, gameteToIdxMap)
    }

    /**
     * Function to build the PS4GData for all the readMappings with optional sorting
     */
    fun buildPS4GData(
        readMappings: Map<List<String>, Int>,
        hapIdToRanges: Map<String, List<ReferenceRange>>,
        graph: HaplotypeGraph,
        gameteCountMap: MutableMap<SampleGamete, Int>,
        gameteToIdxMap: Map<SampleGamete, Int>,
        sortPositions: Boolean = true
    ) : List<PS4GData> {

        val readMappingsForConversion = if(sortPositions) {
            readMappings.keys.sortedBy { it.flatMap { hapId -> hapIdToRanges[hapId]!! }.maxBy { hapId -> hapId } }
        } else {
            readMappings.keys
        }

        return readMappingsForConversion.map { hapIdSet ->
            val count = readMappings[hapIdSet]!!
            createPS4GFileForSingleMapping(hapIdSet, hapIdToRanges, graph, gameteCountMap, count, gameteToIdxMap)
        }
    }

    /**
     * Function to create the PS4GData for a single mapping
     */
    fun createPS4GFileForSingleMapping(
        hapIdSet: List<String>,
        hapIdToRanges: Map<String, List<ReferenceRange>>,
        graph: HaplotypeGraph,
        gameteCountMap: MutableMap<SampleGamete, Int>,
        count: Int,
        gameteToIdxMap: Map<SampleGamete, Int>
    ): PS4GData {
        val mostHitRefRange = hapIdSet.flatMap { hapId -> hapIdToRanges[hapId] ?: listOf() }.maxOf { it }
        //build binned position:
        val pos = Position(mostHitRefRange.contig, mostHitRefRange.start/256) //Need to bin here otherwise it won't be consistent

        //convert the hapIdSet into gametes
        val hapIdForSampleGameteMap = graph.hapIdToSampleGametes(mostHitRefRange)

        val sampleGameteIdxSorted = hapIdSet.flatMap { hapId -> hapIdForSampleGameteMap[hapId] ?: listOf() }
            .mapNotNull { sampleGamete ->
                gameteCountMap[sampleGamete] = gameteCountMap.getOrDefault(sampleGamete, 0) + count
                gameteToIdxMap[sampleGamete]
            }
            .sorted()

        return PS4GData(sampleGameteIdxSorted, pos, count)
    }

}