package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.api.HaplotypeGraph
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


    val binSize by option(help = "Size of the bins to use for the Position Support Gamete File")
        .int()
        .default(256)

    val hvcfDir by option(help = "Directory containing the hvcf files")
        .required()


    override fun run() {
        logCommand(this)

        val graph = HaplotypeGraph(hvcfDir)

        convertReadMappingFile(readMappingFile, outputDir, binSize, graph)

    }

    fun convertReadMappingFile(readMappingFile: String, outputDir: String, binSize: Int, graph: HaplotypeGraph) {

        val hapIdToPosMap = buildHapIdToPosMap(graph)

        val (header,readMappings) = readInReadMappingFile(readMappingFile)

        val cliCommand = headerCommand(this)

        val ps4GData = convertReadMappingDataToPS4G(readMappings, hapIdToPosMap, binSize)

//        writeOutPS4GFile(ps4GData, outputDir, cliCommand)

    }

    fun buildHapIdToPosMap(graph: HaplotypeGraph) : Map<String,List<Position>> {
        return graph.hapIdToRefRangeMap()
            .map { Pair(it.key, it.value.map { refRange -> Position(refRange.contig,refRange.start) }) }.toMap()
    }

    fun readInReadMappingFile(readMappingFile: String) : Pair<List<String>, Map<Array<String>,Int>> {
        val lines = bufferedReader(readMappingFile).readLines()
        val headerLines = mutableListOf<String>()
        val readMappingCounts = mutableMapOf<Array<String>, Int>()
        for(line in lines) {
            if(line.startsWith("#")) {
                headerLines.add(line)
            } else if(!line.startsWith("HapIds")) {
                val splitLine = line.split("\t")
                readMappingCounts[splitLine.toTypedArray()] = splitLine[1].toInt()
            }
        }

        return Pair(headerLines, readMappingCounts)
    }

    fun convertReadMappingDataToPS4G(readMappings: Map<Array<String>,Int>, hapIdToPosMap: Map<String,List<Position>>, binSize: Int) : Map<String,Map<Position,Int>> {
        TODO()

        val ps4GData = mutableMapOf<String,MutableMap<Position,Int>>()
        for((hapIds,count) in readMappings) {
            for(hapId in hapIds) {
                val posList = hapIdToPosMap[hapId]?: listOf()
                for(pos in posList) {
                    val binPos = Position(pos.contig, pos.position / binSize)
                    val posEncoded = encodePosition(binPos)
                    val posMap = ps4GData.getOrPut(hapId) { mutableMapOf() }
                    posMap[binPos] = posMap.getOrDefault(binPos,0) + count
                }
            }
        }
        return ps4GData
    }

    fun encodePosition(pos: Position) : String {
        //Pack into an Int
        return "${pos.contig}:${pos.position}"
    }
}