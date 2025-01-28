package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.headerCommand
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.utils.Position
import org.apache.logging.log4j.LogManager


data class PS4GData(val gameteList: List<Int>, val pos: Int, val count: Int)
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


        val (header,readMappings) = readInReadMappingFile(readMappingFile)

        val cliCommand = headerCommand(this)


        val (ps4GData, sampleGameteCount, gameteToIdxMap) = convertReadMappingDataToPS4G(readMappings, binSize, graph)

        val outputFile = buildOutputFile(readMappingFile, outputDir)

        writeOutPS4GFile(ps4GData, sampleGameteCount, gameteToIdxMap, outputFile, header, cliCommand)

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

    fun convertReadMappingDataToPS4G(readMappings: Map<Array<String>,Int>,
                                     binSize: Int,graph: HaplotypeGraph ) : Triple<List<PS4GData>, Map<SampleGamete,Int>,Map<SampleGamete,Int>> {

        val contigToIdxMap = graph.contigs.mapIndexed { index, s -> Pair(s,index)  }.toMap()

        val gameteToIdxMap = graph.sampleGametesInGraph().mapIndexed { index, s -> Pair(s,index) }.toMap()

        val hapIdToRanges = graph.hapIdToRefRangeMap()

        val gameteCountMap = mutableMapOf<SampleGamete,Int>()

        val ps4GData = readMappings.map { (hapIdSet,count) ->
            val mostHitRefRange = hapIdSet.flatMap { hapId -> hapIdToRanges[hapId]?: listOf() }.maxOf { it }
            //build binned position:
            val binPos = Position(mostHitRefRange.contig, mostHitRefRange.start / binSize)
            val posEncoded = encodePosition(binPos, contigToIdxMap)

            //convert the hapIdSet into gametes
            val hapIdForSampleGameteMap = graph.hapIdToSampleGametes(mostHitRefRange)

            val sampleGameteIdxSorted = hapIdSet.flatMap { hapId -> hapIdForSampleGameteMap[hapId]?: listOf() }
                .mapNotNull { sampleGamete ->
                    gameteCountMap[sampleGamete] = gameteCountMap.getOrDefault(sampleGamete,0) + count
                    gameteToIdxMap[sampleGamete]
                }
                .sorted()

            PS4GData(sampleGameteIdxSorted, posEncoded, count)
        }

        return Triple(ps4GData, gameteCountMap, gameteToIdxMap)
    }

    fun encodePosition(pos: Position, contigOrdering: Map<String, Int>) : Int {
        //Pack into an Int
        val idx = contigOrdering[pos.contig]?: throw IllegalArgumentException("Contig ${pos.contig} not found in contigOrdering")

        //pack last 4 bits of idx into first 8 bits of output then pack the position minus 8 bits into the last 28 bits
        val idxBits = idx and 0xFF //If there are more than 256 contigs this will have unexpected issues
        val posBits = pos.position/256 // div 256 effectively bitshifts by 8

        return (idxBits shl 28) or posBits //we dont care if its negative as we arent comparing them
    }

    fun buildOutputFile(readMappingFile: String, outputDir: String) : String {
        val fileName = readMappingFile.split("/").last().removeSuffix(".txt")
        return "$outputDir/${fileName}_ps4g.txt"
    }

    fun writeOutPS4GFile(pS4GData: List<PS4GData>, sampleGameteCount: Map<SampleGamete,Int>,gameteToIdxMap: Map<SampleGamete,Int> ,outputFile: String, header: List<String>,cliCommand: String) {
        bufferedWriter(outputFile).use { writer ->
            writer.write("#PS4G\n")
            writer.write("#ReadMappingHeader:\n")
            header.forEach { writer.write("$it\n") }
            writer.write("#Command: $cliCommand\n")
            writer.write("#gamete\tgameteIndex\tcount\n")
            sampleGameteCount.forEach { (sampleGamete, count) ->
                writer.write("#$sampleGamete\t${gameteToIdxMap[sampleGamete]}\t$count\n")
            }
            writer.write("gameteSet\tpos\tcount\n")
            pS4GData.forEach { (gameteList, pos, count) ->
                writer.write("${gameteList.joinToString(",")}\t$pos\t$count\n")
            }

        }
    }
}