package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import org.apache.logging.log4j.LogManager
import java.io.File

class MappingCountTableQc: CliktCommand("Read mapping count QC") {

    val myLogger = LogManager.getLogger(ReadMappingCountQc::class.java)

    val hvcfDir by option(help = "Directory with Haplotype VCF files")
        .required()

    val readMappingDir by option(help = "Read mapping directory")
        .required()

    val outputDir by option(help = "Output directory")
        .required()


    override fun run() {
        myLogger.info("Read mapping count QC")

        buildMappingTables(hvcfDir, File(readMappingDir), File(outputDir))
    }

    fun buildMappingTables(hvcfDir: String, readMappingDir: File, outputDir: File) {
        myLogger.info("Building mapping tables")
        myLogger.info("Building Haplotype Graph")
        val graph = HaplotypeGraph(hvcfDir)

        val hapIdAndRangeToSampleGametes = graph.hapIdAndRefRangeToSampleGametes()
        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()


        //Walk through each readMapping file in ReadMappingDir
        readMappingDir.walkTopDown().filter { it.isFile }
            .filter { it.name.endsWith("_readMapping.txt") }
            .forEach { processSingleReadMappingFile(graph,it, outputDir, hapIdAndRangeToSampleGametes, hapIdToRefRangeMap) }
    }

    fun processSingleReadMappingFile(graph: HaplotypeGraph, readMappingFile: File, outputDir: File,
                                     hapIdAndRangeToSampleGametes: Map<Pair<ReferenceRange,String>,List<SampleGamete>>,
                                     hapIdToRefRangeMap: Map<String,List<ReferenceRange>>) {
        myLogger.info("Processing read mapping file: ${readMappingFile.name}")

        val outputFile = buildCountOutputFile(outputDir.absolutePath, readMappingFile)

        val readMapping = AlignmentUtils.importReadMapping(readMappingFile.absolutePath)

        val hapIdCounts = convertReadMappingToHapIdCounts(readMapping, hapIdAndRangeToSampleGametes, hapIdToRefRangeMap)

        writeOutCounts(hapIdCounts, graph,outputFile)
    }

    fun buildCountOutputFile(outputDir: String, readMappingFile: File): File {
        //Remove the _readMapping.txt and add _mappingCounts.txt
        val outputFileName = readMappingFile.nameWithoutExtension.replace("_readMapping", "_mappingCounts") + ".txt"
        return File(outputDir, outputFileName)
    }

    fun convertReadMappingToHapIdCounts(readMappings: Map<List<String>,Int>,
                                        hapIdAndRangeToSampleGametes: Map<Pair<ReferenceRange,String>,List<SampleGamete>>,
                                        hapIdToRefRangeMap: Map<String,List<ReferenceRange>>): Map<Pair<ReferenceRange, SampleGamete>, Int> {
        val finalCountMap = mutableMapOf<Pair<ReferenceRange, SampleGamete>, Int>()

        //go through each hapIdSet from the readMappings.
        // Determine which reference range it belongs to.
        // Then walk through the hapIds in the set and increase the count for that RefRange+hapId pair
        for((hapIdSet, count) in readMappings) {
            val topRefRangeHit = findRefRangeForHapIdSet(hapIdSet, hapIdToRefRangeMap)
            if(topRefRangeHit.contig == "UNKNOWN") {
                myLogger.warn("Haplotype $hapIdSet with count $count unable to resolve reference range.  " +
                        "Make sure the HVCFs used match what were used for building the ropebwt index")
                continue
            }

            //build the list of sample gametes to increase
            hapIdSet.flatMap { hapIdAndRangeToSampleGametes[Pair(topRefRangeHit, it)]!! }
                .forEach {
                    //increment the count map
                    finalCountMap[Pair(topRefRangeHit, it)] = finalCountMap.getOrDefault(Pair(topRefRangeHit, it), 0) + count
                }

        }
        return finalCountMap
    }

    fun findRefRangeForHapIdSet(hapIdSet: List<String>, hapIdToRefRange: Map<String,List<ReferenceRange>>): ReferenceRange {
        //Simple algorithm is just to walk through and count how many times each refRange is hit.
        // Take the max and check to see if the count equals hapIdSet size
        val refRangeCounts = hapIdSet.flatMap { hapId ->
            if(!hapIdToRefRange.containsKey(hapId)) {
               myLogger.info("Haplotype $hapId not in hapIdToRefRange map")
            }
            hapIdToRefRange[hapId]!!
        }.groupingBy { it }.eachCount()

        //get the max
        val mostHitRefRange = refRangeCounts.maxByOrNull { it.value }!!
        return if(mostHitRefRange.value != hapIdSet.size) {
            ReferenceRange("UNKNOWN", -1,-1)
        } else {
            mostHitRefRange.key
        }
    }

    fun writeOutCounts(hapIdCounts: Map<Pair<ReferenceRange, SampleGamete>, Int>, graph: HaplotypeGraph, outputFile: File) {
        val sampleGametes = graph.sampleGametesInGraph()
        val ranges = graph.ranges().sorted()
        bufferedWriter(outputFile.absolutePath).use { writer ->
            //writing out refRange and all the sampleGametes
            writer.write("RefRange\t${sampleGametes.joinToString("\t")}\n")
            for(range in ranges) {
                val countsForRange = sampleGametes.map { sampleGamete ->
                    hapIdCounts.getOrDefault(Pair(range, sampleGamete), 0)
                }
                writer.write("$range\t${countsForRange.joinToString("\t")}\n")
            }
        }
    }


}