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

class ReadMappingCountQc : CliktCommand("Read mapping count QC") {
    val myLogger = LogManager.getLogger(ReadMappingCountQc::class.java)

    val hvcfDir by option(help = "Directory with Haplotype VCF files")
        .required()

    val readMappingFile by option(help = "Read mapping file")
        .required()

    val targetSampleName by option(help = "Target sample name")
        .required()

    val outputDir by option(help = "Output directory")
        .required()

    override fun run() {
        myLogger.info("Read mapping count QC")
        processReadMappingCounts(hvcfDir, readMappingFile, targetSampleName, outputDir)
    }
    fun processReadMappingCounts(hvcfDir: String, readMappingFile: String, targetSampleName: String, outputDir: String) {
        myLogger.info("Processing read mapping counts")
        val graph = HaplotypeGraph(hvcfDir)
        val hapIdToSampleGamete = graph.hapIdsToSampleGametes()
        val rangeToHapIds = graph.refRangeToHapIdList()

        val outputFile = buildCountOutputFile(outputDir, readMappingFile ,targetSampleName)

        val readMapping = AlignmentUtils.importReadMapping(readMappingFile)

        val hapIdCounts = convertReadMappingToHapIdCounts(readMapping)

        writeOutCounts(hapIdCounts, outputFile, targetSampleName, rangeToHapIds, hapIdToSampleGamete)
    }



    fun convertReadMappingToHapIdCounts(readMappings: Map<List<String>,Int>) : Map<String,Int> {
        //convert the hapIds into counts
        return readMappings.flatMap { (hapIds,count) ->
            hapIds.map { hapId -> Pair(hapId, count) }
        }
            .groupBy({ it.first }, { it.second })
            .mapValues { it.value.sum() }
    }

    fun buildCountOutputFile(outputDir: String, readMappingFile: String, targetSampleName: String) : String {
        //Get just the file name from the Mapping file without the extension
        return "$outputDir/${File(readMappingFile).nameWithoutExtension}_${targetSampleName}_counts.txt"
    }

    fun writeOutCounts(hapIdCounts: Map<String,Int>,
                       outputFile: String,
                       targetSampleName: String,
                       rangeToHapId: Map<ReferenceRange, List<String>>,
                       hapIdToSampleGamete: Map<String, List<SampleGamete>>) {
        myLogger.info("Writing out counts")
        bufferedWriter(outputFile).use { writer ->
            writer.write("refRange\t${targetSampleName}_HapID\t${targetSampleName}_HapCount\tHighestAltCount\tDifference\tOtherHapCounts\n")
            rangeToHapId.keys.sorted().forEach { range ->
                val hapIds = rangeToHapId[range]!!
                val targetHapId = hapIds.find { hapId -> hapIdToSampleGamete[hapId]!!.map { sg -> sg.name }.contains("Zm-B73-REFERENCE-NAM-5.0") }

                val counts = hapIds.filter { it != targetHapId }
                    .map { hapId ->Pair(hapId,hapIdCounts.getOrDefault(hapId, 0)) }
                    .sortedByDescending { it.second }
                    .map { Pair(it.second,it.first) }

                val highestAlt = if(counts.isEmpty()) 0 else counts.first().first

                val targetCount = hapIdCounts.getOrDefault(targetHapId,0)
                writer.write("$range\t$targetHapId\t${targetCount}\t${highestAlt}\t${targetCount - highestAlt}\t" +
                        "${counts.joinToString(", ") { "${it.first}_${it.second}" }}\n")
            }
        }
    }
}