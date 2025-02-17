package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.logCommand
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
        logCommand(this)
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


    /**
     * Function to convert the read mappings into hapId -> count map so we could compare the target hapId to the actual counts
     */
    fun convertReadMappingToHapIdCounts(readMappings: Map<List<String>,Int>) : Map<String,Int> {
        //convert the hapIds into counts
        return readMappings.flatMap { (hapIds,count) ->
            hapIds.map { hapId -> Pair(hapId, count) }
        }
            .groupBy({ it.first }, { it.second })
            .mapValues { it.value.sum() }
    }

    /**
     * Function to create an output file name for exporting the counts
     */
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
                val outputString = buildOutputStringForHapIdsInRefRange(rangeToHapId, range, hapIdToSampleGamete, hapIdCounts, targetSampleName)
                writer.write(outputString)
            }
        }
    }

    fun buildOutputStringForHapIdsInRefRange(
        rangeToHapId: Map<ReferenceRange, List<String>>,
        range: ReferenceRange,
        hapIdToSampleGamete: Map<String, List<SampleGamete>>,
        hapIdCounts: Map<String, Int>,
        targetSampleName: String
    ) : String {
        val hapIds = rangeToHapId[range]!!
        val targetHapId = findTargetHapIdForSample(hapIds, hapIdToSampleGamete, targetSampleName)

        val nonTargetCounts = computeCountsForNonTargetHapids(hapIds, targetHapId, hapIdCounts)

        val highestAlt = if (nonTargetCounts.isEmpty()) 0 else nonTargetCounts.first().first

        val targetCount = hapIdCounts.getOrDefault(targetHapId, 0)
        return "$range\t$targetHapId\t${targetCount}\t${highestAlt}\t${targetCount - highestAlt}\t" +
                "${nonTargetCounts.joinToString(", ") { "${it.first}_${it.second}" }}\n"

    }

    fun computeCountsForNonTargetHapids(
        hapIds: List<String>,
        targetHapId: String,
        hapIdCounts: Map<String, Int>
    ) :List<Pair<Int,String>> {
        return hapIds.filter { it != targetHapId }
            .toSet()
            .map { hapId -> Pair(hapId, hapIdCounts.getOrDefault(hapId, 0)) }
            .sortedByDescending { it.second }
            .map { Pair(it.second, it.first) }

    }

    fun findTargetHapIdForSample(
        hapIds: List<String>,
        hapIdToSampleGamete: Map<String, List<SampleGamete>>,
        targetSampleName: String
    ) : String {
        return hapIds.find { hapId -> hapIdToSampleGamete[hapId]!!.map { sg -> sg.name }.contains(targetSampleName) }?:""
    }
}