package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.UsageError
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
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

    val targetSampleName by option(help = "Target sample name. Required unless the --by-sample flag is set.").default("")

    val outputDir by option(help = "Output directory")
        .required()

    val bySample by option(help = "Create a report with a read count for each sample by reference range. " +
            "The default is a report with only one row per reference range. --by-sample report does not use --target-sample-name.").flag()

    override fun run() {
        logCommand(this)
        myLogger.info("Read mapping count QC")
        if (!bySample && targetSampleName.isNullOrBlank()) {
            throw UsageError("--target-sample-name is required unless --by-sample is set")
        }
        processReadMappingCounts(hvcfDir, readMappingFile, targetSampleName, outputDir, bySample)
    }

    fun processReadMappingCounts(hvcfDir: String, readMappingFile: String, targetSampleName: String, outputDir: String, reportBySample: Boolean) {
        myLogger.info("Processing read mapping counts")
        val graph = HaplotypeGraph(hvcfDir)
        val hapIdToSampleGamete = graph.hapIdsToSampleGametes()
        val rangeToHapIds = graph.refRangeToHapIdList()

        val readMapping = AlignmentUtils.importReadMapping(readMappingFile)

        val hapIdCounts = convertReadMappingToHapIdCounts(readMapping)

        val referenceRangeCounts = getReadMappingCountsByReferenceRange(graph, readMapping)

        if (reportBySample) {
            val bysampleOutFile = "$outputDir/${File(readMappingFile).nameWithoutExtension}_counts_bysample.txt"
            writeDetailReport(hapIdCounts, bysampleOutFile, rangeToHapIds, hapIdToSampleGamete, referenceRangeCounts)
        } else {
            val outputFile = buildCountOutputFile(outputDir, readMappingFile ,targetSampleName)
            writeOutCounts(hapIdCounts, outputFile, targetSampleName, rangeToHapIds, hapIdToSampleGamete, referenceRangeCounts)
        }

    }

    /**
    * Function that writes a report with a separate row for each SampleName and ReferenceRange.
    */
    fun writeDetailReport(hapIdCounts: Map<String, Int>,
                          outputFile: String,
                          rangeToHapIds: Map<ReferenceRange,List<String>>,
                          hapIdToSampleGamete: Map<String, List<SampleGamete>>,
                          referenceRangeCounts: Map<ReferenceRange, Int>) {

        val referenceRangeList = referenceRangeCounts.keys.toList().sorted()

        bufferedWriter(outputFile).use { writer ->
            writer.write("RefRange\tSampleName\tHapID\tHapCount\tTotalCount\n")
            referenceRangeList.forEach { refRange ->
                val hapids = rangeToHapIds[refRange]
                if (hapids != null) {
                    val totalCount = referenceRangeCounts[refRange]
                    val hapidSamplePairs = hapids.flatMap {hapid -> hapIdToSampleGamete[hapid]!!.map { sample -> Pair(hapid, sample)}}
                    hapidSamplePairs.forEach { (hapid, sample) -> writer.write("$refRange\t$sample\t$hapid\t${hapIdCounts[hapid]?:0}\t$totalCount\n")}
                }
            }
        }
    }

    /**
    * Function that returns a map of reference range to read mapping counts for that reference range.
    * */
    fun getReadMappingCountsByReferenceRange(graph: HaplotypeGraph, readMappings: Map<List<String>, Int>): Map<ReferenceRange, Int> {
        val readsByReferenceRange = readMappingByRange(readMappings, graph)
        return readsByReferenceRange.entries.associate {(refrange, hapidcounts) -> refrange to hapidcounts.values.sum()}
    }

    private fun readMappingByRange(readCounts: Map<List<String>, Int>, graph: HaplotypeGraph): Map<ReferenceRange, Map<List<String>, Int>> {
        val hapid2Refrange = graph.hapIdToRefRangeMap()
        //groups readCounts for the entire genome into separate maps for each reference range
        //since some hapids map to more than one reference range, assign each hapid set to the reference range with the most hapids in the set

        return readCounts.entries.groupBy { referenceRangeForHapidList(it.key, hapid2Refrange) }
            .mapValues { (_,mapEntries) -> mapEntries.associateBy({it.key}, {it.value}) }
    }

    private fun referenceRangeForHapidList(hapidList: List<String>, hapid2Refrange: Map<String, List<ReferenceRange>>): ReferenceRange {
        //this handles the case where a hapid is in more than one reference range
        val referenceRangeList = hapidList.mapNotNull { hapid2Refrange[it] }.flatten()
        val referenceRangeCount = referenceRangeList.groupingBy { it }.eachCount()
        return referenceRangeCount.maxBy {it.value}.key
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
                       hapIdToSampleGamete: Map<String, List<SampleGamete>>,
                       referenceRangeCounts: Map<ReferenceRange, Int>) {
        myLogger.info("Writing out counts")
        bufferedWriter(outputFile).use { writer ->
            writer.write("refRange\t${targetSampleName}_HapID\t${targetSampleName}_HapCount\tHighestAltCount\tDifference\tTotalCount\tOtherHapCounts\n")
            rangeToHapId.keys.sorted().forEach { range ->
                val outputString = buildOutputStringForHapIdsInRefRange(rangeToHapId, range, hapIdToSampleGamete, hapIdCounts, targetSampleName, referenceRangeCounts)
                writer.write(outputString)
            }
        }
    }

    fun buildOutputStringForHapIdsInRefRange(
        rangeToHapId: Map<ReferenceRange, List<String>>,
        range: ReferenceRange,
        hapIdToSampleGamete: Map<String, List<SampleGamete>>,
        hapIdCounts: Map<String, Int>,
        targetSampleName: String,
        referenceRangeCounts: Map<ReferenceRange, Int>
    ) : String {
        val hapIds = rangeToHapId[range]!!
        val targetHapId = findTargetHapIdForSample(hapIds, hapIdToSampleGamete, targetSampleName)

        val nonTargetCounts = computeCountsForNonTargetHapids(hapIds, targetHapId, hapIdCounts)

        val highestAlt = if (nonTargetCounts.isEmpty()) 0 else nonTargetCounts.first().first

        val targetCount = hapIdCounts.getOrDefault(targetHapId, 0)
        return "$range\t$targetHapId\t${targetCount}\t${highestAlt}\t${targetCount - highestAlt}\t" +
                "${referenceRangeCounts[range]?:0}\t" +
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