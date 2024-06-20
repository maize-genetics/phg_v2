package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.getBufferedWriter
import java.io.File

class MlReadInputArray : CliktCommand(help = "write an ML input array from a read mapping file") {

    val hvcfDir by option("--hvcf-dir", help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        .required()

    val readFile by option(help = "Name of the read mapping file to convert into an ML input array")
        .required()

    val outputFile by option(help = "The to which the ML will be written")
        .required()

    val shuffle by option(help = "if the shuffle flag is set, the order of reads with a reference range are shuffled")
        .flag()

    override fun run() {

        //build graph
        val hvcfFilenameList = File(hvcfDir).listFiles()
            .filter { it.name.endsWith(".h.vcf") ||it.name.endsWith(".h.vcf.gz")}
            .map{it.path}

        val graph = HaplotypeGraph(hvcfFilenameList)

        //read in readMappings
        val readMappings = AlignmentUtils.importReadMapping(readFile)

        //create a map of ReferenceRange -> Map<List<String>, Int> from the readMappings
        val readMappingByRange = readMappingByRange(readMappings, graph)

        //get gamete sample names sorted by most likely parent
        val parentStats = MostLikelyParents(graph).findMostLikelyParents(readMappingByRange, graph.numberOfSamples(), 1.0)

        //write the ML input array
        writeMLArray(readMappingByRange, parentStats, graph)

    }

    private fun readMappingByRange(readCounts: Map<List<String>, Int>, graph: HaplotypeGraph): Map<ReferenceRange, Map<List<String>, Int>> {
        val hapid2Refrange = graph.hapIdToRefRangeMap()
        //groups readCounts for the entire genome into separate maps for each reference range
        //since some hapids map to more than one reference range, assign each hapid set to the reference range with the most hapids in the set

        return readCounts.entries.groupBy { referenceRangeForHapidList(it.key, hapid2Refrange) }
            .mapValues { (_,mapEntries) -> mapEntries.associateBy({it.key}, {it.value}) }
    }

    private fun referenceRangeForHapidList(hapidList: List<String>, hapid2Refrange: Map<String, List<ReferenceRange>>): ReferenceRange {
        val referenceRangeList = hapidList.mapNotNull { hapid2Refrange[it] }.flatten()
        val referenceRangeCount = referenceRangeList.groupingBy { it }.eachCount()
        return referenceRangeCount.maxBy {it.value}.key
    }

    fun writeMLArrayxxx(countMap: Map<ReferenceRange, Map<List<String>, Int>>, parentStats: List<MostLikelyParents.ParentStats>, graph: HaplotypeGraph) {
        //sampleName to column map
        val sep = ","
        val sampleGameteIndex = parentStats.mapIndexed { index, parentStats ->  Pair(parentStats.parent, index)}.toMap().toMutableMap()

        //add any sample gametes from graph not in parentStats
        for (gamete in graph.sampleGametesInGraph()) {
            if (!sampleGameteIndex.keys.contains(gamete)) sampleGameteIndex[gamete] = sampleGameteIndex.size
        }
        val numberOfGametes = sampleGameteIndex.size

        val IndexToSampleGamete = sampleGameteIndex.entries.associateBy({it.value},{it.key})

        getBufferedWriter(outputFile).use { mlWriter ->
            //each line of the output will be a single read
            //the content will be 0's and 1's indicating which gametes are in the hapid set
            //the column labels will be the sample names
            //the row labels will be the reference range followed by a sequential number

            //header
//            mlWriter.write(parentStats.map { it.parent.name }.joinToString(sep, prefix = "name$sep", postfix = "\n"))
            val headerNames = Array(numberOfGametes) {IndexToSampleGamete[it]}
            mlWriter.write("name$sep${headerNames.joinToString(sep)}\n")

            var rowIndex = 0
            for (refrange in graph.ranges()) {
                val refrangeCountMap = countMap[refrange]
                if (refrangeCountMap != null) {
                    val hapidToGamete = graph.hapIdToSampleGametes(refrange)
                    for ((hapidList, count) in refrangeCountMap) {
                        val gameteArray = IntArray(numberOfGametes) {0}
                        //change gameteArray value to 1 for all gametes that have hapids in the list
                        for (hapid in hapidList) {
                            val gameteList = hapidToGamete[hapid]
                            check(gameteList != null) {"$hapid was in count map but not in hapid to gamete map"}
                            for (gamete in gameteList) {
                                val ndx = sampleGameteIndex[gamete]!!
                                gameteArray[ndx] = 1
                            }
                        }
                        val outputString = gameteArray.joinToString(sep)
                        repeat(count) {
                            mlWriter.write("$refrange:$rowIndex$sep$outputString\n")
                            rowIndex++
                        }
                    }
                }
            }

        }
    }

    fun writeMLArray(countMap: Map<ReferenceRange, Map<List<String>, Int>>, parentStats: List<MostLikelyParents.ParentStats>, graph: HaplotypeGraph) {
        //sampleName to column map
        val sep = ","
        val sampleGameteIndex = parentStats.mapIndexed { index, parentStats ->  Pair(parentStats.parent, index)}.toMap().toMutableMap()

        //add any sample gametes from graph not in parentStats
        for (gamete in graph.sampleGametesInGraph()) {
            if (!sampleGameteIndex.keys.contains(gamete)) sampleGameteIndex[gamete] = sampleGameteIndex.size
        }
        val numberOfGametes = sampleGameteIndex.size

        val IndexToSampleGamete = sampleGameteIndex.entries.associateBy({it.value},{it.key})

        getBufferedWriter(outputFile).use { mlWriter ->
            //each line of the output will be a single read
            //the content will be 0's and 1's indicating which gametes are in the hapid set
            //the column labels will be the sample names
            //the row labels will be the reference range followed by a sequential number

            //header
//            mlWriter.write(parentStats.map { it.parent.name }.joinToString(sep, prefix = "name$sep", postfix = "\n"))
            val headerNames = Array(numberOfGametes) {IndexToSampleGamete[it]}
//            mlWriter.write("range$sep${headerNames.joinToString(sep)}\n")
            mlWriter.write("Name\tTokenized\n")

            var rowIndex = 0
            for (refrange in graph.ranges()) {
                val outputRows = mutableListOf<String>()
                val refrangeCountMap = countMap[refrange]
                if (refrangeCountMap != null) {
                    val hapidToGamete = graph.hapIdToSampleGametes(refrange)
                    for ((hapidList, count) in refrangeCountMap) {
                        val gameteArray = IntArray(numberOfGametes) {0}
                        //change gameteArray value to 1 for all gametes that have hapids in the list
                        //or list the gametes
                        for (hapid in hapidList) {
                            val gameteList = hapidToGamete[hapid]
                            check(gameteList != null) {"$hapid was in count map but not in hapid to gamete map"}
                            for (gamete in gameteList) {
                                val ndx = sampleGameteIndex[gamete]!!
                                gameteArray[ndx] = 1
                            }
                        }
                        val dummyNames = hapidList.mapNotNull { hapid -> hapidToGamete[hapid] }.flatten()
                            .map { "gamete_${sampleGameteIndex[it]}" }.joinToString(sep, prefix = "\"", postfix = "\"")
                        repeat(count) {
                            outputRows.add(dummyNames)
                        }
                    }
                    if (shuffle) outputRows.shuffle()
                    outputRows.forEach {
                        mlWriter.write("Read$rowIndex\t$it\n")
                        rowIndex++
                    }
                }
            }

        }
    }
}