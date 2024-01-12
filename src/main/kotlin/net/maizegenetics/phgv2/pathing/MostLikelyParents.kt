package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete

/**
 * Given a set of read mappings, finds the most likely parents in a stepwise fashion. The parent with the most reads
 * is selected as a likely parent. Then, the subset of reads not mapping to that parent is created.
 * The parent with the most reads from that subset is then added to the list of likely parents.
 * The process is repeated until the number of likely parents = maxParents, coverage is equal to or greater than minCoverage,
 * or until all parents have been added to the list, whichever occurs first. Coverage is calculated as the number
 * of reads mapping to any of the likely parents divided by the total number of reads.
 */
class MostLikelyParents(val hapGraph: HaplotypeGraph) {
    //list of parents, which are the sampleGametes present in the HaplotypeGraph
    private val myParentList: List<SampleGamete> = hapGraph.sampleGametesInGraph().toList()

    data class ParentStats(val parent: SampleGamete, val readCount: Int, val coverage: Double)

    /**
     * Finds the most likely parents for a set of hapid counts from read mappings.
     */
    fun findMostLikelyParents(refRangeToHapIdListCounts: Map<ReferenceRange, Map<List<String>, Int>>, maxParents: Int, minCoverage: Double) : List<ParentStats> {
        //convert refRangeToHapIdSetCounts to a Map<ReferenceRange,List<HapIdSetCount>>
        //use only the counts for ranges that are present in myParentToHapidMap
        val rangesInGraph = hapGraph.ranges()
        var filteredCounts = refRangeToHapIdListCounts.filter { (refrange, _) ->
            rangesInGraph.contains(refrange) }
        val bestParentList = mutableListOf<ParentStats>()
        var iteration = 0
        var coverage = 0.0

        //total Count of reads, which is the sum of the read counts for each reference range (outer sum)
        val totalCount = filteredCounts.map {(_,setCounts) -> setCounts.values }.flatten().sum()
        var cumulativeCount = 0

        while (iteration < maxParents && coverage < minCoverage) {
            iteration++
            var bestParent = SampleGamete("none", 0)
            var highestCount = 0
            for (parent in myParentList) {
                //get count for this parent
                val parentCount =
                    filteredCounts.keys.sumOf { refrange ->
                        val hapid = hapGraph.sampleToHapId(refrange, parent)
                        if (hapid == null) 0 else
                            filteredCounts[refrange]!!.filter { it.key.contains(hapid) }.entries.sumOf { it.value }
                    }
                if (parentCount > highestCount) {
                    highestCount = parentCount
                    bestParent = parent
                }
            }

            cumulativeCount += highestCount
            coverage = cumulativeCount / totalCount.toDouble()
            bestParentList.add(ParentStats(bestParent, highestCount, coverage))

            //rebuild the filteredList using only hapIdSetCounts that do not contain the best parent
            //this will be used for the next round
            filteredCounts = filteredCounts.keys.map { refrange ->
                //for each refrange count parent use
                val hapid = hapGraph.sampleToHapId(refrange, bestParent)
                val filteredList = filteredCounts[refrange]!!.filter { !it.key.contains(hapid) }
                Pair(refrange, filteredList)
            }.filter { it.second.isNotEmpty() }.toMap()
        }

        return bestParentList
    }

}

