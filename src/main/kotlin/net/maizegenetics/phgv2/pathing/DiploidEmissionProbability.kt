package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.api.ReferenceRange

class DiploidEmissionProbability(val readMap: Map<ReferenceRange, Map<List<String>, Int>>) {
    private var myCurrentRange = ReferenceRange("NA", 0, 0)
    private var rangeLnProbabilities = mutableMapOf<Pair<String,String>, Double>()
    private var defaultProbability = 1e-10
    private var missingHaplotypeProbability = 1e-10

    /**
     * Returns the emission probability for a given state index and range index.
     * The range index is the zero-based index into the ordered list of ranges for a HaplotypeGraph.
     * The state index is the zero-based index into an ordered list of states, where each state is an ordered pair
     * of HaplotypeNodes returned by the method [nodePairs]. The emission probability for a given node pair is calculated
     * by the method [nodePairProbabilityList].
     */
    fun lnProbObsGivenState(state: Pair<String, String>, refrange: ReferenceRange): Double {
        if (myCurrentRange != refrange) {
            myCurrentRange = refrange
            assignRangeProbabilities()
        }
        return rangeLnProbabilities[state] ?: defaultProbability
    }

    private fun assignRangeProbabilities() {
        //Todo evaluate whether haplotype pair with count = 0 is the same as no haplotype
        val readCounts = readMap[myCurrentRange] ?: mapOf()
        require(readCounts.isNullOrEmpty()) {"No haplotypes in $myCurrentRange"}

//        val currentNodes: List<HaplotypeNode> = rangeToNodes.get(myCurrentRange)
//            ?: throw IllegalArgumentException("Range at chr ${myCurrentRange.chromosome().name}, pos ${myCurrentRange.start()} has no nodes.")

        //get the list of read node matches for this reference range
//        val currentReadMatches = readHapids.get(myCurrentRange)

        //convert that list to counts as Long -> Int, where long is two hapids encoded
        //create a map of hapid pair to count of reads for that pair

//        val nodePairCounts: Map<Long, Int> = countsFromHapidList(currentReadMatches)
        //get the list of node pairs
//        val myNodePairs = nodePairs(currentNodes)
        val totalCount = readCounts.values.sum()

        //if total count = 0, assign probability of 1/n to all nodes
        // otherwise, use the method nodePairProbabilityList(...)
        if (totalCount == 0) {
            val prob = 1.0 / myNodePairs.size.toDouble()
            rangeLnProbabilities = DoubleArray(myNodePairs.size, { i -> prob })
        } else {
            val nodeProbList = nodePairProbabilityList(currentNodes, currentReadMatches.toList(), probabilityCorrect)
            rangeLnProbabilities = nodeProbList.toDoubleArray()
        }
    }

    private fun countOfHapids(hapidPair: Pair<String,String>, readCounts: Map<List<String>, Int>): Int {
        val hapids = listOf(hapidPair.first, hapidPair.second)
        return readCounts.entries.filter { (hapidList,_) -> hapidList.containsAll(hapids) }.sumOf { it.value }
    }

    fun getLnMissingHaplotypeProbability(): Double = missingHaplotypeProbability

}