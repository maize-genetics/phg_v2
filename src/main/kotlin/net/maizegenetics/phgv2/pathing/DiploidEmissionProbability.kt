package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import org.apache.commons.math3.distribution.BinomialDistribution

/**
 * Calculates the natural log (ln) of the emission probability for a pair of haplotypes.
 */
class DiploidEmissionProbability(val readMap: Map<ReferenceRange, Map<List<String>, Int>>, val graph: HaplotypeGraph, val probabilityCorrect: Double) {
    private var myCurrentRange = ReferenceRange("NA", 0, 0)
    //a map of a pair of haplotypes -> ln probability
    private var rangeLnProbabilities = mutableMapOf<Pair<String,String>, Double>()
    //a map of SampleGamete -> haplotype for myCurrentRange
    private var sampleToHaplotype = mapOf<SampleGamete, String>()
    private var defaultProbability = 1e-10
    private var missingHaplotypeProbability = 1e-10

    /**
     * Returns the natural log of the emission probability for a given state and [ReferenceRange].
     * The emission probability is the probability of observing the read counts given that the generating sample carries
     * the haplotypes of the SampleGametes. The state is a pair of SampleGametes.
     */
    fun lnProbObsGivenState(state: Pair<SampleGamete, SampleGamete>, refrange: ReferenceRange): Double {
        if (myCurrentRange != refrange) {
            myCurrentRange = refrange
            assignRangeProbabilities()
        }
        return rangeLnProbabilities[state] ?: defaultProbability
    }

    private fun assignRangeProbabilities() {
        //Todo evaluate whether haplotype pair with count = 0 is the same as no haplotype

        val sampleToHapidMap = createSampleToHaplotypeMap()

        //This next set is needed because a sample in a reference which has no haplotype is expected to not have any reads
        //whereas a sample with a hapid in a reference range is expected to be hit by many of the reads
        val samplesInRefrange = sampleToHapidMap.keys
        val haplotypesInRefrange = sampleToHapidMap.values.toSet().toList()

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
        //iterate over all possible pairs of haplotypes in myCurrentRange
        //The list of haplotypes has to come from the graph, since not all haplotypes will have reads mapped to them.

        val haplotypePairProbability = mutableMapOf<UnorderedHaplotypePair, Double>()
        for (ndx1 in haplotypesInRefrange.indices) {
            for (ndx2 in ndx1 until haplotypesInRefrange.size) {

            }
        }

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

    private fun createSampleToHaplotypeMap(): Map<SampleGamete, String> {
        //converts a map of hapid -> list of SampleGametes to a map of
        //SampleGamete -> hapid. This works because each SampleGamete has only one hapid.
        return graph.hapIdToSampleGametes(myCurrentRange).entries
            .flatMap{ (hapid, sampleList) -> sampleList.map { Pair(it, hapid) } }.toMap()
    }

    fun getLnMissingHaplotypeProbability(): Double = missingHaplotypeProbability

    private fun haplotypePairProbability(haplotypes: UnorderedHaplotypePair, readCounts: Map<Set<String>, Int>): Double {
        //readCounts keys are a set rather than a list for faster contains method
        val hapPair = haplotypes.haplotypePair
        if (haplotypes.haplotypePair.first == haplotypes.haplotypePair.second) {
            val totalCount = readCounts.values.sum()
            val firstCount = readCounts.filter { (hapset, _) -> hapset.contains(hapPair.first) }.values.sum()
            BinomialDistribution(totalCount, probabilityCorrect).probability(firstCount)
        } else {
            val firstNotSecondCount = readCounts.filter { (hapset, _) -> hapset.contains(hapPair.first) && !hapset.contains(hapPair.second)}
                .values.sum()
            val secondNotFirstCount = readCounts.filter { (hapset, _) -> !hapset.contains(hapPair.first) && hapset.contains(hapPair.second)}
                .values.sum()
            val firstAndSecondCount = readCounts.filter { (hapset, _) -> hapset.contains(hapPair.first) && hapset.contains(hapPair.second)}
                .values.sum()
            val neitherFirstNorSecondCount = readCounts.filter { (hapset, _) -> !hapset.contains(hapPair.first) && !hapset.contains(hapPair.second)}
                .values.sum()

            (0..firstAndSecondCount).map { multinomialProbability(intArrayOf(firstNotSecondCount + it, secondNotFirstCount + firstAndSecondCount - it, neitherFirstNorSecondCount),
                doubleArrayOf(halfProb, halfProb, pErr))}.sum()
        }

    }

    }

    data class UnorderedHaplotypePair(val haplotypePair: Pair<String,String>) {
        override fun equals(other: Any?): Boolean {
            return if (other is UnorderedHaplotypePair) {
                if (other.haplotypePair == haplotypePair) true
                else if (other.haplotypePair.first == haplotypePair.second && other.haplotypePair.second == haplotypePair.first) true
                else false
            }
            else false
        }
    }
}

