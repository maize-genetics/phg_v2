package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import org.apache.commons.math3.distribution.BinomialDistribution
import kotlin.math.ln

/**
 * Calculates the natural log (ln) of the emission probability for a pair of haplotypes.
 */
class DiploidEmissionProbability(val readMap: Map<ReferenceRange, Map<List<String>, Int>>, val graph: HaplotypeGraph, val probabilityCorrect: Double) {
    private var myCurrentRange = ReferenceRange("NA", 0, 0)
    //a map of a pair of haplotypes -> ln probability
    private var rangeLnProbabilities = mapOf<UnorderedHaplotypePair, Double>()
    //a map of SampleGamete -> haplotype for myCurrentRange
    private var sampleToHaplotype = mapOf<SampleGamete, String>()
    private var defaultProbability = ln(1e-10)
    private val sampleGametesInGraph = graph.sampleGametesInGraph()
    private val noReadCountProbability = 0.0
    private var rangeHasReads = true

    /**
     * Returns the natural log of the emission probability for a given state and [ReferenceRange].
     * The emission probability is the probability of observing the read counts given that the generating sample carries
     * the haplotypes of the SampleGametes. The state is a pair of SampleGametes.
     *
     * If the range has no reads that mapped to it, then all SampleGamete pairs have the same probability. So,
     * return noReadCountProbability when there are no reads. Since the emission probability only needs to be
     * proportional to the probability of the observations (read counts) given the state, returning 0.0 for all pairs works fine.
     */
    fun lnProbObsGivenState(state: Pair<SampleGamete, SampleGamete>, refrange: ReferenceRange): Double {

        if (myCurrentRange != refrange) {
            myCurrentRange = refrange
            val readsForRange = readMap[refrange]
            rangeHasReads = if (readsForRange == null) false else readsForRange.values.sum() > 0
            if (rangeHasReads) assignRangeProbabilities()
        }

        if (!rangeHasReads) return noReadCountProbability

        val haplotypePair = Pair(sampleToHaplotype[state.first], sampleToHaplotype[state.second])

        //convert the pair of sample gametes to an UnorderedHaplotypePair
        val unorderedHaplotypes = UnorderedHaplotypePair(haplotypePair)
        return rangeLnProbabilities[unorderedHaplotypes] ?: defaultProbability
    }

    private fun assignRangeProbabilities() {
        sampleToHaplotype = createSampleToHaplotypeMap()

        val haplotypesInRefrange = sampleToHaplotype.values.toSet().toList()

        val readCounts = readMap[myCurrentRange] ?: mapOf()

        require(!readCounts.isNullOrEmpty()) {"No haplotypes in $myCurrentRange"}
        val readSetCounts = readCounts.mapKeys { (haplist, _) -> haplist.toSet() }

        val probabilityMap = mutableMapOf<UnorderedHaplotypePair, Double>()
        for (ndx1 in haplotypesInRefrange.indices) {
            for (ndx2 in ndx1 until haplotypesInRefrange.size) {
                val haplotypePair = UnorderedHaplotypePair(Pair(haplotypesInRefrange[ndx1], haplotypesInRefrange[ndx2]))
                probabilityMap.put(haplotypePair, haplotypePairProbability(haplotypePair, readSetCounts))
            }
        }

        //if there are any null haplotypes in this reference range add pairs for each (haplotype,null) and (null,null)
        val anyNullHaplotypes = sampleGametesInGraph.any { sampleToHaplotype[it] == null }
        if (anyNullHaplotypes) {
            for (haplotype in haplotypesInRefrange) {
                val haplotypePair = UnorderedHaplotypePair(Pair(haplotype, null))
                probabilityMap.put(haplotypePair, haplotypePairProbability(haplotypePair, readSetCounts))
            }
            val haplotypePair = UnorderedHaplotypePair(Pair(null, null))
            probabilityMap.put(haplotypePair, haplotypePairProbability(haplotypePair, readSetCounts))
        }

        //take the natural log of the probabilities
        rangeLnProbabilities = probabilityMap.mapValues { (_, pr) -> ln(pr) }
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

    private fun haplotypePairProbability(haplotypes: UnorderedHaplotypePair, readCounts: Map<Set<String>, Int>): Double {
        //todo deal with null haplotypes
        val halfProb = probabilityCorrect / 2
        val pErr = 1 - probabilityCorrect

        //readCounts keys are a set rather than a list for faster contains method
        val hapPair = haplotypes.haplotypePair
        return if (haplotypes.haplotypePair.first == haplotypes.haplotypePair.second) {
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

    /**
     * @param [counts] The counts of each of set of classes
     * @param [probabilities] The probability of each class
     * @return The probability of an array of counts of some classes given the probability of each class,
     * The size of the counts array and the probabilities array are expected to be equal.
     */
    private fun multinomialProbability(counts: IntArray, probabilities: DoubleArray): Double {
        if (counts.size != probabilities.size) throw java.lang.IllegalArgumentException("multinomialProbability error: counts and probabilities arrays do not have the same size.")
        val N = counts.sum()

        val logprod = counts.indices.map { counts[it] * Math.log(probabilities[it]) }.sum()
        val numerator = logFactorial(N)
        val denom = counts.map { logFactorial(it) }.sum()
        val logprob = numerator - denom + logprod
        return Math.pow(Math.E, logprob)
    }

    /**
     *  Calculates the log factorial of any positive integer using the exact value for 0 to 10 and
     *  Stirlings approximation for integers greater than 10. The formula is taken from the
     *  Wikipedia article for Stirlings approximation
     *  @param [intval] an integer
     *  @return   the natural log of the factorial of intval
     */
    private fun logFactorial(intval: Int): Double {
        if (intval <= 10) return smallFactorials[intval] else {
            val n = intval.toDouble()
            return n * ln(n) + 0.5 * ln(2.0 * Math.PI * n) - n
        }
    }

    //factorials of 0 to 10
    private val smallFactorials: DoubleArray = doubleArrayOf(1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0)
        .map { ln(it) }.toDoubleArray()


    data class UnorderedHaplotypePair(val haplotypePair: Pair<String?, String?>) {
        /**
         * Kotlin equality: null == null is true. So, equal works with null haplotypes as intended here.
         */
        override fun equals(other: Any?): Boolean {
            return if (other is UnorderedHaplotypePair) {
                if (other.haplotypePair == haplotypePair) true
                else if (other.haplotypePair.first == haplotypePair.second && other.haplotypePair.second == haplotypePair.first) true
                else false
            }
            else false
        }

        override fun hashCode(): Int {
            return haplotypePair.first.hashCode() + haplotypePair.second.hashCode()
        }
    }
}

