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
    //minimum probability must be greater then 0 since we need to take the ln of it
    private val minProbability = Double.MIN_VALUE
    private val lnMinProbability = ln(minProbability)

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

    /**
     * Generates a map of all possible SampleGamete order pairs to their emission probabilities. This is done
     * once per reference range to avoid repeating the calculations.
     */
    private fun assignRangeProbabilities() {
        sampleToHaplotype = graph.sampleGameteToHaplotypeId(myCurrentRange)

        val haplotypesInRefrange = sampleToHaplotype.values.toSet().toList()

        val readCounts = readMap[myCurrentRange] ?: mapOf()

        require(readCounts.isNotEmpty()) {"No haplotypes in $myCurrentRange"}
        val readSetCounts = readCounts.mapKeys { (haplist, _) -> haplist.toSet() }

        val probabilityMap = mutableMapOf<UnorderedHaplotypePair, Double>()
        for (ndx1 in haplotypesInRefrange.indices) {
            for (ndx2 in ndx1 until haplotypesInRefrange.size) {
                val haplotypePair = UnorderedHaplotypePair(Pair(haplotypesInRefrange[ndx1], haplotypesInRefrange[ndx2]))
                probabilityMap[haplotypePair] = lnHaplotypePairProbability(haplotypePair, readSetCounts)
            }
        }

        //if there are any null haplotypes in this reference range add pairs for each (haplotype,null) and (null,null)
        val anyNullHaplotypes = sampleGametesInGraph.any { sampleToHaplotype[it] == null }
        if (anyNullHaplotypes) {
            for (haplotype in haplotypesInRefrange) {
                val haplotypePair = UnorderedHaplotypePair(Pair(haplotype, null))
                probabilityMap[haplotypePair] = lnHaplotypePairProbability(haplotypePair, readSetCounts)
            }
            val haplotypePair = UnorderedHaplotypePair(Pair(null, null))
            probabilityMap[haplotypePair] = lnHaplotypePairProbability(haplotypePair, readSetCounts)
        }

        //assign natural log of the probabilities to rangeLnProbabilities
        rangeLnProbabilities = probabilityMap
    }

    private fun lnHaplotypePairProbability(haplotypes: UnorderedHaplotypePair, readCounts: Map<Set<String>, Int>): Double {
        val halfProb = probabilityCorrect / 2.0
        val pErr = 1.0 - probabilityCorrect
        //probability the a read mapping to A maps also maps to B, should be user settable or potentially has a different value by reference range
        val pBoth = 0.3
        //assign a low probability for now, maybe should be user settable?
        val diploidProb = halfProb - pBoth - pErr

        //readCounts keys are a Set rather than a list for faster contains method
        val hapPair = haplotypes.haplotypePair
        return if ( hapPair.first != null && (hapPair.second == hapPair.first || hapPair.second == null)) {
            //If the genotype is homozygous or one of the haplotypes is a null haplotype, then the probability of
            // observing these read counts is the same as the haploid genotype case.
            // That is P(obs|homozygous genotype) = Binom(totalCount, # of reads with this haplotype, probabilityCorrect).
            // In words, it is the probability of observing m reads that map to this haplotype out of total of n reads
            // when the probability of a success equals probabilityCorrect
            val totalCount = readCounts.values.sum()
            val firstCount = readCounts.filter { (hapset, _) -> hapset.contains(hapPair.first) }.values.sum()
            BinomialDistribution(totalCount, probabilityCorrect).probability(firstCount)
        } else if (hapPair.first == null && hapPair.second != null) {
            //the same case as the first if. This could have been included in that condition but this is
            // easier to understand. (and code correctly)
            val totalCount = readCounts.values.sum()
            val firstCount = readCounts.filter { (hapset, _) -> hapset.contains(hapPair.second) }.values.sum()
            BinomialDistribution(totalCount, probabilityCorrect).probability(firstCount)
        } else {
            //firstAndSecondCount is the number of reads mapping to both haplotype. firstNotSecondCount is
            // the number of reads mapping to the first but not the second, etc. Note that these classes are mutually
            // exclusive and that they sum to the total counts.
            val firstNotSecondCount = readCounts.filter { (hapset, _) -> hapset.contains(hapPair.first) && !hapset.contains(hapPair.second)}
                .values.sum()
            val secondNotFirstCount = readCounts.filter { (hapset, _) -> !hapset.contains(hapPair.first) && hapset.contains(hapPair.second)}
                .values.sum()
            val firstAndSecondCount = readCounts.filter { (hapset, _) -> hapset.contains(hapPair.first) && hapset.contains(hapPair.second)}
                .values.sum()
            val neitherFirstNorSecondCount = readCounts.filter { (hapset, _) -> !hapset.contains(hapPair.first) && !hapset.contains(hapPair.second)}
                .values.sum()

            /*
            *Consider two haplotypes A and B. Given that the reads come from a heterozygous AB individual,
            * let a = the number of reads mapping to A but not B
            * b = the number of reads mapping to B but not A
            * c = the number of reads mapping to both A and B
            * n = the number of reads mapping to neither A nor B
            * then the emission probability P(data | AB) = multinomialProbability((a,b,c,n), (PA, PB, PC, PD))
            * where PC = pBoth
            * PN = pErr
            * PA = PB = diploidProb
            * This is the simplest way to model PA and PB. Methods based on data could be used.
            * */
            lnMultinomialProbability(intArrayOf(firstNotSecondCount, secondNotFirstCount, firstAndSecondCount, neitherFirstNorSecondCount),
                doubleArrayOf(diploidProb, diploidProb, pBoth, pErr))

        }

    }

    /**
     * @param [counts] The counts of each of set of classes
     * @param [probabilities] The probability of each class
     * @return The probability of an array of counts of some classes given the probability of each class,
     * The size of the counts array and the probabilities array are expected to be equal.
     */
    fun lnMultinomialProbability(counts: IntArray, probabilities: DoubleArray): Double {
        if (counts.size != probabilities.size) throw java.lang.IllegalArgumentException("multinomialProbability error: counts and probabilities arrays do not have the same size.")
        val totalCount = counts.sum()

        val logprod = counts.indices.sumOf { counts[it] * ln(probabilities[it]) }
        val numerator = logFactorial(totalCount)
        val denominator = counts.sumOf { logFactorial(it) }
        return numerator - denominator + logprod
    }

    /**
     *  Calculates the (natural) log factorial of any positive integer using the exact value for 0 to 10 and
     *  Stirlings approximation for integers greater than 10. The formula is taken from the
     *  Wikipedia article for Stirlings approximation
     *  @param [intval] an integer
     *  @return   the natural log of the factorial of intval
     */
    private fun logFactorial(intval: Int): Double {
        return if (intval <= 10) smallFactorials[intval] else {
            val n = intval.toDouble()
            n * ln(n) + 0.5 * ln(2.0 * Math.PI * n) - n
        }
    }

    //factorials of 0 to 10
    private val smallFactorials: DoubleArray = doubleArrayOf(1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0)
        .map { ln(it) }.toDoubleArray()


    data class UnorderedHaplotypePair(val haplotypePair: Pair<String?, String?>) {
        /**
         * Kotlin's equality: null == null is true. So, equal works with null haplotypes as intended here.
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

