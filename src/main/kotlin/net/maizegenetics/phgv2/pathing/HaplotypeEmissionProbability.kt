package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.api.ReferenceRange
import java.lang.IllegalArgumentException
import org.apache.commons.math3.distribution.BinomialDistribution
import kotlin.math.ln

/**
 * @author Peter Bradbury
 *
 * The HaplotypeEmissionProbability class extends EmissionProbability. It provides methods to calculate an
 * emission probability that is used in a hidden Markov model (HMM) for imputing haplotypes. The constructor takes the
 * following arguments:
 * @param refRangeToHapListMap  A map of ReferenceRange -> list of haplotypes in that range
 * @param readMap   A map of ReferenceRange -> list of HaplotypeListCounts (a count of the number of reads mapping to all
 * of the haplotypes in the list)
 * @param pCorrect  The probability that a read has been aligned correctly
 *
 */

class HaplotypeEmissionProbability(val refRangeToHapListMap : Map<ReferenceRange, List<String>>,
                                   val readMap : Map<ReferenceRange, Map<List<String>, Int>>,
                                   val pCorrect : Double) {
    private var currentRefRange = ReferenceRange("none",0,0)
    private var currentEmissionProbabilities = mapOf<String, Double>()
    private var nullProbability : Double = -10.0

    /**
     * Returns the natural log of the probability of observing the read mapping counts for this range
     * given the [haplotype]
     */
    fun getLnProbObsGivenState(haplotype: String?, refRange: ReferenceRange): Double {
        //refrange is a 0..n index into nodeTree.keys
        //haplotype is an index into the List<HaplotypeNode>> for this ReferenceRange
        if (currentRefRange != refRange) {
            currentRefRange = refRange
            currentEmissionProbabilities = calculateLnHaplotypeProbabilities()

        }

        return currentEmissionProbabilities[haplotype] ?: nullProbability
    }


    private fun calculateLnHaplotypeProbabilities() : Map<String, Double> {
        /*
        * P_emission = P(obs|haplotype)
        * = binom(#trials, #successes, P_success)
        * where #trials = number of reads
        *   #successes = number of reads hitting the haplotype
        *   P_success = pcorrect
        *
        * assign minProbability for haplotype index = -1
        * */

        val myRefRange = currentRefRange
        val haplotypes = refRangeToHapListMap[myRefRange]
        check(!haplotypes.isNullOrEmpty()) {IllegalArgumentException("$myRefRange contains no haplotypes.")}

        val myHapMappings = readMap[myRefRange]

        if (myHapMappings.isNullOrEmpty()) {
            //if there no reads then all haplotypes get assigned the same probability, it does not matter what it is.
            //note that 0.0 is ln(1), so this actually assigns a probability of 1 not zero. One interpretation is that
            //if there are no reads then every haplotype (including null haplotypes) will have 0 reads with a probability of 1.
            nullProbability = 0.0
            return haplotypes.associateWith { 0.0 }

        } else {
            val numberOfReads = myHapMappings.map { it.value }.sum()
            val binom = BinomialDistribution(numberOfReads, pCorrect)

            //this next line maps hapid -> sum of counts for sets containing that hapid
            val hapidCountMap = myHapMappings.flatMap { hapCounts -> hapCounts.key.map { Pair(it, hapCounts.value)} }
                .groupingBy { it.first }
                .fold(0) { sum, pr -> sum + pr.second }

            //since a null haplotype should have 0 counts, assign it binom.probability(0)
            nullProbability = ln( binom.probability( 0) )
            return haplotypes.associateWith { hapid ->
                val prob = binom.probability(hapidCountMap[hapid] ?: 0).coerceAtLeast(Double.MIN_VALUE)
                ln( prob )
            }

        }

    }

}

