package net.maizegenetics.phgv2.pathing

import com.google.common.collect.*
import net.maizegenetics.phgv2.api.ReferenceRange
import java.lang.IllegalArgumentException
import org.apache.commons.math3.distribution.BinomialDistribution

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
                                   val readMap : Map<ReferenceRange, List<HaploidPathFinding.HaplotypeListCount>>,
                                   val pCorrect : Double) {
    var currentRefRange = ReferenceRange("none",0,0)
    var currentEmissionProbabilities = mapOf<String, Double>()
    var minProbability : Double = 0.0

    fun getLnProbObsGivenState(haplotype: String, refRange: ReferenceRange): Double {
        //refrange is a 0..n index into nodeTree.keys
        //haplotype is an index into the List<HaplotypeNode>> for this ReferenceRange

        if (currentRefRange != refRange) {
            currentRefRange = refRange
            currentEmissionProbabilities = calculateHaplotypeProbabilities()

        }

        return currentEmissionProbabilities[haplotype] ?: minProbability
    }


    fun calculateHaplotypeProbabilities() : Map<String, Double> {
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

        val numberOfHaplotypes = haplotypes.size
        val myHapMappings = readMap[myRefRange]

        if (myHapMappings == null || myHapMappings.size == 0) {
            //minProbability has to be set because it will be used for hap index < 0 in getProbObsGivenState()
            //since there are no reads mapped all haplotypes have the same probability
            minProbability = 1.0 / numberOfHaplotypes.toDouble()
            return haplotypes.associateWith { minProbability }

        } else {
            val numberOfReads = myHapMappings.map { it.count }.sum()
            val binom = BinomialDistribution(numberOfReads, pCorrect)
            minProbability = binom.probability(0) //probability of 0 reads hitting a haplotype

            //this next line maps hapid -> sum of counts for sets containing that hapid
            val hapidCountMap = myHapMappings.flatMap { hapCounts -> hapCounts.haplotypeList.map { Pair(it, hapCounts.count)} }
                .groupingBy { it.first }
                .fold(0) { sum, pr -> sum + pr.second }

            return haplotypes.associateWith { hapid ->
                binom.probability(hapidCountMap[hapid] ?: 0)
            }

        }

    }

}

