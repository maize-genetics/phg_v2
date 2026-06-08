package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import org.apache.commons.math3.distribution.BinomialDistribution
import kotlin.collections.associateWith
import kotlin.collections.isNullOrEmpty
import kotlin.math.ln

class HaploidPS4GEmissionProbability(val readMap: Map<Int, MutableList<Ps4gGameteSet>>, val genomeIndexSet: Set<Int>, val pCorrect: Double) {
    private var currentPosition = -1
    private var currentEmissionProbabilities = mapOf<Int, Double>()
    private var nullProbability : Double = -10.0
    private val minProbability = Double.MIN_VALUE

    fun getLnProbObsGivenState(genomeIndex: Int, position: Int): Double {
        //refrange is a 0..n index into nodeTree.keys
        //haplotype is an index into the List<HaplotypeNode>> for this ReferenceRange
        if (currentPosition != position) {
            currentPosition = position
            currentEmissionProbabilities = calculateLnHaplotypeProbabilities()

        }

        return currentEmissionProbabilities[genomeIndex] ?: nullProbability
    }

    private fun calculateLnHaplotypeProbabilities() : Map<Int, Double> {
        /*
        * P_emission = P(obs|haplotype)
        * = binom(#trials, #successes, P_success)
        * where #trials = number of reads
        *   #successes = number of reads hitting the haplotype
        *   P_success = pcorrect
        *
        * assign minProbability for haplotype index = -1
        * */

        val myPosition = currentPosition
        val myReadMappings = readMap[myPosition]

        if (myReadMappings.isNullOrEmpty()) {
            //if there no reads then all haplotypes get assigned the same probability, it does not matter what it is.
            //note that 0.0 is ln(1), so this actually assigns a probability of 1 not zero. One interpretation is that
            //if there are no reads then every haplotype (including null haplotypes) will have 0 reads with a probability of 1.
            nullProbability = 0.0
            return genomeIndexSet.associateWith { 0.0 }

        } else {
            val numberOfReads = myReadMappings.map { it.count }.sum()
            val binom = BinomialDistribution(numberOfReads, pCorrect)

            //this next line maps hapid -> sum of counts for sets containing that hapid
            val indexCountMap = myReadMappings.flatMap { indexCounts -> indexCounts.gameteIndices.map { Pair(it, indexCounts.count)} }
                .groupingBy { it.first }
                .fold(0) { sum, pr -> sum + pr.second }

            //since a null haplotype should have 0 counts, assign it min
            nullProbability = ln( binom.probability( 0).coerceAtLeast(minProbability) )
            return genomeIndexSet.associateWith { index ->
                val prob = binom.probability(indexCountMap[index] ?: 0).coerceAtLeast(minProbability)
                ln( prob )
            }

        }

    }
}