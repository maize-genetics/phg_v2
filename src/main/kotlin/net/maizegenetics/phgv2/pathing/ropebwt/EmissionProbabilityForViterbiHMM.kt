package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import org.apache.commons.math3.distribution.BinomialDistribution

/**
 * Calculates emission probabilities based on observations stored in the [readMap]. It has methods that return an array
 * of probabilities for a position. The values in the array are indexed a sorted list of parents in the [parentSet].
 * [probCorrect] is the probability that a read coming from parent A maps to parent A.
 */
class EmissionProbabilityForViterbiHMM(readMap: Map<Int, MutableList<Ps4gGameteSet>>,
                                       parentSet: Set<Int>, probCorrect: Double)
    : AbstractEmissionProbability(readMap, parentSet, probCorrect) {


    override fun binomialProbability(
        distribution: BinomialDistribution,
        numberOfSuccesses: Int
    ): Double {
        return distribution.logProbability(numberOfSuccesses)
    }


}