package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import org.apache.commons.math3.distribution.BinomialDistribution

class EmissionProbabilityForForwardBackward(readMap: Map<Int, MutableList<Ps4gGameteSet>>,
                                            parentSet: Set<Int>, probCorrect: Double)
    : AbstractEmissionProbability(readMap, parentSet, probCorrect) {

    override fun binomialProbability(
        distribution: BinomialDistribution,
        numberOfSuccesses: Int
    ): Double {
        return distribution.probability(numberOfSuccesses)
    }

}