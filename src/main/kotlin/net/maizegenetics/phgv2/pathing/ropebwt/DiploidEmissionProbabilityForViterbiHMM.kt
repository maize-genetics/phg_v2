package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import org.apache.commons.math3.distribution.BinomialDistribution

class DiploidEmissionProbabilityForViterbiHMM(val readMap: Map<Int, MutableList<Ps4gGameteSet>>,
                                              parentSet: Set<Int>, val probCorrect: Double) {
    val parentList = parentSet.sorted()
    val positionList = readMap.keys.sorted()
    val nParents = parentList.size
    val numberOfPairs = nParents * nParents
    val maxNumberOfTrials = 50
    val lnProbabilityArray = cacheSomeProbabilities()

    fun getLnEmissionProbabilityArray(positionIndex: Int) : DoubleArray {
        //get the list of gamete sets for this position
        val gameteSetList = readMap[positionList[positionIndex]]!!
        val numberOfReads = gameteSetList.sumOf { it.count }

        //get probabilities for all ordered pairs in the parentList
        val probabilityArray  = DoubleArray(numberOfPairs)
        var arrayPtr = 0
        if (numberOfReads > maxNumberOfTrials) {
            val distribution = BinomialDistribution(numberOfReads, probCorrect)
            for (parent1 in parentList) {
                for (parent2 in parentList) {
                    probabilityArray[arrayPtr++] = if (parent1 == parent2) {
                        val indexClassCounts = indexCountsForOneIndex(gameteSetList, parent1)
                        distribution.logProbability(indexClassCounts[0])
                    } else {
                        val indexClassCounts = indexCountsForTwoIndices(gameteSetList, parent1, parent2)
                        distribution.logProbability(indexClassCounts[0])
                    }
                }

            }
        } else {
            for (parent1 in parentList) {
                for (parent2 in parentList) {
                    probabilityArray[arrayPtr++] = if (parent1 == parent2) {
                        val indexClassCounts = indexCountsForOneIndex(gameteSetList, parent1)
                        lnProbabilityArray[indexClassCounts[1]][indexClassCounts[0]]
                    } else {
                        val indexClassCounts = indexCountsForTwoIndices(gameteSetList, parent1, parent2)
                        lnProbabilityArray[indexClassCounts[1]][indexClassCounts[0]]
                    }
                }
            }
        }

        return probabilityArray
    }

    /**
     * For a list of gameteSet counts and an index, returns an int array of (total count for sets containing index, total  count)
     */
    private fun indexCountsForOneIndex(gameteList: MutableList<Ps4gGameteSet>, index: Int): IntArray {
        val indexCount = gameteList.filter{gameteSet -> index in gameteSet.gameteIndices}.sumOf { gameteSet -> gameteSet.count }
        val total = gameteList.sumOf { gameteSet -> gameteSet.count }
        return intArrayOf(indexCount,total)
    }

    private fun indexCountsForTwoIndices(gameteList: MutableList<Ps4gGameteSet>, index1: Int, index2: Int): IntArray {
        val indexCount = gameteList.filter{gameteSet -> index1 in gameteSet.gameteIndices || index2 in gameteSet.gameteIndices}.sumOf { gameteSet -> gameteSet.count }
        val total = gameteList.sumOf { gameteSet -> gameteSet.count }
        return intArrayOf(indexCount,total)
    }

    private fun cacheSomeProbabilities(): Array<DoubleArray> {
        val lnBinomialProbabilities = Array<DoubleArray>(maxNumberOfTrials) {i -> DoubleArray(i + 1)}
        lnBinomialProbabilities[0][0] = 0.0
        for (numberOfTrials in 1 until maxNumberOfTrials) {
            val distribution = BinomialDistribution(numberOfTrials, probCorrect)
            for (numberOfSuccesses in 0 .. numberOfTrials)
                lnBinomialProbabilities[numberOfTrials][numberOfSuccesses] = distribution.logProbability(numberOfSuccesses)
        }
        return lnBinomialProbabilities
    }


}