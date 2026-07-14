package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import org.apache.commons.math3.distribution.BinomialDistribution

/**
 * Calculates emission probabilities based on observations stored in the [readMap]. It has methods that return an array
 * of probabilities for a position. The values in the array are indexed a sorted list of parents in the [parentSet].
 * [probCorrect] is the probability that a read coming from parent A maps to parent A.
 */
class EmissionProbabilityForViterbiHMM(val readMap: Map<Int, MutableList<Ps4gGameteSet>>,
                                       parentSet: Set<Int>, val probCorrect: Double) {
    val parentList = parentSet.sorted()
    val positionList = readMap.keys.sorted()
    val nParents = parentList.size
    val numberOfPairs = nParents * nParents
    val maxNumberOfTrials = 50
    val lnProbabilityArray = cacheSomeProbabilities()

    /**
     * Returns a double array of the natural log of emission probabilities at the input [positionIndex].
     * The length of the array is the number of parents in the parent set used to create the class. The order within the
     * array corresponds to a sorted list of the original indices of the parents.
     */
    fun getLnHaploidEmissionProbabilityArray(positionIndex: Int) : DoubleArray {
        //get the list of gamete sets for this position
        val gameteSetList = readMap[positionList[positionIndex]]!!
        val numberOfReads = gameteSetList.sumOf { it.count }

        //get probabilities for the parentList
        val probabilityArray  = DoubleArray(nParents)
        var arrayPtr = 0
        if (numberOfReads > maxNumberOfTrials) {
            val distribution = BinomialDistribution(numberOfReads, probCorrect)
            for (parentIndex in parentList) {
                val indexClassCounts = indexCountsForOneIndex(gameteSetList, parentIndex)
                probabilityArray[arrayPtr++] = distribution.logProbability(indexClassCounts)
            }
        } else {
            for (parentIndex in parentList) {
                val indexClassCounts = indexCountsForOneIndex(gameteSetList, parentIndex)
                probabilityArray[arrayPtr++] = lnProbabilityArray[numberOfReads][indexClassCounts]
            }
        }

        return probabilityArray
    }

    /**
     * Returns a double array of the natural log of emission probabilities at the input [positionIndex].
     * The length of the array is the square of the number of parents in the parent set used to create the class. The
     * resulting array is indexed on ordered pairs of the parent indices.
     */
    fun getLnDiploidEmissionProbabilityArray(positionIndex: Int) : DoubleArray {
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
                        distribution.logProbability(indexClassCounts)
                    } else {
                        val indexClassCounts = indexCountsForTwoIndices(gameteSetList, parent1, parent2)
                        distribution.logProbability(indexClassCounts)
                    }
                }

            }
        } else {
            for (parent1 in parentList) {
                for (parent2 in parentList) {
                    probabilityArray[arrayPtr++] = if (parent1 == parent2) {
                        val indexClassCounts = indexCountsForOneIndex(gameteSetList, parent1)
                        lnProbabilityArray[numberOfReads][indexClassCounts]
                    } else {
                        val indexClassCounts = indexCountsForTwoIndices(gameteSetList, parent1, parent2)
                        lnProbabilityArray[numberOfReads][indexClassCounts]
                    }
                }
            }
        }

        return probabilityArray
    }

    /**
     * For a list of gameteSet counts and an index, returns an int array of (total count for sets containing index, total  count)
     */
    private fun indexCountsForOneIndex(gameteList: MutableList<Ps4gGameteSet>, index: Int): Int {
        val indexCount = gameteList.filter{gameteSet -> index in gameteSet.gameteIndices}.sumOf { gameteSet -> gameteSet.count }
        return indexCount
    }

    private fun indexCountsForTwoIndices(gameteList: MutableList<Ps4gGameteSet>, index1: Int, index2: Int): Int {
        val indexCount = gameteList.filter{gameteSet -> index1 in gameteSet.gameteIndices || index2 in gameteSet.gameteIndices}.sumOf { gameteSet -> gameteSet.count }
        return indexCount
    }

    private fun cacheSomeProbabilities(): Array<DoubleArray> {
        val lnBinomialProbabilities = Array(maxNumberOfTrials + 1) {i -> DoubleArray(i + 1)}
        lnBinomialProbabilities[0][0] = 0.0
        for (numberOfTrials in 1 .. maxNumberOfTrials) {
            val distribution = BinomialDistribution(numberOfTrials, probCorrect)
            for (numberOfSuccesses in 0 .. numberOfTrials)
                lnBinomialProbabilities[numberOfTrials][numberOfSuccesses] = distribution.logProbability(numberOfSuccesses)
        }
        return lnBinomialProbabilities
    }


}