package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import org.apache.commons.math3.distribution.BinomialDistribution

class DiploidEmissionProbabilityForLikelyParents(val readMap: Map<Int, MutableList<Ps4gGameteSet>>,
                                                 val genomeIndexSet: Set<Int>,
                                                 val likelyParentsList: List<Int>,
                                                 val pCorrect: Double) {

    private val pBoth = pCorrect * 0.5  //make this user settable
    private val pNeither = 1 - pCorrect
    private var currentPosition = -1
    private val indexClassProbabilities: DoubleArray
    private val probabilityArray: DoubleArray
    private val nGenomes = genomeIndexSet.size
    private val nParents = likelyParentsList.size

    init {
        val pOnlyOne = (1 - pBoth - pNeither)/2
        indexClassProbabilities = doubleArrayOf(pOnlyOne, pOnlyOne, pBoth, pNeither)

        val arraySize = (nGenomes * (nGenomes + 1)) / 2
        probabilityArray = DoubleArray(arraySize)
    }

    fun parentIndicesToArrayIndex(index1: Int, index2: Int): Int {
        //order of genome indices does not matter: ndx1,ndx2 and ndx2,ndx1 should give the same array index
        //do this by always having index1 <= index2
        var orderedIndex1 = index1
        var orderedIndex2 = index2
        if (index1 > index2) {
            orderedIndex1 = index2
            orderedIndex2 = index1
        }

        //index = i(2n - i + 1)/2 + j - i
        //where i <= j (orderedIndex1 = i, orderedIndex2 = j)
        val index = orderedIndex1 * (2 * nGenomes - orderedIndex1 + 1) / 2 + orderedIndex2 - orderedIndex1
        return index
    }

    fun getLnProbObsGivenState(genomeIndex1: Int, genomeIndex2: Int, position: Int): Double {
        if (currentPosition != position) {
            currentPosition = position
            calculateLnDiploidProbabilities()
        }

        val arrayIndex = parentIndicesToArrayIndex(genomeIndex1, genomeIndex2)
        return probabilityArray[arrayIndex]
    }

    fun calculateLnDiploidProbabilities() {
        val gameteSetList = readMap[currentPosition]

        if (gameteSetList == null) {
            //if there are no counts set all probabilities equal to 0.0.
            probabilityArray.fill(0.0)
        } else {

            for (ndx1 in 0 until nParents) {
                val parent1 = likelyParentsList[ndx1]
                //when both indices are the same (homozygous locus) use the binomial distribution
                val indexClassCounts = indexCountsForOneIndex(gameteSetList, parent1)
                val arrayIndex = parentIndicesToArrayIndex(parent1, parent1)
                probabilityArray[arrayIndex] = BinomialDistribution(indexClassCounts[1], pCorrect)
                    .logProbability(indexClassCounts[0])

                for (ndx2 in ndx1 + 1 until nParents) {
                    val parent2 = likelyParentsList[ndx2]
                    val indexClassCounts = indexCountsForTwoIndices(gameteSetList, parent1, parent2)
                    val arrayIndex = parentIndicesToArrayIndex(parent1, parent2)
                    probabilityArray[arrayIndex] = BinomialDistribution(indexClassCounts[1], pCorrect)
                        .logProbability(indexClassCounts[0])
                }
            }
        }
    }

    /**
     * For a list of gameteSet counts and an index, returns an int array of (total count for sets containing index, total  count)
     */
    fun indexCountsForOneIndex(gameteList: MutableList<Ps4gGameteSet>, index: Int): IntArray {
        val indexCount = gameteList.filter{gameteSet -> index in gameteSet.gameteIndices}.sumOf { gameteSet -> gameteSet.count }
        val total = gameteList.sumOf { gameteSet -> gameteSet.count }
        return intArrayOf(indexCount,total)
    }

    fun indexCountsForTwoIndices(gameteList: MutableList<Ps4gGameteSet>, index1: Int, index2: Int): IntArray {
        val indexCount = gameteList.filter{gameteSet -> index1 in gameteSet.gameteIndices || index2 in gameteSet.gameteIndices}.sumOf { gameteSet -> gameteSet.count }
        val total = gameteList.sumOf { gameteSet -> gameteSet.count }
        return intArrayOf(indexCount,total)
    }
}