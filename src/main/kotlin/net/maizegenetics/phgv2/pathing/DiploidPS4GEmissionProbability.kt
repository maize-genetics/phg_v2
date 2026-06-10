package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import org.apache.commons.math3.distribution.BinomialDistribution

class DiploidPS4GEmissionProbability(val readMap: Map<Int, MutableList<Ps4gGameteSet>>,
                                     val genomeIndexSet: Set<Int>,
                                     val pCorrect: Double) {
    private var currentPosition = -1
    private var currentEmissionProbabilities : Array<DoubleArray>
    private val minProbability = Double.MIN_VALUE

    init {
        val arraySize = genomeIndexSet.max() + 1
        currentEmissionProbabilities = Array(arraySize) { DoubleArray(arraySize) }
    }

    fun getLnProbObsGivenState(genomeIndex1: Int, genomeIndex2: Int, position: Int): Double {
        //refrange is a 0..n index into nodeTree.keys
        //haplotype is an index into the List<HaplotypeNode>> for this ReferenceRange
        if (currentPosition != position) {
            currentPosition = position
            calculateLnDiploidProbabilities()
        }

        return currentEmissionProbabilities[genomeIndex1][genomeIndex2]
    }

    fun calculateLnDiploidProbabilities() {
        val gameteSetList = readMap[currentPosition]

        if (gameteSetList == null) {
            //if there are no counts set all probabilities equal to 0.0.
            for (subArray in currentEmissionProbabilities) { subArray.fill(0.0)}
        } else {
            //for each genome index pair:
            //get counts of index1 only, index2 only, both, neither
            val indexSize = currentEmissionProbabilities.size
            for (ndx1 in 1..indexSize) {
                //when both indices are the same (homozygous locus) use the binomial distribution
                val indexClassCounts = indexCountsForOneIndex(gameteSetList, ndx1)
                currentEmissionProbabilities[ndx1][ndx1] = BinomialDistribution(indexClassCounts[1], pCorrect)
                    .logProbability(indexClassCounts[0])

                for (ndx2 in ndx1 + 1..indexSize) {
                    val indexClassCounts = indexCountsForTwoIndices(gameteSetList, ndx1, ndx2)

                }
            }
        }


    }

    /**
     * For a list of gameteSet counts and two indices, returns an array containg the following:
     * index1 only - the total count of gameteSets containing index1 but not index2,
     * index2 only - the total count of gameteSets containing index2 but not index1,
     * both - the total count of gameteSets containing both index1 and index2,
     * neither - the total count of gameteSets containing neither index nor index2, and
     * total - the sum of all gameteSet counts
     */
    fun indexCountsForTwoIndices(gameteList: List<Ps4gGameteSet>, index1: Int, index2: Int) : IntArray {
        val only1 = gameteList.filter { gameteSet ->  index1 in gameteSet.gameteIndices && index2 !in gameteSet.gameteIndices }
            .sumOf { gameteSet -> gameteSet.count }
        val only2 = gameteList.filter { gameteSet ->  index1 !in gameteSet.gameteIndices && index2 in gameteSet.gameteIndices }
            .sumOf { gameteSet -> gameteSet.count }
        val both = gameteList.filter { gameteSet ->  index1 !in gameteSet.gameteIndices && index2 in gameteSet.gameteIndices }
            .sumOf { gameteSet -> gameteSet.count }
        val neither = gameteList.filter { gameteSet ->  index1 !in gameteSet.gameteIndices && index2 in gameteSet.gameteIndices }
            .sumOf { gameteSet -> gameteSet.count }
        val total = only1 + only2 + both + neither
        return intArrayOf(only1, only2, both, neither, total)
    }

    /**
     * For a list of gameteSet counts and an index, returns an int array of (total count for sets containing index, total  count)
     */
    fun indexCountsForOneIndex(gameteList: MutableList<Ps4gGameteSet>, index: Int): IntArray {
        val index1Count = gameteList.filter{gameteSet -> index in gameteSet.gameteIndices}.sumOf { gameteSet -> gameteSet.count }
        val total = gameteList.sumOf { gameteSet -> gameteSet.count }
        return intArrayOf(index1Count, total)
    }
}