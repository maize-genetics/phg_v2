package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet

class DiploidPS4GEmissionProbability(val readMap: Map<Int, MutableList<Ps4gGameteSet>>, val genomeIndexSet: Set<Int>, val pCorrect: Double {
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

    }
}