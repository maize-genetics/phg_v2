package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import org.apache.commons.math3.distribution.BinomialDistribution
import org.apache.commons.math3.special.Gamma.logGamma
import kotlin.math.ln

class DiploidPS4GEmissionProbability(val readMap: Map<Int, MutableList<Ps4gGameteSet>>,
                                     val genomeIndexSet: Set<Int>,
                                     val pCorrect: Double) {
    private val pBoth = pCorrect * 0.5  //make this user settable
    private val pNeither = 1 - pCorrect
    private var currentPosition = -1
    private var currentEmissionProbabilities : Array<DoubleArray>
    private val minProbability = Double.MIN_VALUE
    private val indexClassProbabilities: DoubleArray

    init {
        val arraySize = genomeIndexSet.max() + 1
        currentEmissionProbabilities = Array(arraySize) { DoubleArray(arraySize) }
        val pOnlyOne = (1 - pBoth - pNeither)/2
        indexClassProbabilities = doubleArrayOf(pOnlyOne, pOnlyOne, pBoth, pNeither)
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
            for (ndx1 in 0 until indexSize) {
                //when both indices are the same (homozygous locus) use the binomial distribution
                val indexClassCounts = indexCountsForOneIndex(gameteSetList, ndx1)
                currentEmissionProbabilities[ndx1][ndx1] = BinomialDistribution(indexClassCounts[1], pCorrect)
                    .logProbability(indexClassCounts[0])

                for (ndx2 in ndx1 + 1 until indexSize) {
                    val indexClassCounts = indexCountsForTwoIndices(gameteSetList, ndx1, ndx2)
                    val prob = lnMultinomialProbability(indexClassCounts, indexClassProbabilities)
                    currentEmissionProbabilities[ndx1][ndx2] = prob
                    currentEmissionProbabilities[ndx2][ndx1] = prob
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
        val both = gameteList.filter { gameteSet ->  index1 in gameteSet.gameteIndices && index2 in gameteSet.gameteIndices }
            .sumOf { gameteSet -> gameteSet.count }
        val neither = gameteList.filter { gameteSet ->  index1 !in gameteSet.gameteIndices && index2 !in gameteSet.gameteIndices }
            .sumOf { gameteSet -> gameteSet.count }

        return intArrayOf(only1, only2, both, neither)
    }

    /**
     * For a list of gameteSet counts and an index, returns an int array of (total count for sets containing index, total  count)
     */
    fun indexCountsForOneIndex(gameteList: MutableList<Ps4gGameteSet>, index: Int): IntArray {
        val index1Count = gameteList.filter{gameteSet -> index in gameteSet.gameteIndices}.sumOf { gameteSet -> gameteSet.count }
        val total = gameteList.sumOf { gameteSet -> gameteSet.count }
        return intArrayOf(index1Count, total)
    }

    companion object {
        fun lnMultinomialProbability(counts: IntArray, probabilities: DoubleArray): Double {
            if (counts.size != probabilities.size) throw java.lang.IllegalArgumentException("multinomialProbability error: counts and probabilities arrays do not have the same size.")
            val totalCount = counts.sum()

            val logprod = counts.indices.sumOf { counts[it] * ln(probabilities[it]) }
            val numerator = logFactorial(totalCount)
            val denominator = counts.sumOf { logFactorial(it) }
            return numerator - denominator + logprod
        }

        val smallFactorials = DoubleArray(1000) {ln(it + 1.0)}.runningFold(0.0) {sum, element -> sum + element}

        /**
         *  Calculates the (natural) log factorial of any positive integer using the exact value for 1 to 1000 and
         *  logGamma(n+1) otherwise. This method produces values identical to the R method dmultinom.
         *  @param [intval] an integer
         *  @return   the natural log of the factorial of intval
         */
        fun logFactorial(intval: Int): Double {

            return if (intval <= 1000) smallFactorials[intval] else {
                val n = intval.toDouble()
                logGamma(n + 1.0)
            }
        }
    }
}