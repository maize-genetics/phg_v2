package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.pathing.ropebwt.GameteSetEmissionProbability
import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import net.maizegenetics.phgv2.pathing.ropebwt.ViterbiHMM
import org.apache.logging.log4j.LogManager
import kotlin.math.expm1
import kotlin.math.ln

/**
 * Finds a haploid or diploid path through a [HaplotypeGraph] from read mappings, using the Viterbi HMM
 * of `impute-path-from-ps4g` ([ViterbiHMM] with [GameteSetEmissionProbability]) in place of
 * [PathFinderWithViterbiHMM]. It takes the same read mappings and returns the same [PathFinderWithViterbiHMM.PathNode]
 * list, so the two are interchangeable in [FindPaths].
 *
 * Positions are the reference ranges of a contig that have reads. Each read's haplotype ID set becomes
 * the set of sample gametes carrying one of those haplotypes in the range, so founders that share a
 * haplotype are hit together, which is what the gameteset emission expects.
 *
 * Transitions are constant per reference range: genic and intergenic haplotypes divide the genome so
 * that crossover density tracks range density. Ranges without reads still count, so the step between
 * two ranges with reads k ranges apart has, for each haplotype, a probability [sameGameteProbability]^k
 * of staying and (1 - [sameGameteProbability]^k) / (n - 1) of switching to any particular other founder.
 *
 * A diploid path is the inbreeding coefficient 0 recursion, with the two haplotypes switching
 * independently; a haploid path is the inbreeding coefficient 1 recursion, where only homozygous
 * states are reachable. A founder with no haplotype in a range is treated as absent there, so a
 * heterozygous state with that founder scores the same as the homozygous state of the other.
 *
 * Every range of a contig with reads is assigned a state. A range without reads takes the state of
 * its neighbors with reads; between two that differ, the ranges up to the midpoint take the
 * earlier state and the rest the later one. A contig with no reads gets no path.
 */
class PathFinderWithFastHMM(
    val graph: HaplotypeGraph,
    val isHaploidPath: Boolean,
    val probCorrect: Double,
    val sameGameteProbability: Double,
    useLikelyParents: Boolean = false,
    val maxParents: Int = Int.MAX_VALUE,
    val minCoverage: Double = 1.0
) {
    private val myLogger = LogManager.getLogger(PathFinderWithFastHMM::class.java)

    private val gametes = graph.sampleGametesInGraph().toList()
    private val gameteToIndex = gametes.withIndex().associate { (index, gamete) -> gamete to index }
    private val rangesByContig = graph.rangesByContig()
    private val findLikelyParents = useLikelyParents && (maxParents < gametes.size || minCoverage < 1.0)
    private val lnSame = ln(sameGameteProbability)

    init {
        require(sameGameteProbability > 0.0 && sameGameteProbability <= 1.0) {
            "sameGameteProbability must be in (0, 1], not $sameGameteProbability"
        }
    }

    /**
     * Returns the path through every contig of the graph, in contig and range order, and the parents
     * chosen when likely parents are used (otherwise an empty list).
     */
    fun findBestPath(readMap: Map<ReferenceRange, Map<List<String>, Int>>):
            Pair<List<PathFinderWithViterbiHMM.PathNode>, List<MostLikelyParents.ParentStats>> {
        val likelyParentList = if (findLikelyParents) MostLikelyParents(graph)
            .findMostLikelyParents(readMap, maxParents = maxParents, minCoverage = minCoverage)
        else listOf()
        val parentSet = if (findLikelyParents) likelyParentList.map { gameteToIndex[it.parent]!! }.toSet()
        else gametes.indices.toSet()

        val nodeList = mutableListOf<PathFinderWithViterbiHMM.PathNode>()
        for (contig in graph.contigs) {
            val startTime = System.nanoTime()
            nodeList.addAll(contigPath(contig, readMap, parentSet))
            myLogger.info("elapsed time for $contig was ${(System.nanoTime() - startTime) / 1e9} sec")
        }
        return Pair(nodeList, likelyParentList)
    }

    private fun contigPath(
        contig: String,
        readMap: Map<ReferenceRange, Map<List<String>, Int>>,
        parentSet: Set<Int>
    ): List<PathFinderWithViterbiHMM.PathNode> {
        val ranges = rangesByContig[contig] ?: return emptyList()
        val contigReads = gameteSetsByRange(ranges, readMap)
        if (contigReads.isEmpty()) {
            myLogger.warn("No reads on $contig; no path written for it")
            return emptyList()
        }

        val positions = contigReads.keys.sorted()
        val parentList = parentSet.sorted()
        val nParents = parentList.size

        val absence = GameteSetEmissionProbability.FounderAbsence { position, absent ->
            val present = graph.sampleGameteToHaplotypeId(ranges[position])
            for (i in parentList.indices) absent[i] = gametes[parentList[i]] !in present
        }
        val emission = GameteSetEmissionProbability(contigReads, parentSet, probCorrect, founderAbsence = absence)

        val (lnNoSwitch, lnSwitch) = rangeTransitions(positions, nParents)
        val (states, lnProbability) = ViterbiHMM(if (isHaploidPath) 1.0 else 0.0, sameGameteProbability, probCorrect)
            .findDiploidStatePath(
                nParents, positions.size,
                emission::getDiploidEmissionProbabilityArray,
                emission::getHomozygousEmissionProbabilityArray,
                lnNoSwitch, lnSwitch
            )

        // the state of each range, filling ranges without reads from their neighbors
        val rangeStates = IntArray(ranges.size)
        var next = 0
        for (rangeIndex in ranges.indices) {
            while (next < positions.size && positions[next] < rangeIndex) next++
            rangeStates[rangeIndex] = when {
                next == 0 -> states[0]
                next == positions.size -> states.last()
                positions[next] == rangeIndex -> states[next]
                rangeIndex <= (positions[next - 1] + positions[next]) / 2 -> states[next - 1]
                else -> states[next]
            }
        }

        return ranges.indices.map { rangeIndex ->
            val state = rangeStates[rangeIndex]
            val first = gametes[parentList[state / nParents]]
            val second = gametes[parentList[state % nParents]]
            val sampleGametes = if (isHaploidPath) listOf(first) else listOf(first, second)
            PathFinderWithViterbiHMM.PathNode(ranges[rangeIndex], null, sampleGametes, lnProbability)
        }
    }

    /**
     * Converts the read mappings of a contig's ranges to gamete sets, keyed by the range's index in
     * [ranges]. Ranges without reads are left out.
     */
    private fun gameteSetsByRange(
        ranges: List<ReferenceRange>,
        readMap: Map<ReferenceRange, Map<List<String>, Int>>
    ): Map<Int, MutableList<Ps4gGameteSet>> {
        val result = mutableMapOf<Int, MutableList<Ps4gGameteSet>>()
        ranges.forEachIndexed { rangeIndex, range ->
            val counts = readMap[range]
            if (counts.isNullOrEmpty() || counts.values.sum() == 0) return@forEachIndexed
            val hapIdToGametes = graph.hapIdToSampleGametes(range)
            result[rangeIndex] = counts.map { (hapIds, count) ->
                val indices = hapIds.flatMap { hapIdToGametes[it] ?: emptyList<SampleGamete>() }
                    .map { gameteToIndex[it]!! }.distinct().sorted().toIntArray()
                Ps4gGameteSet(indices, count)
            }.toMutableList()
        }
        return result
    }

    /**
     * Per-haplotype log transition probabilities into each position: k ranges after the previous
     * position, ln(p^k) to stay and ln((1 - p^k) / (n - 1)) to switch to a particular other founder.
     * Position 0 has no predecessor and its values are unused.
     */
    private fun rangeTransitions(positions: List<Int>, nParents: Int): Pair<DoubleArray, DoubleArray> {
        val lnNoSwitch = DoubleArray(positions.size)
        val lnSwitch = DoubleArray(positions.size)
        for (index in 1 until positions.size) {
            val lnStay = (positions[index] - positions[index - 1]) * lnSame
            val switch = -expm1(lnStay)
            lnNoSwitch[index] = lnStay
            lnSwitch[index] = if (nParents > 1 && switch > 0.0) ln(switch / (nParents - 1)) else -1.0e6
        }
        return Pair(lnNoSwitch, lnSwitch)
    }
}
