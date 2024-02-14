package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import org.apache.logging.log4j.LogManager
import kotlin.math.E
import kotlin.math.ln
import kotlin.math.log

/**
 * This class finds the most likely path through a graph from read mapping data for any number of samples using the Viterbi algorithm.
 *
 * @param graph The [HaplotypeGraph] used to impute paths. The reads must have been mapped against this graph (or a superset of it).
 * @param probCorrect The probability that a read was mapped to the correct range.
 * @param sameGameteProbability The probability that a parent in any given range is the same as the parent in the previous range.
 * @param minGametesPerRange Ranges with fewer than this many sampleGametes in the graph will not be imputed.
 * @param minReadsPerRange Ranges with fewer than this many mapped reads will not be imputed.
 * @param maxReadsPerKB If the number of reads for a range is greater than maxReadsPerKB, then those reads will not be used.
 * @param inbreedCoef The average inbreeding coefficient (F) of the submitted samples. Only used for diploid path finding.
 * @param useLikelyParents Use only likely parents for path finding. Likely parents are determined for each sample based on the read mappings for that sample.
 * @param maxParents If useLikelyParents is true, use at most this many parents
 * @param minCoverage If useLikelyParents is true, stop selecting parents when combined coverage is greater than or equal to this value.
 *
 * Likely parents (or ancestors) are determined and used in pathfinding only when (1) useLikelyParents is true and (2) either
 * maxParents is less than the number of parents in the graph or minCoverage is less than 1.0.
 */
class PathFinderWithViterbiHMM(
        val graph: HaplotypeGraph,
        val isHaploidPath: Boolean,
        val probCorrect: Double,
        val sameGameteProbability: Double,
        val minGametesPerRange: Int,
        val minReadsPerRange: Int,
        val maxReadsPerKB: Int,
        val inbreedCoef: Double = 1.0,
        useLikelyParents: Boolean = false,
        val maxParents: Int = Int.MAX_VALUE,
        val minCoverage: Double = 1.0
) {

    private val myLogger = LogManager.getLogger(PathFinderWithViterbiHMM::class.java)
    private val findLikelyParents: Boolean
    private val parentFinder: MostLikelyParents?
    private val sampleGametesInGraph = graph.sampleGametesInGraph()
    private val rangeToHaplotypeMap = graph.refRangeToHapIdList()

    data class PathNode(
        val refRange: ReferenceRange,
        val parent: PathNode?,
        val sampleGametes: List<SampleGamete>,
        val totalProbability: Double
    )

    init {
        findLikelyParents = useLikelyParents && (maxParents < graph.sampleGametesInGraph().size || minCoverage < 1.0)
        parentFinder = if (findLikelyParents) MostLikelyParents(graph) else null
    }

    fun findBestPath(readMap: Map<ReferenceRange, Map<List<String>, Int>>): Pair<List<PathNode>, List<MostLikelyParents.ParentStats>> {
        //find likely parents, if required
        val likelyParentList = if (findLikelyParents) parentFinder!!
            .findMostLikelyParents(readMap, maxParents = maxParents, minCoverage = minCoverage)
        else listOf()

        val nodeList: MutableList<PathNode> = mutableListOf()
        graph.contigs.forEach { chr ->
            val finalNode = when {
                isHaploidPath && findLikelyParents -> haploidViterbi(chr, readMap, likelyParentList.map { it.parent }.toSet())
                isHaploidPath && !findLikelyParents -> haploidViterbi(chr, readMap)
                !isHaploidPath && findLikelyParents -> diploidViterbi(chr, readMap, likelyParentList.map { it.parent }.toSet())
                else -> diploidViterbi(chr, readMap)
            }

            if (finalNode != null) {
                val reverseNodeList = mutableListOf(finalNode)
                var currentNode = finalNode!!
                while (currentNode.parent != null) {
                    currentNode = currentNode.parent!!
                    reverseNodeList.add(currentNode)
                }
                nodeList.addAll(reverseNodeList.reversed())
            }
        }
        return Pair(nodeList, likelyParentList)

    }
    /**
     * This method takes a graph and  a single chromosome and read mappings and returns a list of HaplotypeNodes that
     * represents the most likely path through the graph for that chromosome given the read mappings. The graph must have been constructed so that
     * there is only a single taxon for each node and all nodes have the same taxa. This condition is not checked because
     * the same graph may be used to process multiple samples. It is the responsibility of the calling code to ensure the condition
     * is satisfied.
     *
     * The input parameter [gameteSet], which defaults to all the SampleGametes in the graph, can be used
     * to restrict the imputation to a subset of potential ancestors. Only the SampleGametes in gameteList will be
     * used for imputation.
     *
     * @param chrom a chromosome name
     * @param readMap   a map of read mappings by ReferenceRange for graph.
     * @param gameteSet a list of the SampleGametes used for imputation
     * @return  a list of HaplotypeNodes representing the most likely path
     */
    private fun haploidViterbi(
        chrom: String,
        readMap: Map<ReferenceRange, Map<List<String>, Int>>,
        gameteSet: Set<SampleGamete> = sampleGametesInGraph
    ): PathNode? {
        myLogger.info("Finding path for chromosome $chrom using haploidViterbi")
        val switchProbability = 1 - sameGameteProbability
        val numberOfSampleGametes = sampleGametesInGraph.size
        val chrRanges = graph.ranges().filter { it.contig == chrom }
        val logSwitch = log(switchProbability / (numberOfSampleGametes.toDouble() - 1.0), E)
        val logNoSwitch = log(1.0 - switchProbability, E)

        //diagnostic counters for discarded ranges:
        //countTooFewReads, countTooManyReadsPerKB, countOfTooFewGametes, countDiscardedRanges
        val counters = IntArray(4) { 0 }

        //create emission probability
        val emissionProb = HaplotypeEmissionProbability(rangeToHaplotypeMap, readMap, probCorrect)

        val rangeIterator = chrRanges.iterator()
        var initialRange = rangeIterator.next()


        //logic for use next range
        //get value for useRange(initialRange)
        //if false go to next range and increment range index
        //this way rangeIndex is the index of the selected range
        //if useRange is false for the last range, then return null for pathNode
        while (!useRange(readMap[initialRange],counters,initialRange)) {
            if(rangeIterator.hasNext()) {
                initialRange = rangeIterator.next()
            } else {
                myLogger.debug("Finished processing reads for chromosome $chrom, sample: ${counters[3]} ranges discarded out of ${chrRanges.size}.")
                myLogger.debug("${counters[0]} ranges had too few reads; ${counters[1]} ranges had too many reads; ${counters[2]} ranges too few samples with haplotypes")
                return null
            }
        }

        var paths = ArrayList<PathNode>()

        //make the initial PathNode for each node in the first range
        //the path nodes keep track of all the nodes on a path by holding a link to the parent node of each node
        //first need a map of sampleGamete -> haplotype


        var sampleGameteToHaplotypeMap = graph.hapIdToSampleGametes(initialRange).map { (hap, gameteList) ->
            gameteList.map { Pair(it, hap) }
        }.flatten().toMap()

        for (gamete in sampleGametesInGraph) {
            val myHaplotype = sampleGameteToHaplotypeMap[gamete]
            paths.add(PathNode(initialRange, null, listOf(gamete), emissionProb.getLnProbObsGivenState(myHaplotype, initialRange)))
        }

        //for each reference range update paths with probabilities and new nodes
        while (rangeIterator.hasNext()) {
            //for each node in the next range find the maximum value of path probability * transition
            //actually only two to consider (1-r) * same taxon (switchProbability) and r/(n-1) * all other taxa
            //where n is number of taxa. Since the transition likelihood is the same for all n of the recombinant taxa,
            //only the most likely of those needs to be considered.
            val nextRange = rangeIterator.next()

            sampleGameteToHaplotypeMap = graph.hapIdToSampleGametes(nextRange).map { (hap, gameteList) ->
                gameteList.map { Pair(it, hap) }
            }.flatten().toMap()

            if (!useRange(readMap[nextRange], counters, nextRange)) {
                continue
            }
            val newPaths = ArrayList<PathNode>()

            //TODO revisit whether randomly picking best is a good strategy
            //choose the most probable path from the previous range. If more than one, any one will do.
            //the most likely path in the previous range
            val bestPath = paths.maxByOrNull { it.totalProbability }
            check(bestPath != null) { "no most likely path before range at ${nextRange.contig}:${nextRange.start}" }

            //map of sampleGamete to path node for the previous range
            val gameteToPath = paths.associateBy { it.sampleGametes[0] }

            //iterate over gametes for the new range
            val probSwitch = bestPath.totalProbability + logSwitch
            gameteSet.forEach { sampleGamete ->
                val gameteHaplotype = sampleGameteToHaplotypeMap[sampleGamete]
                if (bestPath.sampleGametes[0] == sampleGamete) {
                    //this is the best path for this node since is also the no switch path (same gamete)
                    newPaths.add(
                        PathNode(
                            nextRange,
                            bestPath,
                            listOf(sampleGamete),
                            bestPath.totalProbability + logNoSwitch +
                                    emissionProb.getLnProbObsGivenState(gameteHaplotype, nextRange)
                        )
                    )
                } else {
                    //since the best path is a switch path (switching gametes), must compare it to the no switch path
                    val samePath = gameteToPath[sampleGamete]
                    val probNoSwitch = samePath!!.totalProbability + logNoSwitch
                    if (probSwitch > probNoSwitch) newPaths.add(
                        PathNode(
                            nextRange,
                            bestPath,
                            listOf(sampleGamete),
                            probSwitch + emissionProb.getLnProbObsGivenState(gameteHaplotype, nextRange)
                        )
                    )
                    else newPaths.add(
                        PathNode(
                            nextRange,
                            samePath,
                            listOf(sampleGamete),
                            probNoSwitch + emissionProb.getLnProbObsGivenState(gameteHaplotype, nextRange)
                        )
                    )
                }
            }

            paths = newPaths
        }

        //countTooFewReads, countTooManyReadsPerKB, countOfTooFewGametes, countDiscardedRanges
        myLogger.info("Finished processing reads for a sample: ${counters[3]} ranges discarded out of ${chrRanges.size}.")
        myLogger.info("${counters[0]} ranges had too few reads; ${counters[1]} ranges had too many reads; ${counters[2]} ranges too few samples with haplotypes")

        return paths.maxBy { it.totalProbability }
    }

    private fun diploidViterbi(
        chrom: String,
        readMap: Map<ReferenceRange, Map<List<String>, Int>>,
        parentList: Set<SampleGamete> = sampleGametesInGraph
    ): PathNode? {
        myLogger.info("Finding path for chromosome $chrom using diploidViterbi")
        val sampleGametePairs =
            parentList.map { firstName -> parentList.map { secondName -> Pair(firstName, secondName) } }
                .flatten()

        //diagnostic counters for discarded ranges:
        //countTooFewReads, countTooManyReadsPerKB, countReadsEqual, countDiscardedRanges
        val counters = IntArray(4) { 0 }

        //create emission probability
        val emissionProb = DiploidEmissionProbability(readMap, graph, probCorrect)

        val myRanges = graph.rangesByContig()[chrom]
        check(!myRanges.isNullOrEmpty()) { "No ranges in contig $chrom" }
        val rangeIterator = myRanges!!.iterator()
        var initialRange = rangeIterator.next()

        while (!useRange(
                readMap[initialRange],
                counters,
                initialRange,
            ) && rangeIterator.hasNext()
        ) {
            initialRange = rangeIterator.next()
        }

        while (!useRange(readMap[initialRange], counters, initialRange)) {
            if (rangeIterator.hasNext()) {
                initialRange = rangeIterator.next()
            } else {
                myLogger.info(
                    "for chromosome $chrom, ${counters[3]} ranges out of ${myRanges.size} were discarded: " +
                            "too few reads - ${counters[0]}, too many reads per kb - ${counters[1]}, all counts equal - ${counters[2]}"
                )
                return null
            }
        }

        var currentNodes = parentList.flatMap { gamete1 ->
            parentList.map { gamete2 ->
                val emissionP = emissionProb.lnProbObsGivenState(Pair(gamete1, gamete2), initialRange)
                PathNode(initialRange, null, listOf(gamete1, gamete2), emissionP)
            }
        }

        while (rangeIterator.hasNext()) {
            //for each node in the next range find the maximum value of path probability * transition
            //actually only two to consider (1-r) * same taxon (switchProbability) and r/(n-1) * all other taxa
            //where n is number of taxa. Since the transition likelihood is the same for all n of the recombinant taxa,
            //only the most likely of those needs to be considered.
            currentNodes = checkNodeProbabilities(currentNodes)

            val nextRange = rangeIterator.next()
            if (!useRange(readMap[nextRange], counters, nextRange)) {
                continue
            }
            val nextNodes = sampleGametePairs
                .map {
                    bestPathToGametePair(
                        it,
                        currentNodes,
                        emissionProb.lnProbObsGivenState(it, nextRange),
                        nextRange
                    )
                }

            currentNodes = nextNodes
        }

        //diagnostic counters for discarded ranges:
        //countTooFewReads (0), countTooManyReadsPerKB (1), countReadsEqual (2), countDiscardedRanges (3)
        if (counters[3] > 0) myLogger.info(
            "for chromosome $chrom, ${counters[3]} ranges out of ${myRanges.size} were discarded: " +
                    "too few reads - ${counters[0]}, too many reads per kb - ${counters[1]}, all counts equal - ${counters[2]}"
        )

        val bestFinalNode = currentNodes.maxBy { it.totalProbability }

        return bestFinalNode
    }

    private fun bestPathToGametePair(gametePair: Pair<SampleGamete, SampleGamete>,
                                     paths: List<PathNode>,
                                     emissionPr: Double,
                                     refRange: ReferenceRange): PathNode {
        val f = inbreedCoef
        val pNoSwitch = sameGameteProbability
        val pSwitch = 1 - pNoSwitch
        val isHomozygous = gametePair.first == gametePair.second
        var bestPreviousPath = paths[0]
        var greatestLnTotalProbability = -Double.MAX_VALUE
        var bestTransitionProbability = 0.0

        for (node in paths) {
            //Let T = the observed transition
            //Pr(T) = Pr(homozygous) * Pr(T | homozygous) + Pr(heterozygous) * Pr(T | heterozygous)
            //Pr(T|homozygous) = if (previous range is homozygous and current range is homozygous) Pr(first path T) else 0
            //Pr(T|heterozygous) = Pr(first path T) * Pr(second path T)
            //P(homozygous) = f; P(heterzygous) = 1 - f
            val previousGametes = node.sampleGametes
            val firstEqual = gametePair.first == previousGametes[0]
            val secondEqual = gametePair.second == previousGametes[1]
            val isPreviousHomozygous = previousGametes[0] == previousGametes[1]

            val firstPathPr = if (firstEqual) pNoSwitch else pSwitch
            val secondPathPr = if (secondEqual) pNoSwitch else pSwitch
            val transitionProbabilityHomozygous = if (isHomozygous && isPreviousHomozygous) firstPathPr else 0.0
            val transitionProbabilityHeterozygous = firstPathPr * secondPathPr
            val transitionProbability = ln(f * transitionProbabilityHomozygous + (1 - f) * transitionProbabilityHeterozygous)
            val candidateTotalProbability = node.totalProbability + transitionProbability
            if (candidateTotalProbability > greatestLnTotalProbability) {
                bestPreviousPath = node
                greatestLnTotalProbability = candidateTotalProbability
                bestTransitionProbability = transitionProbability
            }
        }

        val newSample = listOf(gametePair.first, gametePair.second)
        return PathNode(refRange, bestPreviousPath, newSample, bestPreviousPath.totalProbability + bestTransitionProbability + emissionPr)
    }

    /**
     * If the any of the total probability values in the nodes becomes too small, adjust all the total
     * probabilities by adding (-maximum total probability-0.1).
     */
    private fun checkNodeProbabilities(nodes: List<PathNode>): List<PathNode> {
        val minAllowablePr = -1E200
        val minPr = nodes.minBy { it.totalProbability }.totalProbability
        return if (minPr < minAllowablePr) {
            val maxPr = nodes.maxBy { it.totalProbability }.totalProbability
            val adj = -0.1 - maxPr
            return nodes.map { node -> PathNode(node.refRange, node.parent, node.sampleGametes, node.totalProbability + adj) }
        } else nodes
    }

    /**
     * Determines whether a range meets these conditions:
     * The number of reads >= minReadsPerRange.
     * The number of reads/kB <= maxReadsPerKB.
     * If (minReadsPerRange > 0 and removeEqual is true) not all haplotypes have the same number of reads.
     */
    private fun useRange(
        hapidListCounts: Map<List<String>, Int>?,
        counters: IntArray,
        refRange: ReferenceRange
    ): Boolean {
        //diagnostic counters for discarded ranges:
        //countTooFewReads(0), countTooManyReadsPerKB(1), countOfTooFewGametes(2), countDiscardedRanges(3)
        //
        //data class HapIdSetCount(val hapIdSet : Set<Int>, val count : Int)
        if (minReadsPerRange < 1) return true

        val rangeLength = refRange.end - refRange.start + 1
        val numberOfReads = if (hapidListCounts.isNullOrEmpty()) 0 else {
            hapidListCounts.entries.sumOf { (_, count) -> count }
        }

        if (numberOfReads < minReadsPerRange) {
            counters[0]++
            counters[3]++
            return false
        }
        if (numberOfReads * 1000 / rangeLength > maxReadsPerKB) {
            counters[1]++
            counters[3]++
            return false
        }

        val numberOfSampleGametes = graph.hapIdToSampleGametes(refRange).entries.sumOf { (_, gametes) -> gametes.size }
        if (minGametesPerRange > numberOfSampleGametes) {
            counters[2]++
            counters[3]++
            return false
        }

        return true
    }
}

