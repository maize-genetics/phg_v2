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

        private data class PathNode(
            val parent: PathNode?,
            val sample: SampleGamete,
            val haplotype: String?,
            val totalProbability: Double
        )

        init {
            findLikelyParents = useLikelyParents && (maxParents < graph.sampleGametesInGraph().size || minCoverage < 1.0)
            parentFinder = if (findLikelyParents) MostLikelyParents(graph) else null
        }

        /**
         * For a map (list of hapids -> count of reads hitting the set), find the most likely path through the graph.
         * If useLikelyParents = true, then a likely parents list is generated for each sample and only those parents are
         * candiates for the path nodes for that sample.
         *
         * @param readMap a map of [ReferenceRange] to map(list of hapids to count of reads mapping to that list)
         */
        fun findBestHaploidPath(readMap: Map<ReferenceRange, Map<List<String>, Int>>): Pair<List<String>, List<MostLikelyParents.ParentStats>> {
            val haplotypeList = mutableListOf<String>()
            val likelyParentList = if (findLikelyParents) parentFinder!!
                .findMostLikelyParents(readMap, maxParents = maxParents, minCoverage = minCoverage)
            else listOf()

            graph.contigs.forEach { chr ->
                if (findLikelyParents) {
                    haplotypeList.addAll(haploidViterbi(chr, readMap, likelyParentList.map { it.parent }.toSet()))
                } else {
                    haplotypeList.addAll(haploidViterbi(chr, readMap))
                }

            }
            return Pair(haplotypeList, likelyParentList)
        }

        fun findBestDiploidPath(readMap: Map<ReferenceRange, Map<List<String>, Int>>): List<List<String>> {

            val haplotypeList = listOf(mutableListOf(), mutableListOf<String>())
            graph.contigs.forEach { chr ->
                val start = System.nanoTime()
                val chrLists = diploidViterbi(chr, readMap)
                myLogger.info("Elapsed time for chr $chr: ${(System.nanoTime() - start) / 1e9} sec.")
                haplotypeList[0].addAll(chrLists[0])
                haplotypeList[1].addAll(chrLists[1])
            }
            return haplotypeList
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
        ): List<String> {
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

            var rangeIndex = 0

            val rangeIterator = chrRanges.iterator()
            var initialRange = rangeIterator.next()

            while (!useRange(
                    readMap[initialRange],
                    counters,
                    initialRange
                ) && rangeIterator.hasNext()
            ) {
                initialRange = rangeIterator.next()
                rangeIndex++
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
                paths.add(PathNode(null, gamete, myHaplotype, emissionProb.getLnProbObsGivenState(myHaplotype, initialRange)))
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

                rangeIndex++
                if (!useRange(readMap[nextRange], counters, nextRange)) {
                    continue
                }
                val newPaths = ArrayList<PathNode>()

                //TODO revisit whether randomly picking best is a good strategy
                //TODO figure out how to test lines 217-222
                //choose the most probable path from the previous range. If more than one, any one will do.
                //the most likely path in the previous range
                val bestPath = paths.maxByOrNull { it.totalProbability }
                check(bestPath != null) { "no most likely path before range at ${nextRange.contig}:${nextRange.start}" }

                //map of sampleGamete to path node for the previous range
                val gameteToPath = paths.associateBy { it.sample }

                //iterate over gametes for the new range
                val probSwitch = bestPath.totalProbability + logSwitch
                gameteSet.forEach { sampleGamete ->
                    val gameteHaplotype = sampleGameteToHaplotypeMap[sampleGamete]
                    if (bestPath.sample == sampleGamete) {
                        //this is the best path for this node since is also the no switch path (same gamete)
                        newPaths.add(
                            PathNode(
                                bestPath,
                                sampleGamete,
                                gameteHaplotype,
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
                                bestPath,
                                sampleGamete,
                                gameteHaplotype,
                                probSwitch + emissionProb.getLnProbObsGivenState(gameteHaplotype, nextRange)
                            )
                        )
                        else newPaths.add(
                            PathNode(
                                samePath,
                                sampleGamete,
                                gameteHaplotype,
                                probNoSwitch + emissionProb.getLnProbObsGivenState(gameteHaplotype, nextRange)
                            )
                        )


                    }

                    paths = newPaths

                }
            }

            //countTooFewReads, countTooManyReadsPerKB, countOfTooFewGametes, countDiscardedRanges
            myLogger.info("Finished processing reads for a sample: ${counters[3]} ranges discarded out of ${chrRanges.size}.")
            myLogger.info("${counters[0]} ranges had too few reads; ${counters[1]} ranges had too many reads; ${counters[2]} ranges too few samples with haplotypes")

            //terminate
            //back track the most likely path to get a HaplotypeNode List
                var currentPathnode = paths.maxByOrNull { it.totalProbability }
                val haplotypeList = mutableListOf<String>()
                while (currentPathnode != null) {
                    val node = currentPathnode.parent
                    if (node?.haplotype != null) haplotypeList.add(node.haplotype)
                    currentPathnode = node
                }

                //the resulting node list is last to first, so reverse it
                haplotypeList.reverse()
                return haplotypeList
        }

    data class DiploidPathNode(val gamete1: SampleGamete, val gamete2: SampleGamete, val refRange: ReferenceRange) {
        fun gametePair(): Pair<SampleGamete, SampleGamete> = Pair(gamete1, gamete2)
    }
    data class DiploidPath(val nodeList: MutableList<DiploidPathNode>, val totalProbability: Double) {
        fun terminatesIn(gametePair: Pair<SampleGamete, SampleGamete>) = nodeList.last().gametePair() == gametePair

        fun terminalGametePair(): Pair<SampleGamete, SampleGamete> = nodeList.last().gametePair()
    }


    /**
     * Creates a new path by appending gametePair to the previous path ana adding probability to totalProbability
     */
    private fun newDiploidPath(
        refRange: ReferenceRange,
        path: DiploidPath,
        gametePair: Pair<SampleGamete, SampleGamete>,
        probability: Double
    ): DiploidPath {
        val nodeList = path.nodeList.toMutableList()
        nodeList.add(DiploidPathNode(gametePair.first, gametePair.second, refRange))
        return DiploidPath(nodeList, path.totalProbability + probability)
    }

        private fun diploidViterbi(
            chrom: String,
            readMap: Map<ReferenceRange, Map<List<String>, Int>>
        ): List<List<String>> {
            myLogger.info("Finding path for chromosome $chrom using diploidViterbi")
            val sampleGameteSet = graph.sampleGametesInGraph()
            val sampleGametePairs =
                sampleGameteSet.map { firstName -> sampleGameteSet.map { secondName -> Pair(firstName, secondName) } }.flatten()

            //diagnostic counters for discarded ranges:
            //countTooFewReads, countTooManyReadsPerKB, countReadsEqual, countDiscardedRanges
            val counters = IntArray(4) { 0 }

            //create emission probability
            val emissionProb = DiploidEmissionProbability(readMap, graph, probCorrect)

            //natural logs needed for transition probabilities
            val lnProbabilitySameGamete = ln(sameGameteProbability)
            val lnProbabilityDiffGamete = ln(1 - sameGameteProbability)

            val myRanges = graph.rangesByContig()[chrom]
            check(myRanges.isNullOrEmpty()) {"No ranges in contig $chrom"}
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

            val paths = mutableListOf<DiploidPath>()
            for (gamete1 in sampleGameteSet) {
                for (gamete2 in sampleGameteSet) {
                    val emissionP = emissionProb.lnProbObsGivenState(Pair(gamete1, gamete2), initialRange)
                    val thisNode = DiploidPathNode(gamete1, gamete2, initialRange)
                    paths.add(DiploidPath(mutableListOf(thisNode), emissionP))
                }
            }

            while (rangeIterator.hasNext()) {
                //for each node in the next range find the maximum value of path probability * transition
                //actually only two to consider (1-r) * same taxon (switchProbability) and r/(n-1) * all other taxa
                //where n is number of taxa. Since the transition likelihood is the same for all n of the recombinant taxa,
                //only the most likely of those needs to be considered.
                val nextRange = rangeIterator.next()

                if (!useRange(readMap[nextRange], counters, nextRange)) {
                    continue
                }

                //choose the most probable path from the previous range. If more than one, any one will do.
                val bestPath = paths.maxByOrNull { it.totalProbability }  //the most likely path in the previous range
                check(bestPath != null) { "no most likely path before range at ${nextRange.contig}:${nextRange.start}" }

                val nextPaths = mutableListOf<DiploidPath>()

                for (gametePair in sampleGametePairs) {
                    val pairEmissionProbability = emissionProb.lnProbObsGivenState(gametePair, nextRange)

                    //find the most likely path leading to this node
                    //if the bestPath ends in the same taxa pair as this node, then that is the most likely path.
                    if (bestPath.terminatesIn(gametePair)) {
                        nextPaths.add(
                            newDiploidPath(
                                nextRange,
                                bestPath,
                                gametePair,
                                pairEmissionProbability + lnProbabilitySameGamete + lnProbabilitySameGamete
                            )
                        )
                    } else {
                        //Need to find the path in paths for which total probability + transition probability is greatest
                        //Step 1. Create three lists of paths corresponding to the transition types of same-same, same-diff, and diff-diff
                        //Step 2. Find the path with the highest totalProbability within each transition type
                        //Step 3. For those three, pick the highest new total probability to create a new path
                        //step 1
                        val sameSameList = mutableListOf<DiploidPath>()
                        val sameDiffList = mutableListOf<DiploidPath>()
                        val diffDiffList = mutableListOf<DiploidPath>()
                        for (path in paths) {
                            val lastGametes = path.terminalGametePair()
                            val firstEqual = gametePair.first == lastGametes.first
                            val secondEqual = gametePair.second == lastGametes.second
                            if (firstEqual) {
                                if (secondEqual) sameSameList.add(path)
                                else sameDiffList.add(path)
                            } else {  //!firstEqual
                                if (secondEqual) sameDiffList.add(path)
                                else diffDiffList.add(path)
                            }
                        }

                        //step2
                        check(sameSameList.size == 1) { "Only one path with same-same transition expected in range $nextRange. There were ${sameSameList.size}" }
                        val sameSameBestPath = sameSameList[0]
                        val sameSameProbability =
                            sameSameBestPath.totalProbability + lnProbabilitySameGamete + lnProbabilitySameGamete
                        val sameDiffBestPath = sameDiffList.maxBy { it.totalProbability }
                        val sameDiffProbability =
                            sameDiffBestPath.totalProbability + lnProbabilitySameGamete + lnProbabilityDiffGamete
                        val diffDiffBestPath = diffDiffList.maxBy { it.totalProbability }
                        val diffDiffProbability =
                            diffDiffBestPath.totalProbability + lnProbabilityDiffGamete + lnProbabilityDiffGamete

                        //step 3
                        if (sameSameProbability >= sameDiffProbability && sameSameProbability >= diffDiffProbability)
                            nextPaths.add(
                                newDiploidPath(
                                    nextRange,
                                    sameSameBestPath,
                                    gametePair,
                                    pairEmissionProbability + sameSameProbability
                                )
                            )
                        else if (sameDiffProbability >= diffDiffProbability)
                            nextPaths.add(
                                newDiploidPath(
                                    nextRange,
                                    sameDiffBestPath,
                                    gametePair,
                                    pairEmissionProbability + sameDiffProbability
                                )
                            )
                        else nextPaths.add(
                            newDiploidPath(
                                nextRange,
                                sameDiffBestPath,
                                gametePair,
                                pairEmissionProbability + sameDiffProbability
                            )
                        )

                    }
                }
            }

            val finalBestPath = paths.maxBy { it.totalProbability }
            //convert SampleGametes to haplotypes
            val hapids1 = mutableListOf<String>()
            val hapids2 = mutableListOf<String>()
            finalBestPath.nodeList.forEach { node ->
                val hap1 = graph.sampleToHapId(node.refRange, node.gamete1)
                if (hap1 != null) hapids1.add(hap1)
                val hap2 = graph.sampleToHapId(node.refRange, node.gamete2)
                if (hap2 != null) hapids2.add(hap2)
            }
            return listOf(hapids1, hapids2)
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

//Comment out Diploid code, which will be migrated later
//    class DiploidTransitionProbabilityWithInbreeding(
//        val numberOfTaxaNames: Int,
//        haploidNoSwitchProb: Double,
//        val f: Double
//    ) {
//        private val pNoSwitch: Double
//        private val pSingleSwitch: Double
//        private val pDoubleSwitch: Double
//        private val pAAtoAB: Double
//        private val pAAtoBB: Double
//        private val pAAtoBC: Double
//        private val numberOfPairs = numberOfTaxaNames * numberOfTaxaNames
//        private val probabilityMatrix = Array(numberOfPairs) { DoubleArray(numberOfPairs) { 0.0 } }
//
//        init {
//            val haploidSwitchProb = (1.0 - haploidNoSwitchProb) / (numberOfTaxaNames.toDouble() - 1.0)
//            pNoSwitch = haploidNoSwitchProb * haploidNoSwitchProb
//            pSingleSwitch = haploidNoSwitchProb * haploidSwitchProb
//            pDoubleSwitch = haploidSwitchProb * haploidSwitchProb
//            pAAtoAB = (1.0 - f) * pSingleSwitch
//            pAAtoBB = f * pSingleSwitch + (1.0 - f) * pDoubleSwitch
//            pAAtoBC = (1.0 - f) * pDoubleSwitch
//
//            for (from in 0 until numberOfPairs) {
//                for (to in 0 until numberOfPairs) {
//                    probabilityMatrix[to][from] = log(transitionProbability(pairFromIndex(from), pairFromIndex(to)), E)
//                }
//            }
//        }
//
//        fun pairFromIndex(index: Int): Pair<Int, Int> {
//            return Pair(index / numberOfTaxaNames, index % numberOfTaxaNames)
//        }
//
//        /**
//         * Adds the vector of [fromProbability] to a vector of transition probabilities and returns
//         * the index of the maximum value and the maximum value as a Pair<Int,Double>
//         */
//        fun maxIndexAndProbabilityForTarget(fromProbability: DoubleArray, to: Int): Pair<Int, Double> {
//            val transProbs = probabilityMatrix[to]
//
//            var maxIndex = 0
//            var maxProb = fromProbability[maxIndex] + transProbs[maxIndex]
//            for (ndx in 1 until fromProbability.size) {
//                val prob = fromProbability[ndx] + transProbs[ndx]
//                if (prob > maxProb) {
//                    maxProb = prob
//                    maxIndex = ndx
//                }
//            }
//            return Pair(maxIndex, maxProb)
//        }
//
//        fun transitionProbability(from: Pair<Int, Int>, to: Pair<Int, Int>): Double {
////        diploid transition as a function of f (inbreeding coefficient)
////        The coefficient of inbreeding of an individual is the probability that two alleles at any locus in an individual are identical by descent from the common ancestor(s) of the two parents
////
////        p(A,A -> A,A) = P(no switch) * P(no switch) [same whether ibd or not]
////        p(A,A -> A,B | not ibd) = P(single switch) * P(no switch)
////        p(A,A -> A,B | ibd) = 0
////        p(A,A -> B,B  | not ibd) = P(single switch) * P(single switch)
////        p(A,A -> B,B  | ibd) = P(single switch)
////        p(A,A -> B,C | not ibd) = P(single switch) * P(single switch)
////        p(A,A -> B,C | ibd) = 0
////        transition from het is always not ibd
//
//
//            return if (from.first == from.second) {
//                when {
//                    from.first == to.first && from.second == to.second -> pNoSwitch
//                    from.first == to.first || from.second == to.second -> pAAtoAB
//                    to.first == to.second -> pAAtoBB
//                    else -> pAAtoBC
//                }
//            } else {
//                when {
//                    from.first == to.first && from.second == to.second -> pNoSwitch
//                    from.first == to.first || from.second == to.second -> pSingleSwitch
//                    else -> pDoubleSwitch
//                }
//            }
//        }
//
//        fun lnTransitionProbability(from: Pair<Int, Int>, to: Pair<Int, Int>): Double {
//            return probabilityMatrix[to.first * numberOfTaxaNames + to.second][from.first * numberOfTaxaNames + from.second]
//        }
//
//    }
