package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.pathing.PathFinderWithViterbiHMM.PathNode
import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import net.maizegenetics.phgv2.utils.Position
import org.apache.logging.log4j.LogManager
import kotlin.collections.forEach
import kotlin.math.E
import kotlin.math.log

class PathFinderHMMPS4G(val probCorrect: Double,
                        val sameGameteProbability: Double,
                        val inbreedCoef: Double) {

    private val myLogger = LogManager.getLogger(PathFinderHMMPS4G::class.java)
    private data class HaploidPathNode(val nodePosition: Position, val parent: HaploidPathNode?, val gameteIndex: Int, val pathProbability: Double)
    private data class DiploidPathNode(val nodePosition: Position, val parent: DiploidPathNode?, val gameteIndex1: Int, val gameteIndex2: Int, val pathProbability: Double)

    fun findHaploidPath(contig: String, gameteIndexMap: Map<Int,String>, readMap: Map<Int, MutableList<Ps4gGameteSet>>) : List<Pair<Position, String>> {
        myLogger.info("Finding haploid path for contig $contig")

        //get the gamete index set
        val gameteIndexSet = gameteIndexMap.keys
        val terminalPathNode = haploidViterbi(contig, readMap, gameteIndexSet)
        var activePathNode = terminalPathNode

        val resultList = mutableListOf<Pair<Position,String>>()
        while (activePathNode.nodePosition.position > 0) {
            resultList.add(Pair(activePathNode.nodePosition, gameteIndexMap[activePathNode.gameteIndex]?: "none"))
            activePathNode = activePathNode.parent ?: HaploidPathNode(Position("na", 0), null, 0, 0.0)
        }

        //at this point the result list is the path in reverse order, so it must be reversed.
        resultList.reverse()
        return resultList
    }

    fun findDiploidPath(contig: String, gameteIndexMap: Map<Int,String>, readMap: Map<Int, MutableList<Ps4gGameteSet>>): List<Triple<Position, String, String>> {
        myLogger.info("Finding diploid path for contig $contig")

        //get the gamete index set
        val gameteIndexSet = gameteIndexMap.keys
        val terminalPathNode = diploidViterbi(contig, readMap, gameteIndexSet)
        var activePathNode = terminalPathNode

        val resultList = mutableListOf<Triple<Position, String, String>>()
        while (activePathNode.nodePosition.position > 0) {
            val genome1 = gameteIndexMap[activePathNode.gameteIndex1] ?: "none"
            val genome2 = gameteIndexMap[activePathNode.gameteIndex2] ?: "none"
            resultList.add(Triple(activePathNode.nodePosition, genome1, genome2))
            activePathNode = activePathNode.parent ?: DiploidPathNode(Position("na", 0), null, 0, 0, 0.0)
        }

        //at this point the result list is the path in reverse order, so it must be reversed.
        resultList.reverse()
        return resultList
    }

    private fun haploidViterbi(
        chrom: String,
        readMap: Map<Int, MutableList<Ps4gGameteSet>>,
        gameteIndexSet: Set<Int>
    ): HaploidPathNode {
        myLogger.info("Finding path for chromosome $chrom using haploidViterbi")
        val switchProbability = 1 - sameGameteProbability
        val numberOfSampleGametes = gameteIndexSet.size
        val logSwitch = log(switchProbability / (numberOfSampleGametes.toDouble() - 1.0), E)
        val logNoSwitch = log(sameGameteProbability, E)

        //create emission probability
        val emissionProb = HaploidPS4GEmissionProbability(readMap, gameteIndexSet, probCorrect)

        var paths = ArrayList<HaploidPathNode>()

        val sortedPositionList = readMap.keys.toList().sorted()

        //make the initial PathNode for each node in the first range
        //the path nodes keep track of all the nodes on a path by holding a link to the parent node of each node
        for (index in gameteIndexSet) {
            paths.add(HaploidPathNode(Position("NA", -1), null, index, 0.0))
        }

        //for each reference range update paths with probabilities and new nodes
        for (positionIndex in sortedPositionList) {
            val currentPosition = Position(chrom, positionIndex)
            //for each node at this position find the maximum value of path probability * transition
            //actually only two to consider (1-r) * same taxon (switchProbability) and r/(n-1) * all other taxa
            //where n is number of taxa. Since the transition likelihood is the same for all n of the recombinant taxa,
            //only the most likely of those needs to be considered.
//            val nextRange = rangeIterator.next()

//            sampleGameteToHaplotypeMap = graph.hapIdToSampleGametes(nextRange).map { (hap, gameteList) ->
//                gameteList.map { Pair(it, hap) }
//            }.flatten().toMap()

//            if (!useRange(readMap[nextRange], counters, nextRange)) {
//                continue
//            }

            val newPaths = ArrayList<HaploidPathNode>()

            //TODO revisit whether randomly picking best is a good strategy
            //choose the most probable path from the previous range. If more than one, any one will do.
            //the most likely path in the previous range
            val maxProb = paths.maxOfOrNull { it.pathProbability } ?: 0.0
            val bestGameteIndices = paths.filter { it.pathProbability == maxProb }.map { it.gameteIndex }.toSet()

            //map of gameteIndex to path node for the previous range
            val gameteToPath = paths.associateBy { it.gameteIndex }

            //iterate over gametes for the new range
            //if switching to a new parent, the best path will be from the highest probability path at the previous position
            val probSwitch = maxProb + logSwitch
            gameteIndexSet.forEach { gameteIndex ->
                if (bestGameteIndices.contains(gameteIndex) ) {
                    //this is the best path for this node since is also the no switch path (same gamete)
                    val parentPath = gameteToPath[gameteIndex]
                    newPaths.add(
                        HaploidPathNode(
                            currentPosition,
                            parentPath,
                            gameteIndex,
                            parentPath!!.pathProbability + logNoSwitch + emissionProb.getLnProbObsGivenState(gameteIndex, currentPosition.position)
                        )
                    )
                } else {
                    //since the best path is a switch path (switching gametes), must compare it to the no switch path
                    val samePath = gameteToPath[gameteIndex]
                    val probNoSwitch = samePath!!.pathProbability + logNoSwitch
                    if (probSwitch > probNoSwitch) {
                        //switch to the any of the gametes that has max Prob
                        val parentPath = gameteToPath[bestGameteIndices.random()]
                        newPaths.add(
                            HaploidPathNode(
                                currentPosition,
                                parentPath,
                                gameteIndex,
                                probSwitch + emissionProb.getLnProbObsGivenState(gameteIndex, currentPosition.position)
                            )
                        )
                    }
                    else {
                        val parentPath = gameteToPath[gameteIndex]
                            newPaths.add(
                            HaploidPathNode(
                                currentPosition,
                                samePath,
                                gameteIndex,
                                probNoSwitch + emissionProb.getLnProbObsGivenState(
                                    gameteIndex,
                                    currentPosition.position
                                )
                            )
                        )
                    }
                }
            }

            paths = newPaths
        }

        return paths.maxBy { it.pathProbability }
    }

    private fun diploidViterbi(
        chrom: String,
        readMap: Map<Int, MutableList<Ps4gGameteSet>>,
        gameteIndexSet: Set<Int>
    ): DiploidPathNode {


    }
}