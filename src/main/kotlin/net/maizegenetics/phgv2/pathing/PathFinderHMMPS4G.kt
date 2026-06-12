package net.maizegenetics.phgv2.pathing

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

    /**
     * Find the best haploid path for a contig. Inputs are the contig name, a map of gamete index to name, and
     * a map of int position (bin) to a list of Ps4gGameteSets. A Ps4gGameteSet contains an array of gamete indices and the count
     * of the number of read hits to that set of gametes. The output is a list of pairs: position, genome name.
     */
    fun findHaploidPath(contig: String, gameteIndexMap: Map<Int,String>, readMap: Map<Int, MutableList<Ps4gGameteSet>>) : List<Pair<Position, String>> {
        myLogger.info("Finding haploid path for contig $contig")

        //get the gamete index set
        val gameteIndexSet = gameteIndexMap.keys
        val terminalPathNode = haploidViterbi(contig, readMap, gameteIndexSet)
        var activePathNode = terminalPathNode

        val resultList = mutableListOf<Pair<Position,String>>()
        while (activePathNode.nodePosition.position > -1) {
            resultList.add(Pair(activePathNode.nodePosition, gameteIndexMap[activePathNode.gameteIndex]?: "none"))
            activePathNode = activePathNode.parent ?: HaploidPathNode(Position("na", 0), null, 0, 0.0)
        }

        //at this point the result list is the path in reverse order, so it must be reversed.
        resultList.reverse()
        return resultList
    }

    /**
     * Find the best diploid path for a contig. Inputs are the contig name, a map of gamete index to name, and
     * a map of int position (bin) to a list of Ps4gGameteSets. A Ps4gGameteSet contains an array of gamete indices and the count
     * of the number of read hits to that set of gametes. The output is a list of triples: position, genome1, genome2.
     */
    fun findDiploidPath(contig: String, gameteIndexMap: Map<Int,String>, readMap: Map<Int, MutableList<Ps4gGameteSet>>): List<Triple<Position, String, String>> {
        myLogger.info("Finding diploid path for contig $contig")

        //get the gamete index set
        val gameteIndexSet = gameteIndexMap.keys
        val terminalPathNode = diploidViterbi(contig, readMap, gameteIndexSet)
        var activePathNode = terminalPathNode

        val resultList = mutableListOf<Triple<Position, String, String>>()
        while (activePathNode.nodePosition.position > -1) {
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
            val newPaths = ArrayList<HaploidPathNode>()

            //choose the most probable path from the previous range. If more than one, any one will do.
            //the most likely path in the previous range
            val maxProb = paths.maxOfOrNull { it.pathProbability } ?: 0.0
            val mostProbablePaths = paths.filter { it.pathProbability == maxProb }.toSet()
            val bestGameteIndices = mostProbablePaths.map { it.gameteIndex }
            val selectedSwitchParent = mostProbablePaths.random()
            val probSwitch = maxProb + logSwitch

            //map of gameteIndex to path node for the previous range
            val gameteToPath = paths.associateBy { it.gameteIndex }

            //iterate over gametes for the new range
            //if switching to a new parent, the best path will be from the highest probability path at the previous position
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
                        newPaths.add(
                            HaploidPathNode(
                                currentPosition,
                                selectedSwitchParent,
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
        myLogger.info("Finding path for chromosome $chrom using diploidViterbi")
        val switchProbability = 1 - sameGameteProbability
        val numberOfSampleGametes = gameteIndexSet.size
        val logSwitch = log(switchProbability / (numberOfSampleGametes.toDouble() - 1.0), E)
        val logNoSwitch = log(sameGameteProbability, E)

        //create emission probability
        val emissionProb = DiploidPS4GEmissionProbability(readMap, gameteIndexSet, probCorrect)

        var paths = ArrayList<DiploidPathNode>()
        val sortedPositionList = readMap.keys.toList().sorted()

        //Make the initial PathNode for each node in the first range.
        //The path nodes keep track of all the nodes on a path by holding a link to the parent node of each node
        //The initial node has position -1 since the first (real) bin is 0. The nodes are assigned equal probability.
        //The value of the probability has no meaning by itself, so all the probabilities can be set to 0.0. What matters
        //is the relative probabilities of the different paths. Since the probabilities are ln(probability) (natural log)
        // this means that all initial probabilities are set to 1.
        for (index1 in gameteIndexSet) {
            for(index2 in gameteIndexSet) {
                paths.add(DiploidPathNode(Position("NA", -1), null, index1, index2, 0.0))
            }
        }

        for (positionIndex in sortedPositionList) {
            val currentPosition = Position(chrom, positionIndex)
            val newPaths = ArrayList<DiploidPathNode>()
            val maxProb = paths.maxOfOrNull { it.pathProbability } ?: 0.0
            val mostProbableNodes = paths.filter { it.pathProbability == maxProb }
            //choose one node to use for the switch paths.
            val mostProbableParentNode = mostProbableNodes.random()
            val probSwitch = maxProb + logSwitch //the maximum total probability for any path with a switch

            //create a mapping of indices to nodes
            val indicesToNodeMap = paths.associate { node -> Pair(Pair(node.gameteIndex1, node.gameteIndex2), node) }

            //iterate through all possible ordered pairs of gamete indices
            for (index1 in gameteIndexSet) {
                for(index2 in gameteIndexSet) {
                    val indexPairEmissionProb = emissionProb.getLnProbObsGivenState(index1, index2, positionIndex)
                    //for this ordered pair of indices find the node at the previous position that results in the highest total probability
                    //then make that the parent node of the new node

                    //if this index pair is in the highest probability set, then choose that as the parent node
                    val samePathNode = indicesToNodeMap[Pair(index1, index2)]
                    check(samePathNode != null) {"indices $index1, $index2 did not map to a node"}
                    if (mostProbableNodes.contains(samePathNode)) {
                        val newPathProbability = samePathNode.pathProbability + logNoSwitch + indexPairEmissionProb
                        newPaths.add(DiploidPathNode(currentPosition, samePathNode,
                            index1, index2, newPathProbability))
                    } else {
                        //otherwise, calculate the total probability if the same index pair is the previous parent (use logNoSwitch)
                        //  and the total probability if any of the mostProbableNodes is the parent (equals probSwitch)
                        //  use the parent that gives the highest total probability
                        val samePathProbability = samePathNode.pathProbability + logNoSwitch
                        if (samePathProbability >= probSwitch) {
                            val newPathProbability = samePathProbability + indexPairEmissionProb
                            newPaths.add(DiploidPathNode(currentPosition, samePathNode,
                                index1, index2, newPathProbability))
                        } else {
                            val newPathProbability = probSwitch + indexPairEmissionProb
                            newPaths.add(DiploidPathNode(currentPosition, mostProbableParentNode,
                                index1, index2, newPathProbability))
                        }
                    }

                }
            }
            paths = newPaths
        }
        return paths.maxBy { it.pathProbability }
    }
}