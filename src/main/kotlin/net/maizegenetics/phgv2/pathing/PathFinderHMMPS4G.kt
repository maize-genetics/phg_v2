package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import net.maizegenetics.phgv2.utils.Position
import org.apache.logging.log4j.LogManager
import kotlin.collections.forEach
import kotlin.math.E
import kotlin.math.ln
import kotlin.math.log

class PathFinderHMMPS4G(val probCorrect: Double,
                        val sameGameteProbability: Double,
                        val inbreedCoef: Double) {

    private val myLogger = LogManager.getLogger(PathFinderHMMPS4G::class.java)
    private data class HaploidPathNode(val nodePosition: Position, val parent: HaploidPathNode?, val gameteIndex: Int, val pathProbability: Double) {
        override fun hashCode(): Int {
            return nodePosition.position
        }
    }
    private data class DiploidPathNode(val nodePosition: Position, val parent: DiploidPathNode?, val gameteIndex1: Int, val gameteIndex2: Int, val pathProbability: Double) {
        override fun hashCode(): Int {
            return nodePosition.position
        }
    }

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
        myLogger.info("Finished path for contig $contig")
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
                    } else {
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

    private fun checkPathProbabilitiesForOverflow(paths: ArrayList<HaploidPathNode>) : ArrayList<HaploidPathNode> {
        val minLimit = Double.MAX_VALUE * - 0.5


        val minProb = paths.minOf { it.pathProbability }
        return if (minProb < minLimit) {
            val maxProb = paths.maxOf { it.pathProbability }
            val adjustedPaths = ArrayList<HaploidPathNode>()
            paths.forEach { path -> adjustedPaths.add(PathFinderHMMPS4G.HaploidPathNode(path.nodePosition, path.parent, path.gameteIndex, path.pathProbability - maxProb)) }
            return adjustedPaths
        } else paths
    }

    private fun diploidViterbi(
        chrom: String,
        readMap: Map<Int, MutableList<Ps4gGameteSet>>,
        gameteIndexSet: Set<Int>
    ): DiploidPathNode {
        myLogger.info("Finding path for chromosome $chrom using diploidViterbi")

        val ngenomes = gameteIndexSet.size
        val transition = DiploidTransitionWithInbreeding(sameGameteProbability, inbreedCoef, ngenomes)

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

        //if the inbreeding coefficient is > 0, then homozygous genotypes are more probable then heterozygous ones.
        // Set the homozygous genotype probabilities = inbreedingCoef / ngenomes
        //Set nHeterozygousGenotypes = ngenomes * (ngenomes - 1)

        val homozygoteProb = ln(1.0 / ngenomes.toDouble())
        val heterozygoteProb = ln(1.0 / (ngenomes * ngenomes - ngenomes).toDouble())
        for (index1 in gameteIndexSet) {
            for(index2 in gameteIndexSet) {
                val initProb = if (index1 == index2) homozygoteProb else heterozygoteProb
                paths.add(DiploidPathNode(Position("NA", -1), null, index1, index2, initProb))
            }
        }

        for (positionIndex in sortedPositionList) {
            val currentPosition = Position(chrom, positionIndex)
            val newPaths = ArrayList<DiploidPathNode>()
            val maxProb = paths.maxOfOrNull { it.pathProbability } ?: 0.0
//            val mostProbableNodes = paths.filter { it.pathProbability == maxProb }
            //choose one node to use for the switch paths.
            //val mostProbableParentNode = mostProbableNodes.random()
            //val probSwitch = maxProb + logSwitch //the maximum total probability for any path with a switch

            //create a mapping of indices to nodes
            val indicesToNodeMap = paths.associate { node -> Pair(Pair(node.gameteIndex1, node.gameteIndex2), node) }

            //iterate through all possible ordered pairs of gamete indices
            for (index1 in gameteIndexSet) {
                for(index2 in gameteIndexSet) {
                    val indexPairEmissionProb = emissionProb.getLnProbObsGivenState(index1, index2, positionIndex)
                    //for this ordered pair of indices find the node at the previous position that results in the highest total probability
                    //then make that the parent node of the new node

                    /*
                    * Calculate the same path probabilty = previous node probability + transition probability
                    * Calculate previous position path probability * transition probability for every node at the previous position
                    * Find the maximum new path probability
                    *
                    * If this equals the same path probability, use same path for the new path.
                    * Otherwise, create a list of path nodes that result in the maximum path probability.
                    * Choose (at random) one of the maximum probability nodes as the previous node.
                    * */

                    val samePathNode = indicesToNodeMap[Pair(index1, index2)]
                    check(samePathNode != null) {"indices $index1, $index2 did not map to a node"}
                    val samePathProb = samePathNode.pathProbability +
                            transition.calculate(Pair(samePathNode.gameteIndex1, samePathNode.gameteIndex2), Pair(index1, index2))

                    val candidatePathList = paths
                        .map { path -> Pair(path, path.pathProbability + transition.calculate(Pair(path.gameteIndex1, path.gameteIndex2), Pair(index1, index2))) }
                    val maxProb = candidatePathList.maxOf { it.second }

                    if (samePathProb >= maxProb) {
                        val newPathProb = samePathProb + indexPairEmissionProb
                        newPaths.add(DiploidPathNode(currentPosition, samePathNode, index1, index2, newPathProb))
                    } else {
                        val bestPath = candidatePathList.filter { it.second == maxProb }.random()
                        val newPathProb = bestPath.second + indexPairEmissionProb
                        newPaths.add(DiploidPathNode(currentPosition, bestPath.first, index1, index2, newPathProb))
                    }

                }
            }
            paths = newPaths
        }
        return paths.maxBy { it.pathProbability }
    }

    class DiploidTransition(val pNoSwitch: Double, numberOfGenomes: Int) {
        val pSwitch = (1.0 - pNoSwitch) / (numberOfGenomes - 1)
        val pss = ln(pSwitch * pSwitch)
        val psn = ln(pSwitch * pNoSwitch)
        val pnn = ln(pNoSwitch * pNoSwitch)
        val samePathProbability = pnn

        fun calculate(from: Pair<Int,Int>, to: Pair<Int,Int>): Double {
            return if (from.first == to.first) {
                if (from.second == to.second) pnn else psn
            } else {
                if (from.second == to.second) psn else pss
            }
        }
    }

    class DiploidTransitionWithInbreeding(val pNoSwitch: Double, val inbreedingCoef: Double, numberOfGenomes: Int) {
        val pSwitch = (1.0 - pNoSwitch) / (numberOfGenomes - 1)
        val pss = pSwitch * pSwitch
        val psn = pSwitch * pNoSwitch
        val pnn = pNoSwitch * pNoSwitch

        fun calculate(from: Pair<Int,Int>, to: Pair<Int,Int>): Double {
            val transitionP =  when (inbreedingCoef) {
                0.0 -> probabilityForF0(from, to)
                1.0 -> probabilityForF1(from, to)
                else -> (1.0 - inbreedingCoef) * probabilityForF0(from, to) + inbreedingCoef * probabilityForF1(from, to)
            }

            return ln(transitionP)
        }

        fun probabilityForF0(from: Pair<Int,Int>, to: Pair<Int,Int>):Double {
            return if (from.first == to.first) {
                if (from.second == to.second) pnn else psn
            } else {
                if (from.second == to.second) psn else pss
            }
        }

        fun probabilityForF1(from: Pair<Int,Int>, to: Pair<Int,Int>):Double {
            return if (from.first == from.second)  {  //from is homozygous
                if (to.first == to.second) {
                    if (from.first == to.first) pNoSwitch else pSwitch
                } else 0.0
            } else { //from is heterozygous
                if (to.first == to.second) {
                    if (from.first == to.first || from.second == to.second) pSwitch
                    else pss
                } else 0.0
            }

        }

    }

}