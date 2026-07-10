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