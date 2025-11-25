package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedWriter
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.utils.Position
import kotlin.math.max

data class PS4GData(val gameteList: List<Int>, val refPos: Position, val count: Int,
                    val numMapped: Int,
                    val numMappedOnMainContig: Int,
                    val totalMaxPosDistOnMainContig: Int)


data class GameteIdsWithMappingStats(
    val gameteIdsHit: List<Int>,
    val numMapped: Int,
    val numMappedOnMainContig: Int,
    val maxPosDistOnMainContig: Int
)

data class PS4GCountValue(
    val count: Int,
    val numMapped: Int,
    val numMappedOnMainContig: Int,
    val totalMaxPosDistOnMainContig: Int
)

class PS4GUtils {
    companion object {
        fun buildOutputFileName(inputFile: String, outputDir: String, sampleGamete: String = "") : String {
            val fileName = inputFile.split("/").last().removeSuffix(".txt").removeSuffix(".bed")
            return if (sampleGamete.isNotEmpty()) {
                "$outputDir/${fileName}_${sampleGamete}_ps4g.txt"
            } else {
                "$outputDir/${fileName}_ps4g.txt"
            }
        }
        fun writeOutPS4GFile(pS4GData: List<PS4GData>, sampleGameteCount: Map<SampleGamete,Int>, gameteToIdxMap: Map<SampleGamete,Int>, outputFile: String, header: List<String>, cliCommand: String) {
            bufferedWriter(outputFile).use { writer ->
                writer.write("#PS4G\n")
                writer.write("#version=3.0\n")
                header.forEach { writer.write("#$it\n") }
                writer.write("#Command: $cliCommand\n")
                writer.write("#TotalUniqueCounts: ${pS4GData.map { it.count }.sum()}\n")
                writer.write("#gamete\tgameteIndex\tcount\n")
                gameteToIdxMap.toList().sortedBy { it.second }.forEach { (sampleGamete, index) ->
                    val count = sampleGameteCount.getOrDefault(sampleGamete, 0)
                    writer.write("#$sampleGamete\t$index\t$count\n")
                }
                writer.write("gameteSet\trefContig\trefPosBinned\tcount\tnumMappings\tpropOnTopContig\tavgPosVariation\n")
                pS4GData.forEach { (gameteList, pos, count,numMapped, numMappedOnMainContig, totalMaxPosDistOnMainContig) ->
                    writer.write("${gameteList.joinToString(",")}\t${pos.contig}\t${pos.position}\t$count\t" +
                            "${numMapped}\t" +
                            "${numMappedOnMainContig.toDouble()/max(numMapped,1)}\t" + //Need to have a max to avoid divide by zero
                            "${totalMaxPosDistOnMainContig.toDouble()/max(numMappedOnMainContig,1)}\n")
                }

            }
        }
//        fun convertCountMapToPS4GData(countMap: Map<Pair<Position,List<Int>>, Int>, sortPositions: Boolean = true) : List<PS4GData> {
//            return if(sortPositions) {
//                countMap.map { (pair, count) ->
//                    PS4GData(pair.second.sorted(), pair.first, count)
//                }.sortedBy { it.refPos } //Need to decode it because chromosome might be in an unexpected order
//            }
//            else {
//                countMap.map { (pair, count) ->
//                    PS4GData(pair.second, pair.first, count)
//                }
//            }
//        }

        fun convertCountMapToPS4GData(countMap: Map<Pair<Position,List<Int>>, PS4GCountValue>, sortPositions: Boolean = true) : List<PS4GData> {
            return if(sortPositions) {
                countMap.map { (pair, countValue) ->
                    PS4GData(pair.second.sorted(), pair.first, countValue.count,
                        countValue.numMapped, countValue.numMappedOnMainContig, countValue.totalMaxPosDistOnMainContig)
                }.sortedBy { it.refPos } //Need to decode it because chromosome might be in an unexpected order
            }
            else {
                countMap.map { (pair, countValue) ->
                    PS4GData(pair.second, pair.first, countValue.count,
                        countValue.numMapped, countValue.numMappedOnMainContig, countValue.totalMaxPosDistOnMainContig)
                }
            }
        }

        fun incrementCountValue(countMap: MutableMap<Pair<Position, List<Int>>, PS4GCountValue>, posToGameteSetAndStats: Pair<Position, GameteIdsWithMappingStats>) {
            val key = Pair(posToGameteSetAndStats.first, posToGameteSetAndStats.second.gameteIdsHit)
            val currentPS4GCountValue = countMap.getOrDefault(key, PS4GCountValue(0,0,0,0))
            val newPS4GCountValue = PS4GCountValue(
                currentPS4GCountValue.count + 1,
                currentPS4GCountValue.numMapped + posToGameteSetAndStats.second.numMapped,
                currentPS4GCountValue.numMappedOnMainContig + posToGameteSetAndStats.second.numMappedOnMainContig,
                currentPS4GCountValue.totalMaxPosDistOnMainContig + posToGameteSetAndStats.second.maxPosDistOnMainContig
            )
            countMap[key] = newPS4GCountValue
        }
    }
}