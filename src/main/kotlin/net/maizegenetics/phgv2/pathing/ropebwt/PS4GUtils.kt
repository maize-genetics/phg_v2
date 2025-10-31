package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedWriter
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.utils.Position

data class PS4GData(val gameteList: List<Int>, val refPos: Position, val count: Int)

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
                writer.write("#version=2.0\n")
                header.forEach { writer.write("#$it\n") }
                writer.write("#Command: $cliCommand\n")
                writer.write("#TotalUniqueCounts: ${pS4GData.map { it.count }.sum()}\n")
                writer.write("#gamete\tgameteIndex\tcount\n")
                gameteToIdxMap.toList().sortedBy { it.second }.forEach { (sampleGamete, index) ->
                    val count = sampleGameteCount.getOrDefault(sampleGamete, 0)
                    writer.write("#$sampleGamete\t$index\t$count\n")
                }
                writer.write("gameteSet\trefContig\trefPosBinned\tcount\n")
                pS4GData.forEach { (gameteList, pos, count) ->
                    writer.write("${gameteList.joinToString(",")}\t${pos.contig}\t${pos.position}\t$count\n")
                }

            }
        }
        fun convertCountMapToPS4GData(countMap: Map<Pair<Position,List<Int>>, Int>, sortPositions: Boolean = true) : List<PS4GData> {
            return if(sortPositions) {
                countMap.map { (pair, count) ->
                    PS4GData(pair.second.sorted(), pair.first, count)
                }.sortedBy { it.refPos } //Need to decode it because chromosome might be in an unexpected order
            }
            else {
                countMap.map { (pair, count) ->
                    PS4GData(pair.second, pair.first, count)
                }
            }
        }
    }
}