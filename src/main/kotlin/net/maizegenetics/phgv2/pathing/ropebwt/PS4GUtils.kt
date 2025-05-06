package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedWriter
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.utils.Position

data class PS4GData(val gameteList: List<Int>, val pos: Int, val count: Int)

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
                header.forEach { writer.write("#$it\n") }
                writer.write("#Command: $cliCommand\n")
                writer.write("#TotalUniqueCounts: ${pS4GData.map { it.count }.sum()}\n")
                writer.write("#gamete\tgameteIndex\tcount\n")
                sampleGameteCount.forEach { (sampleGamete, count) ->
                    writer.write("#$sampleGamete\t${gameteToIdxMap[sampleGamete]}\t$count\n")
                }
                writer.write("gameteSet\tpos\tcount\n")
                pS4GData.forEach { (gameteList, pos, count) ->
                    writer.write("${gameteList.joinToString(",")}\t$pos\t$count\n")
                }

            }
        }
        fun encodePosition(pos: Position, contigIndexMap: Map<String, Int>) : Int {
            //Pack into an Int
            val idx = contigIndexMap[pos.contig]?: throw IllegalArgumentException("Contig ${pos.contig} not found in contigIndexMap")

            return encodePositionFromIdxAndPos(idx, pos.position)
        }

        fun encodePositionNoLookup(pos: Position): Int {
            //make sure that the contig is a number
            val idx = pos.contig.toInt()
            return encodePositionFromIdxAndPos(idx, pos.position)
        }

        fun encodePositionFromIdxAndPos(idx: Int, pos: Int) : Int {
            //Pack into an Int
            //pack last 8 bits of idx into first 8 bits of output then pack the position minus 8 bits into the last 24 bits
            val idxBits = idx and 0xFF //If there are more than 256 contigs this will have unexpected issues
            val posBits = pos/256 // div 256 effectively bitshifts by 8

            return (idxBits shl 24) or posBits //we dont care if its negative as we arent comparing them
        }

        //This is a lossy function as we /256 the position during encoding.  So it will be in a bin of 256
        fun decodePosition(encodedPos : Int) : Position {
            val idx = encodedPos shr 24
            val pos = (encodedPos and 0x0FFFFFFF) * 256
            return Position("$idx", pos)
        }

        /**
         * Function to convert the count map to a PS4GData class for easy export
         */
        fun convertCountMapToPS4GData(countMap: Map<Pair<Int,List<Int>>, Int>, sortPositions: Boolean = true) : List<PS4GData> {
            return if(sortPositions) {
                countMap.map { (pair, count) ->
                    PS4GData(pair.second.sorted(), pair.first, count)
                }.sortedBy { PS4GUtils.decodePosition(it.pos) } //Need to decode it because chromosome might be in an unexpected order
            }
            else {
                countMap.map { (pair, count) ->
                    PS4GData(pair.second, pair.first, count)
                }
            }
        }
    }
}