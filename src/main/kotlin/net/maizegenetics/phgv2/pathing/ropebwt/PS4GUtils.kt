package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedWriter
import net.maizegenetics.phgv2.api.SampleGamete

data class PS4GData(val gameteList: List<Int>, val pos: Int, val count: Int)

class PS4GUtils {
    companion object {
        fun buildOutputFileName(inputFile: String, outputDir: String) : String {
            val fileName = inputFile.split("/").last().removeSuffix(".txt")
            return "$outputDir/${fileName}_ps4g.txt"
        }
        fun writeOutPS4GFile(pS4GData: List<PS4GData>, sampleGameteCount: Map<SampleGamete,Int>, gameteToIdxMap: Map<SampleGamete,Int>, outputFile: String, header: List<String>, cliCommand: String) {
            bufferedWriter(outputFile).use { writer ->
                writer.write("#PS4G\n")
                header.forEach { writer.write("#$it\n") }
                writer.write("#Command: $cliCommand\n")
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
    }
}