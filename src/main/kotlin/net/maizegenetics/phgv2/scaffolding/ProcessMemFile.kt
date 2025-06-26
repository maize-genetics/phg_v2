package net.maizegenetics.phgv2.scaffolding

import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.colormath.parse
import net.maizegenetics.phgv2.pathing.ropebwt.MEM
import net.maizegenetics.phgv2.pathing.ropebwt.RopeBWTUtils
import java.io.File

class ProcessMemFile : CliktCommand(help = "Process Mem File counting number of hits") {
    val memFile by option(help = "Path to the MEM file")
        .required()

    val outputDir by option(help = "Output Directory")
        .required()

    override fun run() {
        // Implement the logic to process the MEM file
        // This could involve reading the file, extracting relevant data,
        // and writing results to an output file or directory.

        processMemFile(memFile, outputDir)
    }

    fun processMemFile(memFile: String, outputDir: String) {
        File(outputDir).mkdirs()

        val mems = bufferedReader(memFile).readLines().map { line -> RopeBWTUtils.parseStringIntoMem(line) }

        val memsCollected = collectMemsByReadName(mems)

        println("Total Reads Into Mems: ${mems.size}")
        writeOutAllMemsNumHits(memsCollected, outputDir)


        val keepMemsLowHits = filterMemsByNumHits(memsCollected)
        writeOutNumHits(keepMemsLowHits, "${outputDir}/lowHitMems.txt")

        val keepLongestMemHit = filterMemsByLength(memsCollected)
        writeOutNumHits(keepLongestMemHit, "${outputDir}/longestHitMems.txt")
    }

    fun collectMemsByReadName(mems: List<MEM>): Map<String,List<MEM>> {
        //we need to walk through the mems and rename their read names by how they are grouped in order

        var currentReadName = mems[0].readName
        val currentGroup = mutableListOf<MEM>(mems[0])
        val groupedMems = mutableMapOf<String,List<MEM>>()
        val readNameCount = mutableMapOf<String, Int>()
        for (i in 1 until mems.size) {
            if (mems[i].readName != currentReadName) {
                if(!readNameCount.containsKey(currentGroup.first().readName)) {
                    readNameCount[currentGroup.first().readName] = 1
                }
                else {
                    readNameCount[currentGroup.first().readName] = readNameCount[currentGroup.first().readName]!! + 1
                }
                groupedMems["${currentGroup.first().readName}_${readNameCount[currentGroup.first().readName]!!}"] = currentGroup.toList()

                currentReadName = mems[i].readName
                currentGroup.clear()
            }

            currentGroup.add(mems[i])

        }
        // Add the last group
        if (currentGroup.isNotEmpty()) {
            if(!readNameCount.containsKey(currentGroup.first().readName)) {
                readNameCount[currentGroup.first().readName] = 1
            } else {
                readNameCount[currentGroup.first().readName] = readNameCount[currentGroup.first().readName]!! + 1
            }
            groupedMems["${currentGroup.first().readName}_${readNameCount[currentGroup.first().readName]!!}"] = currentGroup.toList()
        }

        return groupedMems
    }

    fun writeOutAllMemsNumHits(memsCollected: Map<String,List<MEM>>, outputDir: String) {
        val outputFileNameNumHits = "$outputDir/memNumHitsWithRepeats.txt"
        val outputFileNameNumMems = "$outputDir/total_mems.txt"
        bufferedWriter(outputFileNameNumHits).use { numHits ->
            numHits.write("readName\tindex\treadLength\tnumHits\n")
            bufferedWriter(outputFileNameNumMems).use { numMems ->
                numMems.write("readName\tindex\tnumMems\n")
                for (readMems in memsCollected) {
                    val readNameWIndex = readMems.key
                    val index = readNameWIndex.split("_").last()
                    val readName = readNameWIndex.substringBeforeLast("_")
                    numMems.write("${readName}\t${index}\t${readMems.value.size}\n")
                    for( mem in readMems.value) {
                        numHits.write("${readName}\t${index}\t${mem.readEnd - mem.readStart}\t${mem.numHits}\n")
                    }
                }
            }
        }
    }


    fun filterMemsByNumHits(memsCollected: Map<String,List<MEM>>): Map<String,MEM> {
        return memsCollected.map { (readName, mems) ->
            val minHits = mems.minOfOrNull { it.numHits } ?: 0
            val filteredMems = mems.first { it.numHits == minHits }
            Pair(readName, filteredMems)
        }.toMap()
    }

    fun filterMemsByLength(memsCollected: Map<String,List<MEM>>): Map<String,MEM> {
        return memsCollected.map { (readName, mems) ->
            val longestMem = mems.maxByOrNull { it.readEnd - it.readStart } ?: mems.first()
            Pair(readName, longestMem)
        }.toMap()
    }

    fun filterMemsByHitsVsLength(memsCollected: Map<String,List<MEM>>): Map<String,MEM> {
        return memsCollected.map { (readName, mems) ->
            val memsWithScore = mems.map { Pair(it,it.numHits.toDouble() / (it.readEnd - it.readStart)) }
            val minScore = memsWithScore.minOfOrNull { it.second } ?: 0.0
            val filteredMems = memsWithScore.first { it.second == minScore }.first
            Pair(readName, filteredMems)
        }.toMap()
    }

    fun writeOutNumHits(filteredMems: Map<String,MEM>, outputFile: String) {
        bufferedWriter(outputFile).use { writer ->
            writer.write("readName\tindex\treadLength\tnumHits\n")
            for ((readNameWIndex, mems) in filteredMems) {
                val index = readNameWIndex.split("_").last()
                val readName = readNameWIndex.substringBeforeLast("_")
                writer.write("${readName}\t${index}\t${mems.readEnd - mems.readStart}\t${mems.numHits}\n")
            }
        }
    }
}