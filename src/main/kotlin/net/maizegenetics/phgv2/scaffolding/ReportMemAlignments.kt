package net.maizegenetics.phgv2.scaffolding

import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.pathing.ropebwt.MEM
import net.maizegenetics.phgv2.pathing.ropebwt.RopeBWTUtils
import java.io.File

class ReportMemAlignments : CliktCommand(help = "Report the MEM alignments") {

    val gbsSampleName by option( help = "Sample name")
        .required()

    val memFile by option(help = "Path to the MEM alignment file")
        .required()

    val outputFile by option(help = "Output Directory")
        .required()

    val maxNumHits by option(help = "Max Hit Count")
        .int()
        .default(4)

    val infoFileDir by option(help = "Path to the info file directory")
        .required()

    override fun run() {
        // This command is a placeholder for reporting MEM alignments.
        // The actual implementation would involve reading MEM alignment files
        // and generating a report based on the alignments.

        println("Processing MEM alignments from file: $memFile")
        reportMemAlignments(memFile, outputFile,gbsSampleName, maxNumHits, infoFileDir)
    }

    fun reportMemAlignments(memFile: String, outputFile: String, gbsSampleName : String, maxNumHits: Int = 4, infoFileDirectory: String) {
        val mems = bufferedReader(memFile).readLines().map { line -> RopeBWTUtils.parseStringIntoMem(line) }

        println("Total Reads Into Mems: ${mems.size}")

        val memsFilteredByNumHits = mems.filter { it.numHits <= maxNumHits }

        println("NumberOfMems <= ${maxNumHits} hits: ${memsFilteredByNumHits.size}")

        val assemblyUnitigMetadata = buildUnitigMetadataMap(infoFileDirectory)

        processMemsToUnitigs(memsFilteredByNumHits, outputFile, gbsSampleName, assemblyUnitigMetadata)
    }

    fun buildUnitigMetadataMap(infoFileDirectory: String): Map<Pair<String,String>, Triple<String,String,String>> {
        val unitigMetadata = mutableMapOf<Pair<String,String>, Triple<String,String,String>>()
        File(infoFileDirectory).walk().filter { it.isFile }.forEach { file ->

            val currentFileMap = bufferedReader(file.path).readLines()
                .map { it.split(" ") }
                .associate { Pair(it[4], it[0]) to Triple(it[1], it[2], it[3]) }

            unitigMetadata.putAll(currentFileMap)
        }
        return unitigMetadata
    }

    fun processMemsToUnitigs(memsFilteredByNumHits: List<MEM>, outputFile: String, gbsSampleName: String, assemblyUnitigMetadata: Map<Pair<String,String>, Triple<String,String,String>>) {
        bufferedWriter(outputFile).use { writer ->
            writer.write("gbsSampleName\tassemblySampleName\tUnitigId\trefChr\trefPos\tunitigLength\tcount\n")
            var numDup = 0
            val memsRemoveDups = memsFilteredByNumHits.map {
                val memHits = it.listMemHits
                val memHitsNoDups = memHits.distinctBy { hit -> hit.contig }

                if(memHits.size != memHitsNoDups.size) {
                    numDup++
                }
                memHitsNoDups
            }
            println("Number of MEMs with duplicate hits removed: ${numDup}")

            val counts = memsRemoveDups.flatMap{ it }
                .map { it.contig }.groupingBy { it }.eachCount()

            val sampleCounter = mutableMapOf<String, Int>()
            val sampleUnitigCounter = mutableMapOf<String, Int>()

            for ((unitigId, count) in counts) {
                val unitigName= unitigId.substringBeforeLast("_")
                val sampleName = unitigId.substringAfterLast("_")
                val (refChr, refPos, unitigLength) = assemblyUnitigMetadata[Pair(sampleName, unitigName)]
                    ?: throw IllegalArgumentException("No metadata found for unitig $unitigId in assembly metadata.")
                writer.write("${gbsSampleName}\t${sampleName}\t${unitigName}\t${refChr}\t${refPos}\t${unitigLength}\t${count}\n")
                if (!sampleCounter.containsKey(sampleName)) {
                    sampleCounter[sampleName] = 0
                }
                sampleCounter[sampleName] = sampleCounter[sampleName]!! + count

                if (!sampleUnitigCounter.containsKey(sampleName)) {
                    sampleUnitigCounter[sampleName] = 0
                }
                sampleUnitigCounter[sampleName] = sampleUnitigCounter[sampleName]!! + 1

            }

            println("Total Unitigs: ${counts.size}")
            println("Sample Counts:")
            for ((sampleName, count) in sampleCounter) {
                println("$sampleName: $count")
            }

            println("Sample Unitig Counts:")
            for ((sampleName, count) in sampleUnitigCounter) {
                println("$sampleName: $count")
            }


        }
    }


}