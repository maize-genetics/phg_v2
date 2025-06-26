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

    val gbsSampleName by option("--gbs-sample-name", help = "Sample name")
        .required()

    val memFile by option(help = "Path to the MEM alignment file")
        .required()

    val outputFile by option(help = "Output Directory")
        .required()

    val maxNumHits by option(help = "Max Hit Count")
        .int()
        .default(4)

    override fun run() {
        // This command is a placeholder for reporting MEM alignments.
        // The actual implementation would involve reading MEM alignment files
        // and generating a report based on the alignments.

        println("Processing MEM alignments from file: $memFile")
        reportMemAlignments(memFile, outputFile,gbsSampleName, maxNumHits)
    }

    fun reportMemAlignments(memFile: String, outputFile: String, gbsSampleName : String, maxNumHits: Int = 4) {
        val mems = bufferedReader(memFile).readLines().map { line -> RopeBWTUtils.parseStringIntoMem(line) }

        println("Total Reads Into Mems: ${mems.size}")

        val memsFilteredByNumHits = mems.filter { it.numHits <= maxNumHits }

        println("NumberOfMems <= ${maxNumHits} hits: ${memsFilteredByNumHits.size}")

        processMemsToUnitigs(memsFilteredByNumHits, outputFile, gbsSampleName)
    }

    fun processMemsToUnitigs(memsFilteredByNumHits: List<MEM>, outputFile: String, gbsSampleName: String) {
        bufferedWriter(outputFile).use { writer ->
            writer.write("gbsSampleName\tassemblySampleName\tUnitigId\tcount\n")
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
                writer.write("${gbsSampleName}\t${sampleName}\t${unitigName}\t${count}\n")
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