package net.maizegenetics.phgv2.cli

import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.loadRanges
import net.maizegenetics.phgv2.utils.seqFromAGC
import org.apache.logging.log4j.LogManager
import java.io.File

class RefRangeFasta : CliktCommand(help = "Write Reference Range Sequences to Fasta.") {

    private val myLogger = LogManager.getLogger(RefRangeFasta::class.java)

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        .default("")

    val inputDir by option(help = "Full path to input HVCF file directory")
        .required()

    val outputFile by option(help = "Full path to base output Fasta file. The range will be appended to this.")
        .required()

    val rangeBedfile by option(help = "Full path to bedfile specifying the ranges to export")
        .required()

    override fun run() {

        val inputFiles = File(inputDir)
            .walk()
            .filter {
                it.isFile && (it.name.endsWith(".h.vcf") || it.name.endsWith(".h.vcf.gz") ||
                        it.name.endsWith(".hvcf") || it.name.endsWith(".hvcf.gz"))
            }
            .map { it.absolutePath }
            .toList()

        val graph = HaplotypeGraph(inputFiles)

        val ranges = loadRanges(rangeBedfile)

        ranges.forEach { range ->
            writeFasta(graph, range)
        }

    }

    private fun writeFasta(graph: HaplotypeGraph, range: Pair<Position, Position>) {

        val rangeStr = "${range.first.contig}_${range.first.position}-${range.second.position}"

        val filename = "${outputFile}-${rangeStr}.fasta"

        val sampleToHapid = graph.sampleGameteToHaplotypeId(
            ReferenceRange(
                range.first.contig,
                range.first.position,
                range.second.position
            )
        )

        println("sampleToHapid: $sampleToHapid")

        if (sampleToHapid.isEmpty()) {
            myLogger.warn("No haplotypes found for range: $rangeStr")
            return
        }

        bufferedWriter(filename).use { writer ->

            sampleToHapid
                .filter { it.key.gameteId == 0 }
                .map { it.key.name to it.value }
                .forEach { (sample, hapid) ->

                    val seq: String
                    val displayRanges: List<String>
                    try {
                        val temp = seqFromAGC(dbPath, graph, hapid, range)
                        seq = temp.first
                        displayRanges = temp.second
                    } catch (e: Exception) {
                        myLogger.error("Error getting sequence for $sample-$rangeStr")
                        e.printStackTrace()
                        return
                    }

                    writer.write("> ${displayRanges.joinToString(";")}\n")
                    seq.chunked(60)
                        .forEach { chunk ->
                            writer.write(chunk)
                            writer.newLine()
                        }
                    writer.write("\n")
                }

        }

    }

}