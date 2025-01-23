package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.utils.getBufferedReader
import net.maizegenetics.phgv2.utils.getBufferedWriter
import org.apache.logging.log4j.LogManager
import java.io.File
import kotlin.time.DurationUnit
import kotlin.time.measureTimedValue

class KmerIndexStatistics : CliktCommand(help="Write kmer counts by reference range")  {
    private val myLogger = LogManager.getLogger(KmerIndexStatistics::class.java)

    //Todo: add dbPath and make hvcfDir optional once java interface is operational

    val hvcfDir by option("--hvcf-dir", help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB. Required")
        .required()

    val indexFile by option(help = "The full path of the kmer index file. Default = <hvcf-dir>/kmerIndex.txt")
        .default("")

    val outputFile by option(help = "The file (with path) to which the output will be written. Default is <hvcf-dir>/kmerIndexCounts.txt")
        .default("")

    override fun run() {
        logCommand(this)

        val graph = buildHaplotypeGraph()
        val indexFilename = if (indexFile.isNotBlank()) indexFile else File(hvcfDir).resolve("kmerIndex.txt").absolutePath
        val outputFilename = if (outputFile.isNotBlank()) outputFile else File(hvcfDir).resolve("kmerIndexCounts.txt").absolutePath
        val rangeCounts = countKmersByRefrange(indexFilename)

        //write the resulting counts to a file
        getBufferedWriter(outputFilename).use { myWriter ->
            myWriter.write("contig\tstart\tend\tlength\tkmerCount\n")
            for (range in graph.ranges()) {
                val rangeLength = range.end - range.start + 1
                val rangeCount = rangeCounts.getOrElse(range) {0}
                myWriter.write("${range.contig}\t${range.start}\t${range.end}\t$rangeLength\t$rangeCount\n")
            }
        }
    }

    private fun buildHaplotypeGraph(): HaplotypeGraph {
        val timedValue = measureTimedValue {
            if(hvcfDir != "") {
                HaplotypeGraph(hvcfDir)
            }
            else {
                //Load in the TileDB
                TODO("TileDB VCF Reader Not implemented yet.  Please run with --hvcf-dir")
            }
        }

        myLogger.info("HaplotypeGraph built in ${timedValue.duration.toDouble(DurationUnit.MILLISECONDS)} ms.")
        return timedValue.value
    }

    /**
     * Similar to AlignmentUtils.loadKmerMaps but only counts the number of kmers in each reference range
     */
    private fun countKmersByRefrange(filename: String): Map<ReferenceRange, Int> {
        val kmerCounts = mutableMapOf<ReferenceRange, Int>()
        var lineCount = 0
        var totalLineCount = 0
        var refrange = ReferenceRange("NULL", 0, 0)

        getBufferedReader(filename).useLines {
            // Read the first line, then contine
            it.forEach { inputStr ->
                totalLineCount++
                if (inputStr.startsWith("GraphHash")) return@forEach

                lineCount = when (inputStr.first()) {
                    '>' -> 1
                    else -> lineCount + 1
                }

                when (lineCount) {
                    1 -> {
                        //line 1 is the range
                        refrange = ReferenceRange.parse(inputStr.substring(1))
                    }

                    2 -> {}
                    3 -> {
                        val hashCount = inputStr.split(",").count { hashStr -> hashStr.contains("@") }
                        kmerCounts[refrange] = hashCount
                    }

                    else -> {
                        throw IllegalArgumentException("Line count = $lineCount in kmerMap file format at line $totalLineCount in file $filename")
                    }
                }
            }
        }

        return kmerCounts
    }


}