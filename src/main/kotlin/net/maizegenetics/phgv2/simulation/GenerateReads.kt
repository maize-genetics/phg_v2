package net.maizegenetics.phgv2.simulation

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.getBufferedWriter
import net.maizegenetics.phgv2.utils.retrieveAgcContigs
import org.apache.logging.log4j.LogManager
import java.io.File
import kotlin.random.Random

class GenerateReads: CliktCommand(help="Generate simulated reads") {
    private val myLogger = LogManager.getLogger(GenerateReads::class.java)

    val sampleNames by option(help = "sample names").default("")

    override fun run() {
        for (sampleName in sampleNames.split(",")) {
            singleReadsFromSample(sampleName, "/tiledb_maize")
        }
    }

    fun singleReadsFromSample(sampleName: String, dbPath: String) {
        //read header will be @<sampleName>:<chr>:<start>:<end>
        //then sequence, +, string of "E" same length as sequence
        val readLength = 100
        val outputFile = File("/workdir/simulation/reads/$sampleName-single-reads.fq")

        myLogger.info("Generating single reads for $sampleName")
        getBufferedWriter(outputFile).use { fastqWriter ->
            for (chr in (1..10)) {
                val sequence = retrieveAgcContigs(dbPath, listOf("chr${chr.toString()}@$sampleName"))
                val nucseq = sequence.entries.first().value
                val chrlen = nucseq.size()
                val numberOfReads = chrlen / 1000
                val startLimit = chrlen - 100
                repeat(numberOfReads) {
                    val start = Random.nextInt(startLimit)
                    val end = start + readLength - 1
                    val readSeq = nucseq[start..end].seq()
                    fastqWriter.write(readSeq + "\n")
                    fastqWriter.write("+\n")
                    fastqWriter.write("${"E".repeat(readLength)}\n")
                }
            }
        }

    }

    fun makeDiploidReads() {
        val readLength = 100
        val outputFile = File("/workdir/simulation/reads/CML247XOh43-single-reads.fq")
        val breakpointFile = File("/workdir/simulation/reads/CML247XOh43-breakpoints.txt")

        val graph = HaplotypeGraph(listOf("/workdir/simulation/hvcf_files/CML247.h.vcf.gz","/workdir/simulation/hvcf_files/Oh43.h.vcf.gz"))
        val ranges = graph.ranges()

        myLogger.info("Generating single reads for CML247 X Oh43")
        getBufferedWriter(breakpointFile).use {breakPointWriter ->
            getBufferedWriter(outputFile).use { fastqWriter ->
                for (chr in (1..10)) {
                    //pick the breakpoints
                    val chrRanges = ranges.filter { it.contig =="chr$chr" }
                    var breakpointRange1 = chrRanges.random()


                    val sequence = retrieveAgcContigs("/tiledb_maize", listOf("chr$chr@CML247, chr$chr@Oh43"))
                    val nucseq = sequence.entries.first().value
                    val chrlen = nucseq.size()
                    val numberOfReads = chrlen / 1000
                    val startLimit = chrlen - 100
                    repeat(numberOfReads) {
                        val start = Random.nextInt(startLimit)
                        val end = start + readLength - 1
                        val readSeq = nucseq[start..end].seq()
                        fastqWriter.write(readSeq + "\n")
                        fastqWriter.write("+\n")
                        fastqWriter.write("${"E".repeat(readLength)}\n")
                    }
                }
            }

        }

    }

    fun breakpointRanges(graph: HaplotypeGraph, ranges: List<ReferenceRange>): List<ReferenceRange> {
        //don't use ranges that are not at least 10 from the end
        var badRange = true
        while (badRange) {
            val candidateRange = ranges.random()
            val hapidSampleMap = graph.hapIdToSampleGametes(candidateRange)
            var allRegionsOK = true
            for (hapid in hapidSampleMap.keys) {
                val regions = graph.altHeader(hapid)!!.regions
                val isRegionOK = regions.size == 1 && //only 1 region maps to this range
                        regions[0].first.contig == candidateRange.contig &&  //the mapped sequence is from the same range
                        regions[0].first.position < regions[0].second.position  //the mapped sequence is not an inversion
                allRegionsOK = allRegionsOK && isRegionOK
            }
        }

        return listOf()
    }

}

