package net.maizegenetics.phgv2.simulation

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import htsjdk.samtools.fastq.FastqReader
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.pathing.AlignmentUtils
import net.maizegenetics.phgv2.utils.getBufferedWriter
import net.maizegenetics.phgv2.utils.retrieveAgcContigs
import org.apache.logging.log4j.LogManager
import java.io.File
import kotlin.math.abs
import kotlin.math.max
import kotlin.random.Random

class GenerateReads: CliktCommand(help="Generate simulated reads") {
    private val myLogger = LogManager.getLogger(GenerateReads::class.java)

    val sampleNames by option(help = "sample names").default("")

    override fun run() {
//        for (sampleName in sampleNames.split(",")) {
//            singleReadsFromSample(sampleName, "/tiledb_maize")
//        }
        testReadMapping()
    }

    fun singleReadsFromSample(sampleName: String, dbPath: String) {
        //read header will be @<sampleName>:<chr>:<start>:<end>
        //then sequence, +, string of "E" same length as sequence
        val readLength = 100
        val outputFile = File("/workdir/simulation/fastq/$sampleName-single-reads.fq")

        myLogger.info("Generating single reads for $sampleName")
        getBufferedWriter(outputFile).use { fastqWriter ->
            for (chr in listOf(10,9,8,7,6,5,4,3,2,1)) {
                val sequence = retrieveAgcContigs(dbPath, listOf("chr${chr.toString()}@$sampleName"))
                val nucseq = sequence.entries.first().value
                val chrlen = nucseq.size()
                val numberOfReads = chrlen / 1000
                val startLimit = chrlen - 100
                repeat(numberOfReads) {
                    val start = Random.nextInt(startLimit)
                    val end = start + readLength - 1
                    val readSeq = nucseq[start..end].seq()
                    fastqWriter.write("@$sampleName.$it\n")
                    fastqWriter.write(readSeq + "\n")
                    fastqWriter.write("+\n")
                    fastqWriter.write("${"E".repeat(readLength)}\n")
                }
            }
        }

    }

    fun makeDiploidReads() {
        val readLength = 100
        val outputFile = File("/workdir/simulation/fastq/CML247XOh43-single-reads.fq")
        val breakpointFile = File("/workdir/simulation/fastq/CML247XOh43-breakpoints.txt")

        val graph = HaplotypeGraph(listOf("/workdir/simulation/hvcf_files/CML247.h.vcf.gz","/workdir/simulation/hvcf_files/Oh43.h.vcf.gz"))
        val ranges = graph.ranges()

        myLogger.info("Generating single reads for CML247 X Oh43")
        //generate the F2 name to use in the fastq file
        val headerName = "CML247XOh43F2"
        getBufferedWriter(breakpointFile).use {breakPointWriter ->
            getBufferedWriter(outputFile).use { fastqWriter ->
                for (chr in (1..10)) {
                    //pick the breakpoints
                    val chrRanges = ranges.filter { it.contig =="chr$chr" }
                    val breakpointRanges1 = breakpointRanges(graph, chrRanges)
                    val breakpointRanges2 = breakpointRanges(graph, chrRanges)

                    //need to generate reads from 2 chromosomes as they are not identical for a heterozygous diploid
                    //for chromosome 1 get sequence for CML247 before breakpoint1 and after breakpoint3, otherwise Oh43
                    //for chromosome 2 get sequence for Oh43 before breakpoint1 and after breakpoint3, otherwise CML247
                    val sequence = retrieveAgcContigs("/tiledb_maize", listOf("chr$chr@CML247, chr$chr@Oh43"))
                    val nucseq = sequence.entries.first().value
                    val chrlen = nucseq.size()
                    val numberOfReads = chrlen / 1000
                    val startLimit = chrlen - 100
                    repeat(numberOfReads) {
                        val start = Random.nextInt(startLimit)
                        val end = start + readLength - 1
                        val readSeq = nucseq[start..end].seq()
                        fastqWriter.write("@$headerName.$it\n")
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
        val minDistance = 500
        val breakPointRanges = mutableListOf<ReferenceRange>()
        val refrangeToIndexMap = graph.refRangeToIndexMap()

        while (breakPointRanges.size < 2) {
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

            if (allRegionsOK && breakPointRanges.size == 0) breakPointRanges.add(candidateRange)
            else if (allRegionsOK) {
                //only use if range indices are different by minDistance or more
                val existingRangeIndex = refrangeToIndexMap[breakPointRanges[0]]!!
                val candidateRangeIndex = refrangeToIndexMap[candidateRange]!!
                if (abs(candidateRangeIndex - existingRangeIndex) >= minDistance) breakPointRanges.add(candidateRange)
            }

        }

        return breakPointRanges
    }

    fun testReadMapping() {
        val minProportionOfMaxCount = 1.0
        val minSameReferenceRange = 0.9

        val hvcfFilenames = File("/workdir/simulation/hvcf_files").listFiles()
            .filter { it.name.endsWith(".h.vcf.gz") }.map { it.absolutePath }
        val graph = HaplotypeGraph(hvcfFilenames)
        val kmerIndexMap = AlignmentUtils.loadKmerMaps("/workdir/simulation/hvcf_files/kmerIndex.txt", graph)
//        val fastqFiles = Pair(keyFileRecord.file1, keyFileRecord.file2)


//        net.maizegenetics.phgv2.pathing.myLogger.info("reading records from the fastq file(s): ${fastqFiles.first}, ${fastqFiles.second}")


        val rangeIdToBitsetMap =
            AlignmentUtils.convertRefRangeToIdBitsetMap(kmerIndexMap.rangeToBitSetMap, graph.refRangeToIndexMap())

        //Setting up the channels for the coroutines
//        val readStringChannel = Channel<Pair<String, String>>(100)
//        val sortedHapidListChannel = Channel<List<String>>(100)
        val numberOfMappingThreads = 10

        val mySequence = retrieveAgcContigs("/tiledb_maize", listOf("chr1@B73"))

        val hapidToRange = graph.hapIdToRefRangeMap()
        val b73SampleGamete = SampleGamete("B73")

        //Run the blocking coroutine to read the fastq files and process the reads
        for (refrange in graph.ranges()) {
            if (refrange.contig == "chr1") {

                val rangeSequence = mySequence.get(Pair("B73", "chr1"))!!.get(refrange.start..refrange.end).seq()
                val seq2 = ""
                val rangeHapidSet = mutableSetOf<String>()
                var readsMapped = 0
                for (seq in rangeSequence.windowed(100, 100)) {
                    val hapids = AlignmentUtils.readToHapidSet(
                        seq, minProportionOfMaxCount, minSameReferenceRange,
                        kmerIndexMap.kmerHashToLongMap, rangeIdToBitsetMap, graph.refRangeIdToHapIdMap()
                    )
                    if (hapids.isNotEmpty()) {
                        rangeHapidSet.addAll(hapids)
                        readsMapped++
                    }

                }

                //test whether hapids come from the correct range and are B73 hapids

                val hapidToSampleGameteMap = graph.hapIdToSampleGametes(refrange)
                val rangeHapids = rangeHapidSet.filter { hapidToRange[it]!!.contains(refrange) }
                val nonRangeHapids = rangeHapidSet.filter { !hapidToRange[it]!!.contains(refrange) }
                val rangeSampleGametes = rangeHapids.mapNotNull { hapidToSampleGameteMap[it] }
                val otherRanges = nonRangeHapids.map { hapidToRange[it] }

                println("$refrange - $readsMapped reads mapped, samples: $rangeSampleGametes")
                println("$refrange - other: $otherRanges")

            }


        }

    }
}

