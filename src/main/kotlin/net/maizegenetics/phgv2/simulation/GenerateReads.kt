package net.maizegenetics.phgv2.simulation

import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.pathing.AlignmentUtils
import net.maizegenetics.phgv2.pathing.FindPaths
import net.maizegenetics.phgv2.utils.getBufferedReader
import net.maizegenetics.phgv2.utils.getBufferedWriter
import net.maizegenetics.phgv2.utils.retrieveAgcContigs
import org.apache.logging.log4j.LogManager
import java.io.BufferedWriter
import java.io.File
import kotlin.math.abs
import kotlin.random.Random

class GenerateReads: CliktCommand(help="Generate simulated reads") {
    private val myLogger = LogManager.getLogger(GenerateReads::class.java)

    val sampleNames by option(help = "sample names").default("")

    override fun run() {
//        for (sampleName in sampleNames.split(",")) {
//            singleReadsFromSample(sampleName, "/tiledb_maize")
//        }
//        testReadMapping()

//        makeDiploidReads()

//        singleReadsWithinRangesOnly("/tiledb_maize")

//        decodeReadMapping()
        testHapidsFromGraph()
    }

    fun singleReadsFromSample(sampleName: String, dbPath: String) {
        //read header will be @<sampleName>:<chr>:<start>:<end>
        //then sequence, +, string of "E" same length as sequence
        val readLength = 100
        val outputFile = File("/workdir/simulation/fastq/$sampleName-single-reads.fq")

        myLogger.info("Generating single reads for $sampleName")
        getBufferedWriter(outputFile).use { fastqWriter ->
            for (chr in 1..10) {
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

    fun singleReadsWithinRangesOnly(dbPath: String) {
        //read header will be @<sampleName>:<chr>:<start>:<end>
        //then sequence, +, string of "E" same length as sequence
        val readLength = 100
        val outputFile = File("/workdir/simulation/fastq/B73-single-reads-byrange.fq")
        val graph = HaplotypeGraph(listOf("/workdir/simulation/hvcf_files/B73.h.vcf.gz"))
        val contigToRefranges = graph.rangesByContig()

        myLogger.info("Generating single reads for B73 by range")
        getBufferedWriter(outputFile).use { fastqWriter ->
            for (chr in 1..10) {
                val ranges = contigToRefranges["chr$chr"]!!
                val sequence = retrieveAgcContigs(dbPath, listOf("chr${chr.toString()}@B73"))
                val nucseq = sequence.entries.first().value
//                val chrlen = nucseq.size()
//                val numberOfReads = chrlen / 1000
//                val startLimit = chrlen - 100
                for (range in ranges) {
                    val myLength = range.end - range.start + 1
                    val startLimit = myLength - 100
                    val numberOfReads = myLength / 1000 + 1
                    repeat(numberOfReads) {
                        val start = Random.nextInt(startLimit) + range.start
                        val end = start + readLength - 1
                        val readSeq = nucseq[start..end].seq()
                        fastqWriter.write("@B73.$it\n")
                        fastqWriter.write(readSeq + "\n")
                        fastqWriter.write("+\n")
                        fastqWriter.write("${"E".repeat(readLength)}\n")
                    }

                }
            }
        }

    }

    fun makeDiploidReads() {
        val readLength = 100
        val outputFile = File("/workdir/simulation/fastq/CML247XOh43-single-reads.fq")
        val breakpointFile = File("/workdir/simulation/fastq/CML247XOh43-breakpoints.txt")

        val graph = HaplotypeGraph(listOf("/workdir/simulation/hvcf_files/CML247.h.vcf.gz","/workdir/simulation/hvcf_files/Oh43.h.vcf.gz"))
        val rangesByChrom = graph.rangesByContig()

        myLogger.info("Generating single reads for CML247 X Oh43 F2")
        //generate the F2 name to use in the fastq file
        val headerName = "CML247XOh43F2"
        getBufferedWriter(breakpointFile).use {breakPointWriter ->
            breakPointWriter.write("index\tchr\tstart\tend\n")
            getBufferedWriter(outputFile).use { fastqWriter ->
                for (chr in (1..10)) {
                    println("generating reads for chromosome $chr")
                    //pick the breakpoints
                    val chrRanges = rangesByChrom["chr$chr"]
                    check(chrRanges != null) {"chrRanges null for chr$chr"}
                    val breakpointRangeArray = Array<List<ReferenceRange>>(2) {breakpointRanges(graph, chrRanges)}
                    for (index in breakpointRangeArray.indices) {
                        for (refrange in breakpointRangeArray[index]) breakPointWriter.write("$index\t${refrange.contig}\t${refrange.start}\t${refrange.end}\n")
                    }

                    //need to generate reads from 2 chromosomes as they are not identical for a heterozygous diploid
                    //for chromosome 1 get sequence for CML247 before breakpoint1 and after breakpoint2, otherwise Oh43
                    //for chromosome 2 get sequence for Oh43 before breakpoint1 and after breakpoint2, otherwise CML247
                    val sequence = retrieveAgcContigs("/tiledb_maize", listOf("chr$chr@CML247","chr$chr@Oh43"))
                    println("nucseqs retrieved from agc:")
                    for (namePair in sequence.keys) println("${namePair.first}, ${namePair.second}")
                    //generate reads
                    for (index in breakpointRangeArray.indices) {
                        for (segment in 1..3) {
                            val sampleName = when(segment) {
                                1 -> if (index == 0) "CML247" else "Oh43"
                                2 -> if (index == 0) "Oh43" else "CML247"
                                3 -> if (index == 0) "CML247" else "Oh43"
                                else -> throw IllegalArgumentException("invalid segment number $segment")
                            }
                            val sampleGamete = SampleGamete(sampleName)
                            val nucseq = sequence[Pair(sampleName, "chr$chr")]!!
                            val endpoints = segmentEndpoints(segment, sampleGamete, breakpointRangeArray[index], graph, nucseq)
                            writeSegmentedReads(endpoints[0], endpoints[1], readLength, nucseq, fastqWriter, headerName, sampleGamete)
                        }
                    }

                }
            }

        }

    }

    /**
     * @param whichSegment segment (1, 2, or 3) for which to return endpoints
     * @param sampleGamete the SampleGamete
     * @param breakPointRanges the ranges used for the breakpoints
     * @param graph the HaplotypeGraph
     * @param nucseq the nucseq for the sample gamete, used to get the chromosome length
     * @return the segment endpoints in the sample gamete coordinates
     */
    fun segmentEndpoints(whichSegment: Int, sampleGamete: SampleGamete, breakPointRanges: List<ReferenceRange>, graph: HaplotypeGraph, nucseq: NucSeq): IntArray {
        val endpoints = IntArray(2)

        val hapid1 = graph.sampleToHapId(breakPointRanges[0], sampleGamete)
        check(hapid1 != null) {"hapid is null for $sampleGamete at ${breakPointRanges[0]}"}
        val region1 = graph.altHeader(hapid1)
        check(region1 != null) {"region is null for $hapid1, $sampleGamete, ${breakPointRanges[0]}"}

        val hapid2 = graph.sampleToHapId(breakPointRanges[1], sampleGamete)
        check(hapid2 != null) {"hapid is null for $sampleGamete at ${breakPointRanges[1]}"}
        val region2 = graph.altHeader(hapid2)
        check(region2 != null) {"region is null for $hapid2, $sampleGamete, ${breakPointRanges[1]}"}


        endpoints[0] = when (whichSegment) {
            1 -> 1
            2 -> graph.altHeader(graph.sampleToHapId(breakPointRanges[0], sampleGamete)!!)!!.regions[0].second.position
            3 -> graph.altHeader(graph.sampleToHapId(breakPointRanges[1], sampleGamete)!!)!!.regions[0].second.position
            else -> throw IllegalArgumentException("invalid segment number")
        }
        endpoints[1] = when (whichSegment) {
            1 -> graph.altHeader(graph.sampleToHapId(breakPointRanges[0], sampleGamete)!!)!!.regions[0].second.position
            2 -> graph.altHeader(graph.sampleToHapId(breakPointRanges[1], sampleGamete)!!)!!.regions[0].second.position
            3 -> {
                nucseq.size()
            }
            else -> throw IllegalArgumentException("invalid segment number")
        }
        return endpoints
    }

    fun writeSegmentedReads(start: Int, end: Int, readLength: Int, nucseq: NucSeq, fastqWriter: BufferedWriter, headerName: String, sampleGamete: SampleGamete) {
        val segmentLength = end - start + 1
        val numberOfReads = segmentLength / 2000
        val startLimit = segmentLength - 100

        println("writing reads for $sampleGamete, $start - $end")
        repeat(numberOfReads) {
            val readStart = Random.nextInt(startLimit) + start
            val readEnd = readStart + readLength - 1
            val readSeq = nucseq[readStart..readEnd].seq()
            fastqWriter.write("@$headerName.$it\n")
            fastqWriter.write(readSeq + "\n")
            fastqWriter.write("+\n")
            fastqWriter.write("${"E".repeat(readLength)}\n")
        }
    }

    /**
     * select a pair of breakpoints for a single chromosome.
     * @param graph     a HaplotypeGraph
     * @param ranges    a list of the ReferenceRanges for a single chromosome
     * @return      a list of two ranges selected as breakpoints
     */
    fun breakpointRanges(graph: HaplotypeGraph, ranges: List<ReferenceRange>): List<ReferenceRange> {
        //don't use ranges that are not at least 10 from the end
        val minDistance = 500
        val breakPointRanges = mutableListOf<ReferenceRange>()
        val refrangeToIndexMap = graph.refRangeToIndexMap()
        val requiredSampleGametes = listOf(SampleGamete("CML247"), SampleGamete("Oh43"))


        while (breakPointRanges.size < 2) {
            val candidateRange = ranges.random()
            //make sure that both CML247 and Oh43 have haplotypes in this range
            val hapidSampleMap = graph.hapIdToSampleGametes(candidateRange)
            val sampleGametesInCandidateRange = hapidSampleMap.flatMap { it.value }
            if (!sampleGametesInCandidateRange.containsAll(requiredSampleGametes)) continue

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


        for (refrange in graph.ranges()) {
            if (refrange.contig == "chr1" && refrange.start in 5769996..6162051) {

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

    fun decodeReadMapping() {
        getBufferedWriter("/workdir/simulation/reads/B73-single-reads_decoded.txt").use {myWriter ->
            myWriter.write("contig\tstart\tend\tsamples\tcount\n")
            val fileList = File("/workdir/simulation/hvcf_files").listFiles().filter { it.name.endsWith(".h.vcf.gz") }.map { it.absolutePath }
            val graph = HaplotypeGraph(fileList)
            val readMappings = AlignmentUtils.importReadMapping("/workdir/simulation/reads/B73-single-reads_readMapping.txt")
            val mappingsByRefrange = FindPaths().readMappingByRange(readMappings, graph)

            for ((refrange, mappings) in mappingsByRefrange) {
                val hapidToSample = graph.hapIdToSampleGametes(refrange)
                val hapidCounts = mappings.entries.map {(hapidList, count) -> hapidList.map { Pair(it, count) }}
                    .flatten().groupBy({it.first},{it.second}).mapValues { it.value.sum() }
                val sampleCounts = hapidCounts.mapKeys { hapidToSample[it.key]!! }
                for (entry in sampleCounts) {
                    myWriter.write("${refrange.contig}\t${refrange.start}\t${refrange.end}\t${entry.key}\t${entry.value}\n")
                }
            }
        }

    }

    fun testHapidsFromGraph() {
        val fileList = File("/workdir/simulation/hvcf_files").listFiles().filter { it.name.endsWith(".h.vcf.gz") }.map { it.absolutePath }
        println("file list :")
        fileList.forEach { println(it) }
        val graph = HaplotypeGraph(fileList)
        val b73SampleGamete = SampleGamete("B73")
        for (range in graph.ranges()) {
            val b73hapid = graph.sampleToHapId(range, b73SampleGamete)
            if (b73hapid == null) println("b73hapid is null at $range")
            else {
                val altHeader = graph.altHeader(b73hapid)
                if (altHeader == null) println("altHeader is null for $b73hapid at $range") else {
                    if (altHeader.sampleName() != "B73") {
                        println("For $range and hapid = $b73hapid, altHeader sample name is ${altHeader.sampleName()} and regions = ${altHeader.regions}")
                        val hapidMap = graph.hapIdToSampleGametes(range)
                        println(" hapidMap = $hapidMap")
                    }
                }
            }
        }
    }

    fun testAgcContent() {

    }
}

