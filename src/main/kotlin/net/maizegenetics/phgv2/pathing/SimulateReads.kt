package net.maizegenetics.phgv2.pathing

import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.utils.getBufferedReader
import net.maizegenetics.phgv2.utils.getBufferedWriter
import net.maizegenetics.phgv2.utils.retrieveAgcContigs
import net.maizegenetics.phgv2.utils.retrieveAgcData
import org.apache.logging.log4j.LogManager
import java.io.BufferedWriter
import java.io.File
import kotlin.math.abs
import kotlin.random.Random

sealed class SampleList {
    abstract fun sampleChromosomeList(): List<Pair<String,String>>
    data class SimpleSampleList(val listOfSamples: String): SampleList() {
        override fun sampleChromosomeList() = listOfSamples.split(",").map { Pair(it,"") }
    }

    data class SampleChromList(val filename: String) : SampleList() {
        override fun sampleChromosomeList() = getBufferedReader(filename).use { myReader ->
            myReader.readLines().map {
                val sampleChrom = it.split("\t")
                Pair(sampleChrom[0], sampleChrom[1])
            }
        }
    }
}

/**
 * This class simulates reads from sequence stored in an agc data store. The reads are selected at random from each contig.
 * The class provides for two different ways to specify what samples
 *
 */
class SimulateReads : CliktCommand(help = "Simulate reads from agc sequence") {
    private val myLogger = LogManager.getLogger(SimulateReads::class.java)

    val dbPath by option(help = "Path of directory containing tiledb and asseblies.agc").required()

    val hvcfDir by option(help = "The directory containing the h.vcf files for the samples.").required()

    val sampleInput: SampleList by mutuallyExclusiveOptions<SampleList>(
        option("--sample-names", help = "")
            .convert{ SampleList.SimpleSampleList(it) },
        option("--sample-file", help = "")
            .convert{ SampleList.SampleChromList(it) }
    ).single().required()

    val fastqOutDir by option(help = "Output directory for the fastq files.")
        .required()

    val simulateF2 by option(help = "Flag. If present, simulate reads from an F2 with 2 breakpoints per chromosome. " +
                "The default is to simulate reads for a homozygous diploid or haploid.")
        .flag()

    val coverage by option(help = "Proportion of sequence covered by reads. Default = 0.1")
        .double().default(0.1)

    val readLength by option(help = "Single reads of this length will be simulated. Default = 100")
        .int().default(100)

    val chromosomes by option(help = "A list of the chromosomes to be used to generate reads when only sample names are supplied. Default = no list.")
        .default("")

    //agcData is a map of sample name -> chromosome list
    lateinit var agcChromosomeBySample:Map<String, List<String>>

    override fun run() {
        val sampleNameSet = sampleInput.sampleChromosomeList().map { it.first }.toSet()
        val cmdList = mutableListOf("listctg")
        cmdList.addAll(sampleNameSet)
        agcChromosomeBySample = retrieveAgcData(dbPath, cmdList).associate { it ->
            val splitOnColon = it.split(":")
            val genome = splitOnColon[0]
            val chromList = splitOnColon[1].split(",")
            Pair(genome, chromList)
        }
        if (simulateF2) diploidReads() else singleReads()
    }

    /**
     * Expands a SimpleSampleList to a SampleChromosomeList using the chromosome list or chromosome names from agc
     */
    private fun sampleChromosomePairs(sampleNameSet: Set<String>) :  List<Pair<String,String>> {
        val useChromosomeList = chromosomes.isNotBlank()
        return if (useChromosomeList) {
            sampleNameSet.flatMap { sample -> chromosomes.split(",").map { Pair(sample, it.trim()) } }
        } else {
            sampleNameSet.flatMap { sample -> agcChromosomeBySample[sample]!!.map { Pair(sample, it.trim()) } }
        }
    }

    private fun singleReads() {
        sampleToChromosomes().entries.forEach { (sample, chrList) ->  singleReadsFromSample(sample, chrList)}
    }

    private fun sampleToChromosomes(): Map<String, List<String>> {
        val samplePairList = if (sampleInput is SampleList.SimpleSampleList) {
            sampleChromosomePairs(sampleInput.sampleChromosomeList().map { it.first }.toSet())
        } else sampleInput.sampleChromosomeList()

        return samplePairList.groupBy({it.first}, {it.second})
    }

    fun singleReadsFromSample(sampleName: String, chromosomes: List<String>) {
        //read header will be @<sampleName>:<chr>:<start>:<end>
        //then sequence, +, string of "E" same length as sequence
        val outdir = if (fastqOutDir.endsWith("/")) fastqOutDir else "$fastqOutDir/"
        val outputFile = File("$outdir$sampleName-single-reads.fq")

        myLogger.info("Generating single reads for $sampleName, read length = $readLength, coverage = $coverage")
        getBufferedWriter(outputFile).use { fastqWriter ->
            for (chr in chromosomes) {
                val sequence = retrieveAgcContigs(dbPath, listOf("chr${chr.toString()}@$sampleName"))
                val nucseq = sequence.entries.first().value
                val chrlen = nucseq.size()
                val numberOfReads = (chrlen / readLength * coverage).toInt()
                val startLimit = chrlen - readLength
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

    fun diploidReads() {
        val sampleToChromosomeMap = sampleToChromosomes()
        if (sampleInput is SampleList.SimpleSampleList) {

            val samples = sampleInput.sampleChromosomeList().map { it.first }
            //generate reads for each pair of samples
            samples.windowed(2,2).forEach {
                diploidReadsForSamplePair(it, sampleToChromosomeMap)
            }
        } else {
            //same sample could be in more than one pair, so need to process list of Pair(sample, chromosome) in order
            val sampleChromosomeList = sampleInput.sampleChromosomeList()

            var firstSample = ""
            var secondSample = ""

            for (sampleChrPair in sampleChromosomeList) {
                when {
                    firstSample == "" -> firstSample = sampleChrPair.first
                    firstSample == sampleChrPair.first -> {} //do nothing, just skip
                    secondSample == "" -> secondSample = sampleChrPair.first
                    secondSample == sampleChrPair.first -> {} //do nothing, just skip
                    else -> {
                        diploidReadsForSamplePair(listOf(firstSample, secondSample), sampleToChromosomeMap)
                        firstSample = sampleChrPair.first
                        secondSample = ""
                    }
                }
            }

            //The final pair will not have been processed, so process it now.
            if (secondSample != "") diploidReadsForSamplePair(listOf(firstSample, secondSample), sampleToChromosomeMap)

        }
    }

    fun diploidReadsForSamplePair(samples: List<String>, chromosomesBySample: Map<String, List<String>>) {

        val outdir = if (fastqOutDir.endsWith("/")) fastqOutDir else "$fastqOutDir/"
        val f2name = "${samples[0]}X${samples[1]}F2"
        val outputFile = File("$outdir$f2name-single-reads.fq")
        val breakpointFile = File("$outdir$f2name-breakpoints.txt")

        val candidateList = listOf("${samples[0]}.h.vcf", "${samples[0]}.h.vcf.gz","${samples[1]}.h.vcf", "${samples[1]}.h.vcf.gz")
        val hvcfList = File(hvcfDir).listFiles().filter { candidateList.contains(it.name) }.map { it.absolutePath }
        val graph = HaplotypeGraph(hvcfList)
        val rangesByChrom = graph.rangesByContig()

        myLogger.info("Generating single reads for ${samples[0]} X ${samples[1]} F2")

        //generate the F2 name to use in the fastq file
        val headerName = f2name
        getBufferedWriter(breakpointFile).use { breakPointWriter ->
            breakPointWriter.write("index\tchr\tstart\tend\n")
            getBufferedWriter(outputFile).use { fastqWriter ->
                for (chr in (1..10)) {
                    println("generating reads for chromosome $chr")
                    //pick the breakpoints
                    val chrRanges = rangesByChrom["chr$chr"]
                    check(chrRanges != null) { "chrRanges null for chr$chr" }
                    val breakpointRangeArray = Array<List<ReferenceRange>>(2) { breakpointRanges(graph, chrRanges) }
                    for (index in breakpointRangeArray.indices) {
                        for (refrange in breakpointRangeArray[index]) breakPointWriter.write("$index\t${refrange.contig}\t${refrange.start}\t${refrange.end}\n")
                    }

                    //need to generate reads from 2 chromosomes as they are not identical for a heterozygous diploid
                    //for chromosome 1 get sequence for samples[0] before breakpoint1 and after breakpoint2, otherwise samples[1]
                    //for chromosome 2 get sequence for samples[1] before breakpoint1 and after breakpoint2, otherwise samples[0]
                    val sequence = retrieveAgcContigs("/tiledb_maize", listOf("chr$chr@${samples[0]}", "chr$chr@${samples[1]}"))
                    println("nucseqs retrieved from agc:")
                    for (namePair in sequence.keys) println("${namePair.first}, ${namePair.second}")
                    //generate reads
                    for (index in breakpointRangeArray.indices) {
                        for (segment in 1..3) {
                            val sampleName = when (segment) {
                                1 -> if (index == 0) "samples[0]" else "samples[1]"
                                2 -> if (index == 0) "samples[1]" else "samples[0]"
                                3 -> if (index == 0) "samples[0]" else "samples[1]"
                                else -> throw IllegalArgumentException("invalid segment number $segment")
                            }
                            val sampleGamete = SampleGamete(sampleName)
                            val nucseq = sequence[Pair(sampleName, "chr$chr")]!!
                            val endpoints =
                                segmentEndpoints(segment, sampleGamete, breakpointRangeArray[index], graph, nucseq)
                            writeSegmentedReads(
                                endpoints[0],
                                endpoints[1],
                                readLength,
                                nucseq,
                                fastqWriter,
                                headerName,
                                sampleGamete
                            )
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
    fun segmentEndpoints(
        whichSegment: Int,
        sampleGamete: SampleGamete,
        breakPointRanges: List<ReferenceRange>,
        graph: HaplotypeGraph,
        nucseq: NucSeq
    ): IntArray {
        val endpoints = IntArray(2)

        val hapid1 = graph.sampleToHapId(breakPointRanges[0], sampleGamete)
        check(hapid1 != null) { "hapid is null for $sampleGamete at ${breakPointRanges[0]}" }
        val region1 = graph.altHeader(hapid1)
        check(region1 != null) { "region is null for $hapid1, $sampleGamete, ${breakPointRanges[0]}" }

        val hapid2 = graph.sampleToHapId(breakPointRanges[1], sampleGamete)
        check(hapid2 != null) { "hapid is null for $sampleGamete at ${breakPointRanges[1]}" }
        val region2 = graph.altHeader(hapid2)
        check(region2 != null) { "region is null for $hapid2, $sampleGamete, ${breakPointRanges[1]}" }


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
        check(endpoints[0] < endpoints[1]) { "endpoints not ordered correctly for $sampleGamete, ${breakPointRanges[0]}, ${breakPointRanges[1]}" }
        return endpoints
    }

    fun writeSegmentedReads(
        start: Int,
        end: Int,
        readLength: Int,
        nucseq: NucSeq,
        fastqWriter: BufferedWriter,
        headerName: String,
        sampleGamete: SampleGamete
    ) {
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
     * @param graph     a HaplotypeGraph with only the two samples to be used for breakpoints
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
            //make sure that both samples have haplotypes in this range
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

        return breakPointRanges.sortedBy { it.start }
    }

}
