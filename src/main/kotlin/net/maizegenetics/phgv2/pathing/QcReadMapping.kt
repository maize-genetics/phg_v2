package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import htsjdk.samtools.fastq.FastqReader
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.pathing.AlignmentUtils.Companion.convertRefRangeToIdBitsetMap
import net.maizegenetics.phgv2.pathing.AlignmentUtils.Companion.rangeHapidMapFromKmerHash
import kotlin.math.min

/**
 * Function that will take a pair of fastqs and will extract out the kmers and map those to the haplotypes in the graph.
 * This will output a pair of files where we show each kmer and how well they map to the haplotypes.
 */
class QcReadMapping  : CliktCommand(help="Error check reads against a kmer index") {
    val hvcfDir by option(help = "Directory containing hvcf files used to build the HaplotypeGraph.")
        .required()

    val kmerIndex by option(help = "Kmer index file created by build-kmer-index. Default is <hvcfDir>/kmerIndex.txt.")
        .default("")

    val outputDir by option(help = "Output folder for the read mapping results.")
        .required()

    val readFiles by option(help = "Comma separated list of read files to map.")
        .required()

    val numReads by option(help = "Number of reads to process.")
        .int()
        .default(20)


    override fun run() {
        logCommand(this)

        val graph = HaplotypeGraph(hvcfDir)
        println("Loaded Graph")
        //Load up the index file
        val kmerIndexMap = AlignmentUtils.loadKmerMaps(kmerIndex, graph)

        val readNameToCount = mutableMapOf<String, Int>()

        for(file in readFiles.split(",")) {
            FastqReader(bufferedReader(file)).use { reader ->
                var counter = 0
                while(reader.hasNext()) {
                    if(counter < numReads) {
                        val read = reader.next()
                        val readName = read.readName.replace(":", "_")
                        val currentReadNameCount = readNameToCount.getOrDefault(readName, 0)
                        processReads(graph, kmerIndexMap, read.readString, "${outputDir}/${readName}_${currentReadNameCount+1}.txt")
                        readNameToCount[readName] = currentReadNameCount + 1
                        counter++
                    }
                    else {
                        break
                    }
                }

            }
        }

    }

    private fun processReads(graph: HaplotypeGraph, kmerIndex: KmerMapData, read: String, outputFile: String) {
        val kmerHashOffsetMap = kmerIndex.kmerHashToLongMap
        val refrangeToBitSet = convertRefRangeToIdBitsetMap(kmerIndex.rangeToBitSetMap, graph.refRangeToIndexMap())
        val rangeToHapidIndexMap = graph.refRangeIdToHapIdMap()
        val rangeToHapidMap = mutableMapOf<Int, MutableList<String>>()
        val rangeIdToRefRange = graph.refRangeToIndexMap().map { it.value to it.key }.toMap()

        bufferedWriter(outputFile).use { writer ->
            writer.write("Kmer\tHash1\tHash2\tRefRanges\tHapids\n")

            val splitList = read
                .split("[^ACGT]+".toRegex())
                .filter { it.length > 31 }

            for (sequence in splitList) {
                var previousHash = Pair(0L, 0L)

                //for first 31 nucleotides just update the hash
                for (nucleotide in sequence.subSequence(0..30)) {
                    previousHash = BuildKmerIndex.updateKmerHashAndReverseCompliment(previousHash, nucleotide)
                }
                var counter = 0
                for (nucleotide in sequence.subSequence(31 until sequence.length)) {
                    previousHash = BuildKmerIndex.updateKmerHashAndReverseCompliment(previousHash, nucleotide)
                    val minHash = min(previousHash.first.toULong(), previousHash.second.toULong()).toLong()
                    val hapidsMatched =
                        rangeHapidMapFromKmerHash(minHash, kmerHashOffsetMap, refrangeToBitSet, rangeToHapidIndexMap)

                    val currentSeq = sequence.substring(counter, counter + 32)
                    counter++

                    val refRangesHit = hapidsMatched.keys.map { rangeIdToRefRange[it]!! }.map { "${it.contig}:${it.start}-${it.end}" }

                    val rangesHitString = if (refRangesHit.isEmpty()) "None" else refRangesHit.joinToString(",")

                    writer.write("$currentSeq\t${previousHash.first.toULong()}\t${previousHash.second.toULong()}\t${rangesHitString}\t${hapidsMatched}\n")

                    for (entry in hapidsMatched) {
                        val hapidList = rangeToHapidMap[entry.key]
                        if (hapidList == null) rangeToHapidMap.put(entry.key, entry.value)
                        else hapidList.addAll(entry.value)
                    }
                }
            }
        }
    }

}