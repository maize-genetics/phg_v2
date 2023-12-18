package net.maizegenetics.phgv2.pathing

import it.unimi.dsi.fastutil.longs.Long2LongOpenHashMap
import kotlinx.coroutines.channels.ReceiveChannel
import kotlinx.coroutines.channels.SendChannel
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.getBufferedWriter
import net.maizegenetics.phgv2.utils.getBufferedReader
import org.apache.logging.log4j.LogManager
import java.io.*
import java.util.*
import kotlin.math.ceil
import kotlin.math.min

/**
 * data class to hold the KmerMap information.  From the PHGv1 source code:
 */
//data class KmerMapData(val haplotypeListId: Int, val rangeToBitSetMap: Map<Int, BitSet>, val kmerHashToLongMap: Long2LongOpenHashMap)
data class KmerMapData(val rangeToBitSetMap: Map<ReferenceRange, BitSet>, val kmerHashToLongMap: Long2LongOpenHashMap)

data class KeyFileData(val sampleName: String, val file1: String, val file2: String = "")

private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.utils.AlignmentUtils")


fun alignReadsToHaplotypes(hvcfDir: String, kmerIndexFile:String, keyFileRecords:List<KeyFileData>, paired:Boolean, outputDir:String) {
    //loop through all files in hvcfDir and create a list of hvcf files
    val hvcfFiles = File(hvcfDir).walkTopDown().filter { it.isFile }.filter { it.extension == "h.vcf" }.map { "${it.path}/${it.name}" }.toList()

    //create a HaplotypeGraph from the list of hvcf files
    val graph = HaplotypeGraph(hvcfFiles)

    val kmerIndexMap = loadKmerMaps(kmerIndexFile, graph)


}

/**
 * Loads a kmer hash map file to memory for use in mapping reads.
 */
fun loadKmerMaps(filename: String, graph: HaplotypeGraph): KmerMapData {
    //Load the contents of the file into
    //rangeHapidMap: a map of refRange to an Array<String> of the haplotype ids in the ReferenceRange
    // and the BitSet of all hapid sets
    //kmerHashmap: a map of kmer hash to reference range and offset into its BitSet, encoded as a long
    //These data structures all the reference range and haplotype set to be looked up for a kmer has


    val rangeToBitSetMap = mutableMapOf<ReferenceRange, BitSet>()
    val kmerHashMap = Long2LongOpenHashMap()
    var lineCount = 0
    var totalLineCount = 0
    var refrange = ReferenceRange("NULL", 0, 0)

    val refRangeToIdMap = getReferenceRangeToIndexMap(graph)
    getBufferedReader(filename).useLines {
        it.forEach { inputStr ->
            totalLineCount++

            lineCount = when(inputStr.first()) {
                '>' -> 1
                else -> lineCount + 1
            }

            when (lineCount) {
                1 -> {
                    //line 1 is the range
                    refrange = ReferenceRange.parse(inputStr.substring(1))
                }
                2-> {
                    try {
                        val parsedLine = inputStr.split(",")
                        val myBitset = BitSet.valueOf(parsedLine.map { it.toLong() }.toLongArray())
                        rangeToBitSetMap[refrange] = myBitset
                    } catch(e: Exception) {
                        println("error at line $totalLineCount for input = $inputStr")
                        throw java.lang.IllegalArgumentException(e)
                    }
                }
                3 -> {
                    val refrangeLong = (refRangeToIdMap[refrange]!!.toLong()) shl 32
                    val parsedLine = inputStr.split(",")
                    //values are pairs of hash, offset (long,int)
                    //add the pairs of entries to kmerHashMap
                    for (datapair in parsedLine) {
                        val hashOffset = datapair.split("@")

                        if (hashOffset.size < 2) {
                            println("improperly formatted datapair at line $totalLineCount")
                        }

                        val hash = hashOffset[0].toLong()
                        val offset = refrangeLong or hashOffset[1].toLong()
                        kmerHashMap.put(hash, offset)
                    }
                }
                else -> {
                    throw IllegalArgumentException("Line count = $lineCount in kmerMap file format at line $totalLineCount in file $filename")
                }
            }
        }
    }

    return KmerMapData(rangeToBitSetMap, kmerHashMap)
}

/**
 * Creates a map of ReferenceRange -> (map of hapid -> index) for a [HaplotypeGraph]
 */
fun getRefRangeToHapidMap(graph: HaplotypeGraph) : Map<ReferenceRange,Map<String,Int>>{
    //This creates a map of ReferenceRangeId -> (map of hapid -> index)
    return graph.ranges().associateWith { range ->
        graph.hapIdToSamples(range).keys.toSortedSet()
            .mapIndexed { index, hapid -> hapid to index }.toMap()
    }
}

/**
 * Creates a map of ReferenceRange -> index for a [HaplotypeGraph]. Because the ranges are sorted by the
 * HaplotypeGraph method range(), a graph always returns the same map.
 */
fun getReferenceRangeToIndexMap(graph: HaplotypeGraph) : Map<ReferenceRange,Int> {
    return graph.ranges().mapIndexed { index, range -> range to index }.toMap()
}


/**
 * Function to process a read or a readPair into their corresponding hapIdSet.  This will work with both single and paired end reads.
 * To Work with single ended reads, the second String in the Pair should be an empty String.
 */
suspend fun processReads(reads: ReceiveChannel<Pair<String, String>>,
                               sortedLists: SendChannel<List<Int>>,
                               minProportionOfMaxCount: Double = 1.0,
                               minSameReferenceRange: Double = 0.9,
                               kmerHashOffsetMap: Long2LongOpenHashMap,
                               refrangeToBitSet: Map<Int, BitSet>,
                               rangeToHapidIndexMap: Map<Int, Map<Int, Int>>) {
    for (pairOfReads in reads) {
        val result1 = readToHapidSet(pairOfReads.first, minProportionOfMaxCount, minSameReferenceRange, kmerHashOffsetMap, refrangeToBitSet, rangeToHapidIndexMap)
        val result2 = if(pairOfReads.second.isNotEmpty()) {
            readToHapidSet(pairOfReads.second, minProportionOfMaxCount, minSameReferenceRange, kmerHashOffsetMap, refrangeToBitSet, rangeToHapidIndexMap)
        }
        else {
            setOf<Int>()
        }

        val intersectResult = if(result2.isNotEmpty()) {
            result1.intersect(result2)
        }
        else {
            result1
        }

        if (intersectResult.isNotEmpty()) sortedLists.send(intersectResult.sorted())

    }
}

/**
 * Takes a [read] and generates a list of hapids to which its kmers map. Returns only hapids from a single
 * reference range. If no kmers map or if fewer than [minSameReferenceRange] of the kmers map to a single
 * reference range an empty Set will be returned.
 */
private fun readToHapidSet(read: String,
                           minProportionOfMaxCount: Double = 1.0,
                           minSameReferenceRange: Double = 0.9,
                           kmerHashOffsetMap: Long2LongOpenHashMap,
                           refrangeToBitSet: Map<Int, BitSet>,
                           rangeToHapidIndexMap: Map<Int, Map<Int, Int>>): Set<Int> {
    //generate kmer hash from the read
    val rangeToHapidMap = mutableMapOf<Int, MutableList<Int>>()
    val splitList = read
        .split("[^ACGT]+".toRegex())
        .filter { it.length > 31 }
    for (sequence in splitList) {
        var previousHash = Pair(0L, 0L)

        //for first 31 nucleotides just update the hash
        for (nucleotide in sequence.subSequence(0..30)) {
            previousHash = BuildKmerIndex.updateKmerHashAndReverseCompliment(previousHash, nucleotide)
        }

        //start using kmers starting with the 32nd nucleotide
        //lookup hapids and add to the list
        for (nucleotide in sequence.subSequence(31 until sequence.length)) {
            previousHash = BuildKmerIndex.updateKmerHashAndReverseCompliment(previousHash, nucleotide)
            val minHash = min(previousHash.first, previousHash.second).toLong()
            val hapidsMatched = rangeHapidMapFromKmerHash(minHash, kmerHashOffsetMap, refrangeToBitSet, rangeToHapidIndexMap)
            for (entry in hapidsMatched) {
                val hapidList = rangeToHapidMap[entry.key]
                if (hapidList == null) rangeToHapidMap.put(entry.key, entry.value)
                else hapidList.addAll(entry.value)
            }
        }
    }

    //if no hapids map to this read, return an empty set
    if (rangeToHapidMap.size == 0) return setOf()

    //all hapids should be from the same reference range,
    //but if some are not then only those from the majority reference range should be used, so...
    val filteredHapidList = hapidsFromOneReferenceRange(rangeToHapidMap, minSameReferenceRange)

    //count the hapids
    val hapidCounts = filteredHapidList.groupingBy { it }.eachCount()

    //determine maxcount then create a hapidset from the hapids with maxcount

    val maxcount = hapidCounts.values.maxOrNull()
    return if (maxcount == null) setOf<Int>()
    else {
        val retainCount = ceil(maxcount.toDouble() * minProportionOfMaxCount)
        hapidCounts.entries.filter { (_, hapidCount) -> hapidCount >= retainCount }.map { it.key }.toHashSet()
    }
}

/**
 * For a kmer hash value, [kmerHash], generates a map of reference range id -> haplotype id list
 * for all of the haplotypes containing that kmer. The method returns an empty map if the kmer hash
 * does not map to any haplotypes.
 */
private fun rangeHapidMapFromKmerHash(kmerHash: Long,
                                      kmerHashOffsetMap: Long2LongOpenHashMap,
                                      refrangeToBitSet: Map<Int, BitSet>,
                                      rangeToHapidIndexMap: Map<Int, Map<Int, Int>>): Map<Int, MutableList<Int>> {
    //get the encoded refrangeId, offset from the kmerMap
    val encodedOffset = kmerHashOffsetMap[kmerHash]

    //a Long2LongOpenHashSet returns 0L rather than null when the kmerHash is not in the map.
    //so if encodedOffset = 0, the kmerhash does not map to any haplotypes. Return an empty map.
    if (encodedOffset == 0L) return mapOf()

    //decode the long encoded refrangeId, offset
    val rangeId = (encodedOffset shr 32).toInt()
    val offset = (encodedOffset and 0xFFFFFFFF).toInt()

    val hapidIndex = rangeToHapidIndexMap[rangeId]
    val hapidBitSet = refrangeToBitSet[rangeId]
    if (hapidIndex == null || hapidBitSet == null) return mapOf()

    val hapidList = hapidIndex.entries.filter { (_, index) -> hapidBitSet.get(offset + index) }
        .map { it.key }.toMutableList()
    return mapOf(rangeId to hapidList)
}

/**
 * Takes a map of range -> (haplotype id list) then determines whether at least [minSameReferenceRange] of them
 * map to the same reference range and returns only those hapids mapping to that reference range.
 * Returns an empty list if there is no reference range meeting that criterion.
 */
private fun hapidsFromOneReferenceRange(rangeHapidMap: Map<Int, List<Int>>, minSameReferenceRange: Double = 0.9): List<Int> {
    val numberOfHapids = rangeHapidMap.values.sumOf { it.size }
    val maxHapidCount = rangeHapidMap.values.maxOf { it.size }
    if (maxHapidCount.toDouble() / numberOfHapids.toDouble() < minSameReferenceRange) return listOf()
    val entryWithMaxCount = rangeHapidMap.entries.find { it.value.size == maxHapidCount }
    check(entryWithMaxCount != null) { "No entry has the max count.  This should never happen." }
    return entryWithMaxCount.value
}


/**
 * Function to export the read mapping files to disk.  These files can then be read in to be used in path finding
 */
fun exportReadMapping(outputFileName: String, hapIdMapping: Map<List<String>, Int>, sampleName : String, fastqFiles: Pair<String,String>) {
    getBufferedWriter(outputFileName).use { output ->
        output.write("#sampleName=${sampleName}\n")
        output.write("#filename1=${fastqFiles.first}\n")
        if(fastqFiles.second.isNotEmpty()) {
            output.write("#filename2=${fastqFiles.second}\n")
        }
        output.write("HapIds\tcount\n")
        hapIdMapping.keys.forEach { output.write("${it.joinToString(",")}\t${hapIdMapping[it]}\n") }
    }
}

/**
 * Function to read in the Readmapping file into a Map<List<String>,Int>
 */
fun importReadMapping(inputFileName: String) : Map<List<String>, Int> {
    return getBufferedReader(inputFileName).readLines()
        .filter { !it.startsWith("#") } //remove any of the metadata
        .filter { !it.startsWith("HapIds") }.associate {
            val lineSplit = it.split("\t")
            val hapIds = lineSplit[0].split(",")
                .map { hapId -> hapId }
                .toList()
            val count = lineSplit[1].toInt()

            Pair(hapIds, count)
        }
}