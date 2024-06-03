package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedReader
import htsjdk.samtools.fastq.FastqReader
import it.unimi.dsi.fastutil.longs.Long2LongOpenHashMap
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
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
import kotlin.math.max
import kotlin.math.min

/**
 * data class to hold the KmerMap information.  From the PHGv1 source code:
 */
data class KmerMapData(val rangeToBitSetMap: Map<ReferenceRange, BitSet>, val kmerHashToLongMap: Long2ObjectOpenHashMap<List<RefRangeOffset>>)

data class RefRangeOffset(val refRange: ReferenceRange, val offset: Long)

data class KeyFileData(val sampleName: String, val file1: String, val file2: String = "")

private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.pathing.AlignmentUtils")

class AlignmentUtils {
    companion object {
        /**
         * Function to align the reads coming from a single ended Fastq or a pair of fastq files to the haplotype graph kmer
         * index found in [kmerIndexFile] and create a readMapping file which will be written to the [outputDir].
         */
        fun alignReadsToHaplotypes(
            graph: HaplotypeGraph,
            kmerIndexFile: String,
            keyFileRecords: List<KeyFileData>,
            outputDir: String,
            numThreads: Int = 5,
            minProportionOfMaxCount: Double = 1.0,
            minSameReferenceRange: Double = 0.9,
        ) {
            val kmerIndexMap = loadKmerMaps(kmerIndexFile, graph)

            for (keyFileRecord in keyFileRecords) {
                val hapIdMapping = processReadMappingForKeyFileRecord(keyFileRecord, kmerIndexMap, graph, numThreads, minProportionOfMaxCount, minSameReferenceRange)

                //export the read mapping to disk
                //Use the first file name as the readmapping output name
                val inputFile1 = File(keyFileRecord.file1)
                //Todo does this need to work with .fastq.gz files?
                val outputFileName = "${outputDir}/${inputFile1.nameWithoutExtension}_readMapping.txt"


                exportReadMapping(
                    outputFileName,
                    hapIdMapping,
                    keyFileRecord.sampleName,
                    Pair(keyFileRecord.file1, keyFileRecord.file2)
                )
            }

        }

        /**
         * Loads a kmer hash map file to memory for use in mapping reads.
         */
        fun loadKmerMaps(filename: String, graph: HaplotypeGraph): KmerMapData {
            //Load the contents of the file into
            //rangeHapidMap: a map of refRange to an Array<String> of the haplotype ids in the ReferenceRange
            //    and the BitSet of all hapid sets.  This Bitset keeps track of what haplotypes are a part of a specific
            //    haplotype set in an efficient manner.  It is then used as a lookup for when a read hits the set for each
            //    haplotype to get attributed to having that read.
            //kmerHashmap: a map of kmer hash to reference range and offset into its BitSet, encoded as a long
            //These data structures all the reference range and haplotype set to be looked up for a kmer has


            val rangeToBitSetMap = mutableMapOf<ReferenceRange, BitSet>()
            val kmerToRefRangeAndOffsetMap = Long2ObjectOpenHashMap<MutableList<RefRangeOffset>>()
            var lineCount = 0
            var totalLineCount = 0
            var refrange = ReferenceRange("NULL", 0, 0)

            val refRangeToIdMap = graph.refRangeToIndexMap()
            getBufferedReader(filename).useLines {
                it.forEach { inputStr ->
                    totalLineCount++

                    lineCount = when (inputStr.first()) {
                        '>' -> 1
                        else -> lineCount + 1
                    }

                    when (lineCount) {
                        1 -> {
                            //line 1 is the range
                            refrange = ReferenceRange.parse(inputStr.substring(1))
                        }

                        2 -> {
                            try {
                                val parsedLine = inputStr.split(",")
                                val myBitset = BitSet.valueOf(parsedLine.map { it.toLong() }.toLongArray())
                                rangeToBitSetMap[refrange] = myBitset
                            } catch (e: Exception) {
                                myLogger.error("error at line $totalLineCount for input = $inputStr")
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
                                    myLogger.warn("improperly formatted datapair at line $totalLineCount")
                                }

                                val hash = hashOffset[0].toLong()
                                val offset = refrangeLong or hashOffset[1].toLong()


                                val offsetList = kmerToRefRangeAndOffsetMap.get(hash)
                                if(offsetList == null) {
                                    val newOffsetList = mutableListOf(RefRangeOffset(refrange, offset))
                                    kmerToRefRangeAndOffsetMap.put(hash, newOffsetList)
                                } else {
                                    offsetList.add(RefRangeOffset(refrange, offset))
                                }

                            }
                        }

                        else -> {
                            throw IllegalArgumentException("Line count = $lineCount in kmerMap file format at line $totalLineCount in file $filename")
                        }
                    }
                }
            }

            //convert kmerToRefRangeAndOffsetMap to be immutable
            val immutableKmerToRefRangeMap = Long2ObjectOpenHashMap<List<RefRangeOffset>>()
            for(key in kmerToRefRangeAndOffsetMap.keys) {
                immutableKmerToRefRangeMap.put(key, kmerToRefRangeAndOffsetMap.get(key)!!.toList())
            }

            return KmerMapData(rangeToBitSetMap, immutableKmerToRefRangeMap)
        }


        /**
         * Function to do the read mapping process for a single sample.
         * This works with both single and paired end read files.
         */
        fun processReadMappingForKeyFileRecord(
            keyFileRecord: KeyFileData,
            kmerIndexMap: KmerMapData,
            graph: HaplotypeGraph,
            numThreads: Int = 5,
            minProportionOfMaxCount: Double = 1.0,
            minSameReferenceRange: Double = 0.9,
        ): Map<List<String>, Int> {
            val fastqFiles = Pair(keyFileRecord.file1, keyFileRecord.file2)

            val hapidSetCount = mutableMapOf<List<String>, Int>()
            myLogger.info("reading records from the fastq file(s): ${fastqFiles.first}, ${fastqFiles.second}")
            var readCount = 0

            val rangeIdToBitsetMap =
                convertRefRangeToIdBitsetMap(kmerIndexMap.rangeToBitSetMap, graph.refRangeToIndexMap())

            FastqReader(bufferedReader(fastqFiles.first)).use { reader1 ->
                //Need to do this otherwise we cannot handle single and paired end in an efficient way
                val reader2 = if (fastqFiles.second.isNotEmpty()) {
                    FastqReader(bufferedReader(fastqFiles.second))
                } else {
                    null
                }

                //Setting up the channels for the coroutines
                val readStringChannel = Channel<Pair<String, String>>(100)
                val sortedHapidListChannel = Channel<List<String>>(100)
                val numberOfMappingThreads = max(1, numThreads - 2)

                //Run the blocking coroutine to read the fastq files and process the reads
                runBlocking {
                    launch(Dispatchers.IO) {
                        while (reader1.hasNext() && reader2?.hasNext() ?: true) {
                            val seq1 = reader1.next().readString
                            val seq2 = reader2?.next()?.readString ?: ""

                            readStringChannel.send(Pair(seq1, seq2))
                            readCount++

                        }
                        readStringChannel.close()
                    }

                    val processReadsJobList: MutableList<Job> = mutableListOf()
                    repeat(numberOfMappingThreads) {
                        processReadsJobList.add(launch(Dispatchers.Default) {
                            processReads(
                                readStringChannel,
                                sortedHapidListChannel,
                                minProportionOfMaxCount,
                                minSameReferenceRange,
                                kmerIndexMap.kmerHashToLongMap,
                                rangeIdToBitsetMap,
                                graph.refRangeIdToHapIdMap()
                            )
                        })
                    }

                    val addListsToMapJob = launch {
                        addListsToMap(hapidSetCount, sortedHapidListChannel)
                    }

                    processReadsJobList.joinAll()
                    sortedHapidListChannel.close()
                    addListsToMapJob.join()
                }
                //Need to close the reader2 if it exists
                reader2?.close()
            }

            return hapidSetCount
        }

        /**
         * Function to convert a map of ReferenceRange to BitSet to a map of RefRangeId to BitSet
         */
        fun convertRefRangeToIdBitsetMap(
            rangeToBitSetMap: Map<ReferenceRange, BitSet>,
            refRangeToIndexMap: Map<ReferenceRange, Int>
        ): Map<Int, BitSet> {
            val rangeIdToBitsetMap = mutableMapOf<Int, BitSet>()
            for ((range, bitset) in rangeToBitSetMap) {
                val rangeId = refRangeToIndexMap[range]!!
                rangeIdToBitsetMap[rangeId] = bitset
            }
            return rangeIdToBitsetMap
        }

        /**
         * Takes a [hapidSetCounts] map and a [hapLists] channel and adds the hapLists to the map.
         */
        suspend fun addListsToMap(
            hapidSetCounts: MutableMap<List<String>, Int>,
            hapLists: ReceiveChannel<List<String>>
        ) {
            for (hapList in hapLists) {
                val count = hapidSetCounts[hapList]
                if (count == null) hapidSetCounts[hapList] = 1
                else hapidSetCounts[hapList] = count + 1
            }
        }


        /**
         * Function to process a read or a readPair into their corresponding hapIdSet.  This will work with both single and paired end reads.
         * To Work with single ended reads, the second String in the Pair should be an empty String.
         */
        suspend fun processReads(
            reads: ReceiveChannel<Pair<String, String>>,
            sortedLists: SendChannel<List<String>>,
            minProportionOfMaxCount: Double = 1.0,
            minSameReferenceRange: Double = 0.9,
            kmerHashOffsetMap: Long2ObjectOpenHashMap<List<RefRangeOffset>>,
            refrangeToBitSet: Map<Int, BitSet>,
            rangeToHapidIndexMap: Map<Int, Map<String, Int>>
        ) {
            for (pairOfReads in reads) {
                val result1 = readToHapidSet(
                    pairOfReads.first,
                    minProportionOfMaxCount,
                    minSameReferenceRange,
                    kmerHashOffsetMap,
                    refrangeToBitSet,
                    rangeToHapidIndexMap
                )
                val result2 = if (pairOfReads.second.isNotEmpty()) {
                    readToHapidSet(
                        pairOfReads.second,
                        minProportionOfMaxCount,
                        minSameReferenceRange,
                        kmerHashOffsetMap,
                        refrangeToBitSet,
                        rangeToHapidIndexMap
                    )
                } else {
                    setOf<String>()
                }

                val intersectResult = if (result2.isNotEmpty()) {
                    result1.intersect(result2)
                } else {
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
        fun readToHapidSet(
            read: String,
            minProportionOfMaxCount: Double = 1.0,
            minSameReferenceRange: Double = 0.9,
            kmerHashOffsetMap: Long2ObjectOpenHashMap<List<RefRangeOffset>>,
            refrangeToBitSet: Map<Int, BitSet>,
            rangeToHapidIndexMap: Map<Int, Map<String, Int>>
        ): Set<String> {
            //generate kmer hash from the read
            val rangeToHapidMap = mutableMapOf<Int, MutableList<String>>()
            val splitList = read
                .split("[^ACGT]+".toRegex())
                .filter { it.length > 31 }
            for (sequence in splitList) {
                extractKmersFromSequence(
                    sequence,
                    kmerHashOffsetMap,
                    refrangeToBitSet,
                    rangeToHapidIndexMap,
                    rangeToHapidMap
                )
            }

            //if no hapids map to this read, return an empty set
            if (rangeToHapidMap.size == 0) return setOf()

            //all hapids should be from the same reference range,
            //but if some are not then only those from the majority reference range should be used, so...
            //TODO Make this work with multiple reference ranges
            val filteredHapidList = hapidsFromOneReferenceRange(rangeToHapidMap, minSameReferenceRange)

            //count the hapids
            val hapidCounts = filteredHapidList.groupingBy { it }.eachCount()

            //determine maxcount then create a hapidset from the hapids with maxcount

            val maxcount = hapidCounts.values.maxOrNull()
            return if (maxcount == null) setOf<String>()
            else {
                val retainCount = ceil(maxcount.toDouble() * minProportionOfMaxCount)
                hapidCounts.entries.filter { (_, hapidCount) -> hapidCount >= retainCount }.map { it.key }.toHashSet()
            }
        }

        /**
         * Function to extract kmers from a sequence and add them to the [rangeToHapidMap] for the reference range
         * This will work by doing a single pass over the read sequence and building the kmer hash and then looking up what hapIds belong to that kmer.
         */
        fun extractKmersFromSequence(
            sequence: String,
            kmerHashOffsetMap: Long2ObjectOpenHashMap<List<RefRangeOffset>>,
            refrangeToBitSet: Map<Int, BitSet>,
            rangeToHapidIndexMap: Map<Int, Map<String, Int>>,
            rangeToHapidMap: MutableMap<Int, MutableList<String>>
        ) {
            var previousHash = Pair(0L, 0L)

            //for first 31 nucleotides just update the hash
            for (nucleotide in sequence.subSequence(0..30)) {
                previousHash = BuildKmerIndex.updateKmerHashAndReverseCompliment(previousHash, nucleotide)
            }

            //start using kmers starting with the 32nd nucleotide
            //lookup hapids and add to the list
            for (nucleotide in sequence.subSequence(31 until sequence.length)) {
                previousHash = BuildKmerIndex.updateKmerHashAndReverseCompliment(previousHash, nucleotide)
                val minHash = min(previousHash.first, previousHash.second)
                val hapidsMatched =
                    rangeHapidMapFromKmerHash(minHash, kmerHashOffsetMap, refrangeToBitSet, rangeToHapidIndexMap)
                for (entry in hapidsMatched) {
                    val hapidList = rangeToHapidMap[entry.key]
                    if (hapidList == null) rangeToHapidMap.put(entry.key, entry.value)
                    else hapidList.addAll(entry.value)
                }
            }
        }

        /**
         * For a kmer hash value, [kmerHash], generates a map of reference range id -> haplotype id list
         * for all of the haplotypes containing that kmer. The method returns an empty map if the kmer hash
         * does not map to any haplotypes.
         */
        fun rangeHapidMapFromKmerHash(
            kmerHash: Long,
            kmerHashOffsetMap: Long2ObjectOpenHashMap<List<RefRangeOffset>>,
            refrangeToBitSet: Map<Int, BitSet>,
            rangeToHapidIndexMap: Map<Int, Map<String, Int>>
        ): Map<Int, MutableList<String>> {

            val encodedOffsets = kmerHashOffsetMap.getOrDefault(kmerHash, mutableListOf())     //[kmerHash]

            //Loop through the offsets as we can have multiple offsets for a single kmerHash across multiple reference ranges
            val rangeToHapidMap = mutableMapOf<Int, MutableList<String>>()
            for (encodedOffset in encodedOffsets) {
                val (rangeId, offset) = decodeRangeIdAndOffset(encodedOffset.offset)
                val hapidIndex = rangeToHapidIndexMap[rangeId]
                val hapidBitSet = refrangeToBitSet[rangeId]
                if (hapidIndex == null || hapidBitSet == null) continue

                val hapidList = hapidIndex.entries.filter { (_, index) -> hapidBitSet.get(offset + index) }
                    .map { it.key }.toMutableList()
                rangeToHapidMap[rangeId] = hapidList
            }

            return rangeToHapidMap
        }

        /**
         * Takes a long encoded refrangeId, offset into a Long and returns the rangeId and offset as a Pair<Int,Int>
         */
        fun decodeRangeIdAndOffset(encodedOffset: Long): Pair<Int, Int> {
            //decode the long encoded refrangeId, offset
            val rangeId = (encodedOffset shr 32).toInt()
            val offset = (encodedOffset and 0xFFFFFFFF).toInt()
            return Pair(rangeId, offset)
        }

        /**
         * Takes a map of range -> (haplotype id list) then determines whether at least [minSameReferenceRange] of them
         * map to the same reference range and returns only those hapids mapping to that reference range.
         * Returns an empty list if there is no reference range meeting that criterion.
         */
        fun hapidsFromOneReferenceRange(
            rangeHapidMap: Map<Int, List<String>>,
            minSameReferenceRange: Double = 0.9
        ): List<String> {
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
        fun exportReadMapping(
            outputFileName: String,
            hapIdMapping: Map<List<String>, Int>,
            sampleName: String,
            fastqFiles: Pair<String, String>
        ) {
            getBufferedWriter(outputFileName).use { output ->
                output.write("#sampleName=${sampleName}\n")
                output.write("#filename1=${fastqFiles.first}\n")
                if (fastqFiles.second.isNotEmpty()) {
                    output.write("#filename2=${fastqFiles.second}\n")
                }
                output.write("HapIds\tcount\n")
                hapIdMapping.keys.forEach { output.write("${it.joinToString(",")}\t${hapIdMapping[it]}\n") }
            }
        }

        /**
         * Function to read in the Readmapping file into a Map<List<String>,Int>
         */
        fun importReadMapping(inputFileName: String): Map<List<String>, Int> {
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

        /**
         * Function to merge multiple read mapping lists into a single list by summing the counts
         * for each read mapping set in the lists.
         */
        fun mergeReadMappings(listOfReadMappings: List<Map<List<String>, Int>>): Map<List<String>, Int> {
            if (listOfReadMappings.isEmpty()) return mapOf()
            if (listOfReadMappings.size == 1) return listOfReadMappings[0]
            val mergedMappings = listOfReadMappings[0].toMutableMap()
            listOfReadMappings.drop(1).forEach {
                it.entries.forEach { (haplist, count) ->
                    val newCount = count + (mergedMappings[haplist] ?: 0)
                    mergedMappings.put(haplist, newCount)
                }
            }
            return mergedMappings
        }
    }
}