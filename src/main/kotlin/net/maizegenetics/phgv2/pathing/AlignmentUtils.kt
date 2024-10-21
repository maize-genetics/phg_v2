package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedReader
import htsjdk.samtools.fastq.FastqReader
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

data class ReadIdAndSeq(val readId: String, val seq1: String, val seq2: String = "")

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
            runDiagnostic: Boolean = false
        ) {
            myLogger.info("Loading in Kmer Index Map")
            val startTime = System.nanoTime()
            val kmerIndexMap = loadKmerMaps(kmerIndexFile, graph)
            myLogger.info("Kmer Index Map loaded in ${(System.nanoTime() - startTime) / 1e9} seconds")
            myLogger.info("Starting read mapping process")
            for (keyFileRecord in keyFileRecords) {

                //Build the output file names
                //Use the first file name as the readmapping output name
                val inputFile1 = File(keyFileRecord.file1)

                val inputFile1NoExt = if(inputFile1.name.endsWith(".gz")) {
                    inputFile1.nameWithoutExtension.substringBeforeLast(".")
                } else {
                    inputFile1.nameWithoutExtension
                }

                val outputFileName = "${outputDir}/${inputFile1NoExt}_readMapping.txt"

                val diagnosticFileName = if(runDiagnostic) {
                    "${outputDir}/${inputFile1.nameWithoutExtension}_diagnostic.txt"
                } else {
                    ""
                }

                myLogger.info("Processing reads for sample ${keyFileRecord.sampleName}")
                val hapIdMapping = processReadMappingForKeyFileRecord(
                    keyFileRecord,
                    kmerIndexMap,
                    graph,
                    numThreads,
                    minProportionOfMaxCount,
                    minSameReferenceRange,
                    diagnosticFileName
                )
                //export the read mapping to disk
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

            val readFileStartTime = System.nanoTime()
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

                                if(kmerToRefRangeAndOffsetMap.containsKey(hash)) {
                                    kmerToRefRangeAndOffsetMap[hash]!!.add(RefRangeOffset(refrange, offset))
                                } else {
                                    kmerToRefRangeAndOffsetMap[hash] = mutableListOf(RefRangeOffset(refrange, offset))
                                }
                            }
                        }

                        else -> {
                            throw IllegalArgumentException("Line count = $lineCount in kmerMap file format at line $totalLineCount in file $filename")
                        }
                    }
                }
            }
            myLogger.info("KmerMap read in ${(System.nanoTime() - readFileStartTime) / 1e9} seconds")

            return KmerMapData(rangeToBitSetMap, kmerToRefRangeAndOffsetMap as Long2ObjectOpenHashMap<List<RefRangeOffset>>)
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
            diagnosticFileName: String = ""
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
                val readStringChannel = Channel<ReadIdAndSeq>(100)
                val sortedHapidListChannel = Channel<List<String>>(100)
                val diagnosticChannel = if(diagnosticFileName.isNotEmpty()) Channel<String>(100) //This channel will be used to handle output of diagnostic information
                                        else null
                val numberOfMappingThreads = max(1, numThreads - 2)

                //Run the blocking coroutine to read the fastq files and process the reads
                runBlocking {
                    launch(Dispatchers.IO) {
                        while (reader1.hasNext() && reader2?.hasNext() ?: true) {
                            val reader1CurrentRecord = reader1.next()
                            val seq1 = reader1CurrentRecord.readString
                            val seq2 = reader2?.next()?.readString ?: ""

                            readStringChannel.send(ReadIdAndSeq(reader1CurrentRecord.readName,seq1, seq2))
                            readCount++

                        }
                        readStringChannel.close()
                    }

                    val processReadsJobList: MutableList<Job> = mutableListOf()
                    repeat(numberOfMappingThreads) {
                        processReadsJobList.add(launch(Dispatchers.Default) {
                            processReadsMultipleRefRanges(
                                readStringChannel,
                                sortedHapidListChannel,
                                diagnosticChannel,
                                kmerIndexMap.kmerHashToLongMap,
                                rangeIdToBitsetMap,
                                graph.refRangeIdToHapIdMap(),
                                minProportionOfMaxCount,
                                minSameReferenceRange
                            )
                        })
                    }

                    val addListsToMapJob = launch {
                        addListsToMap(hapidSetCount, sortedHapidListChannel)
                    }

                    val diagnosticJob = if(diagnosticFileName.isNotEmpty()) {
                        //call diagnostic file writer
                        launch {
                            writeDiagnosticFile(diagnosticFileName, diagnosticChannel)
                        }
                    } else {
                        null
                    }

                    processReadsJobList.joinAll()
                    sortedHapidListChannel.close()
                    addListsToMapJob.join()
                    diagnosticChannel?.close()
                    diagnosticJob?.join()
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
         * Function to handle diagnostic output lines so we don't have any multithreading writing issues as only 1 thread has the writer
         */
        suspend fun writeDiagnosticFile(diagnosticFileName: String, diagnosticChannel:ReceiveChannel<String>?) {
            if(diagnosticChannel == null) return
            getBufferedWriter(diagnosticFileName).use { writer ->
                writer.write("ReadName\tMapped\tReferenceRange\tHapIds\tReason\n")
                for (line in diagnosticChannel) {
                    writer.write(line)
                    writer.newLine()
                }
            }
        }

        /**
         * Function to process the reads for multiple ReferenceRanges.  This will work with both paired and single ended reads.
         * To Work with single ended reads, the second String in the Pair should be an empty String.
         */
        private suspend fun processReadsMultipleRefRanges(
            reads: ReceiveChannel<ReadIdAndSeq>,
            sortedLists: SendChannel<List<String>>,
            diagnosticChannel: SendChannel<String>?,
            kmerHashOffsetMap: Long2ObjectOpenHashMap<List<RefRangeOffset>>,
            refrangeToBitSet: Map<Int, BitSet>,
            rangeToHapidIndexMap: Map<Int, Map<String, Int>>,
            minProportionOfMaxCount: Double = 1.0,
            minSameReferenceRange: Double = 0.9
        ) {
            for (pairOfReads in reads) {
                val result = readToHapIdSetForPairedReads(
                    pairOfReads.readId,
                    Pair(pairOfReads.seq1, pairOfReads.seq2),
                    kmerHashOffsetMap,
                    refrangeToBitSet,
                    rangeToHapidIndexMap,
                    minProportionOfMaxCount,
                    minSameReferenceRange,
                    diagnosticChannel
                )

                //Process the mappings for each reference range
                for(entry in result) {
                    val sortedIds = entry.value.sorted()
                    sortedLists.send(sortedIds)
                    diagnosticChannel?.send("${pairOfReads.readId}\ttrue\t${entry.key}\t${sortedIds.joinToString(",")}\tMapped")
                }
            }
        }

        /**
         * Function to process the hapIdHits for a read pair.  This will work with both paired and single ended reads.
         * If the second read is empty, then it will be treated as a single ended read and will return only 1 set of hapids
         */
        private suspend fun processHapIdHits(
            readName: String,
            result1: Map<Int, Set<String>>,
            result2: Map<Int, Set<String>>,
            sortedLists: SendChannel<List<String>>,
            diagnosticChannel: SendChannel<String>?,
        ) {
            for (entry in result1) {
                val rangeId = entry.key
                val hapIds = entry.value

                val intersectResult = if (result2.containsKey(rangeId)) {
                    hapIds.intersect(result2[rangeId]!!)
                } else {
                    hapIds
                }
                if (intersectResult.isNotEmpty()) {
                    val result = intersectResult.sorted()
                    sortedLists.send(result)
                    diagnosticChannel?.send("$readName\ttrue\t$rangeId\t${result.joinToString(",")}\tMapped")
                }
                else diagnosticChannel?.send("$readName\tfalse\tNA\tNA\tNoSharedHapsInPair")
            }
        }

        /**
         * Function to extract the hapIdHits for a read pair.  This will work with both paired and single ended reads.
         */
        private suspend fun extractHapIdHitsForReadPair(
            pairOfReads: ReadIdAndSeq,
            kmerHashOffsetMap: Long2ObjectOpenHashMap<List<RefRangeOffset>>,
            refrangeToBitSet: Map<Int, BitSet>,
            rangeToHapidIndexMap: Map<Int, Map<String, Int>>,
            minProportionOfMaxCount: Double = 1.0,
            limitSingleRefRange: Boolean = true,
            minSameReferenceRange: Double = 0.9,
            diagnosticChannel: SendChannel<String>?
        ): Pair<Map<Int, Set<String>>, Map<Int, Set<String>>> {
            println("File1")
            val result1 = readToHapIdSetMultipleRefRanges(
                pairOfReads.readId,
                pairOfReads.seq1,
                kmerHashOffsetMap,
                refrangeToBitSet,
                rangeToHapidIndexMap,
                minProportionOfMaxCount,
                limitSingleRefRange,
                minSameReferenceRange,
                diagnosticChannel
            )
            val result2 = if (pairOfReads.seq2.isNotEmpty()) {
                println("File2")
                readToHapIdSetMultipleRefRanges(
                    pairOfReads.readId,
                    pairOfReads.seq2,
                    kmerHashOffsetMap,
                    refrangeToBitSet,
                    rangeToHapidIndexMap,
                    minProportionOfMaxCount,
                    limitSingleRefRange,
                    minSameReferenceRange,
                    diagnosticChannel
                )
            } else {
                emptyMap()
            }
            return Pair(result1, result2)
        }

        /**
         * Takes a [read] and generates a list of hapids to which its kmers map. Returns hapids from multiple
         * reference ranges. If no kmers map or if fewer than [minSameReferenceRange] of the kmers map to a single
         * reference range an empty Set will be returned.
         */
        suspend fun readToHapIdSetMultipleRefRanges(
            readId: String,
            read: String,
            kmerHashOffsetMap: Long2ObjectOpenHashMap<List<RefRangeOffset>>,
            refrangeToBitSet: Map<Int, BitSet>,
            rangeToHapidIndexMap: Map<Int, Map<String, Int>>,
            minProportionOfMaxCount: Double = 1.0,
            limitSingleRefRange: Boolean = true,
            minSameReferenceRange: Double = 0.9,
            diagnosticChannel: SendChannel<String>?,

        ): Map<Int, Set<String>> {
            //generate kmer hash from the read
            var rangeToHapIdMap = mutableMapOf<Int, MutableList<String>>()
            extractKmersForSingleRead(
                read,
                kmerHashOffsetMap,
                refrangeToBitSet,
                rangeToHapidIndexMap,
                rangeToHapIdMap
            )

            println("RangeToHapIdMap Outer: $rangeToHapIdMap IsEmpty: ${rangeToHapIdMap.isEmpty()}")
            //if no hapids map to this read, return an empty set
            if (rangeToHapIdMap.isEmpty()) {
                println("RangeToHapIdMap If Statement: $rangeToHapIdMap IsEmpty: ${rangeToHapIdMap.isEmpty()}")
                diagnosticChannel?.send("$readId\tfalse\tNA\tNA\tNoHapIdsMapped")
                return emptyMap()
            }

            //If the user requests singleRefRange Mode we filter down using the [minSameReferenceRange] parameter
            if(limitSingleRefRange) {
                //all hapids should be from the same reference range,
                //but if some are not then only those from the majority reference range should be used, so...
                rangeToHapIdMap = filterHapIdsToOneReferenceRange(rangeToHapIdMap, minSameReferenceRange)
            }

            if(rangeToHapIdMap.isEmpty()) {
                diagnosticChannel?.send("$readId\tfalse\tNA\tNA\tKmersFilteredByRefRange")
            }

            //hapIds are already grouped by refRange so we just need to process each ranges hapIds and make sure that we
            // have good enough coverage
            val rangeIdToHapIdSet = rangeToHapIdMap.map { (rangeId, hapIds) ->
               Pair(rangeId, filterHapIdsByKmerCount(hapIds, minProportionOfMaxCount))
            }.toMap()

            if(rangeIdToHapIdSet.isEmpty()) {
                diagnosticChannel?.send("$readId\tfalse\tNA\tNA\tHapIdsFilteredByKmerCount")
            }

            return rangeIdToHapIdSet
        }

        /**
         * Takes a [read] and generates a list of hapids to which its kmers map. Returns hapids from multiple
         * reference ranges. If no kmers map or if fewer than [minSameReferenceRange] of the kmers map to a single
         * reference range an empty Set will be returned.
         */
        suspend fun readToHapIdSetForPairedReads(
            readId: String,
            reads: Pair<String,String?>,
            kmerHashOffsetMap: Long2ObjectOpenHashMap<List<RefRangeOffset>>,
            refrangeToBitSet: Map<Int, BitSet>,
            rangeToHapidIndexMap: Map<Int, Map<String, Int>>,
            minProportionOfMaxCount: Double = 1.0,
            minSameReferenceRange: Double = 0.9,
            diagnosticChannel: SendChannel<String>?
        ): Map<Int, Set<String>> {
            //generate kmer hash from the read
            var rangeToHapIdMap = mutableMapOf<Int, MutableList<String>>()

            //Extract out the kmer hashes from the pair of reads
            extractKmersForPairedReads(reads, kmerHashOffsetMap, refrangeToBitSet, rangeToHapidIndexMap, rangeToHapIdMap)

            //if no hapids map to this read, return an empty set
            if (rangeToHapIdMap.isEmpty()) {
                diagnosticChannel?.send("$readId\tfalse\tNA\tNA\tNoHapIdsMapped")
                return emptyMap()
            }

            //all hapids should be from the same reference range,
            //but if some are not then only those from the majority reference range should be used, so...
            rangeToHapIdMap = filterHapIdsToOneReferenceRange(rangeToHapIdMap, minSameReferenceRange)

            if(rangeToHapIdMap.isEmpty()) {
                diagnosticChannel?.send("$readId\tfalse\tNA\tNA\tKmersFilteredByRefRange")
                return emptyMap()
            }

            //hapIds are already grouped by refRange so we just need to process each ranges hapIds and make sure that we
            // have good enough coverage
            val rangeIdToHapIdSet = rangeToHapIdMap.map { (rangeId, hapIds) ->
                Pair(rangeId, filterHapIdsByKmerCount(hapIds, minProportionOfMaxCount))
            }.toMap()

            if(rangeIdToHapIdSet.isEmpty()) {
                diagnosticChannel?.send("$readId\tfalse\tNA\tNA\tHapIdsFilteredByKmerCount")
            }

            return rangeIdToHapIdSet
        }

        /**
         * Function to extract out the kmers for a single or Pair of reads
         * If reads.second is null we skip it.
         * This allows for simple processing of the read pair
         */
        private fun extractKmersForPairedReads(
            reads: Pair<String, String?>,
            kmerHashOffsetMap: Long2ObjectOpenHashMap<List<RefRangeOffset>>,
            refrangeToBitSet: Map<Int, BitSet>,
            rangeToHapidIndexMap: Map<Int, Map<String, Int>>,
            rangeToHapIdMap: MutableMap<Int, MutableList<String>>
        ) {
            val (read1, read2) = reads
            extractKmersForSingleRead(read1, kmerHashOffsetMap, refrangeToBitSet, rangeToHapidIndexMap, rangeToHapIdMap)
            if(read2!=null) {
                extractKmersForSingleRead(read2, kmerHashOffsetMap, refrangeToBitSet, rangeToHapidIndexMap, rangeToHapIdMap)
            }
        }

        /**
         * Extract out the kmers for a single read.
         * It first splits the read sequence by Ns and filters out any subsequence smaller than 32 bps
         * Then it extracts the kmers and minHashes.
         */
        private fun extractKmersForSingleRead(
            read1: String,
            kmerHashOffsetMap: Long2ObjectOpenHashMap<List<RefRangeOffset>>,
            refrangeToBitSet: Map<Int, BitSet>,
            rangeToHapidIndexMap: Map<Int, Map<String, Int>>,
            rangeToHapIdMap: MutableMap<Int, MutableList<String>>
        ) {
            val splitList = read1
                .split("[^ACGT]+".toRegex())
                .filter { it.length > 31 }
            for (sequence in splitList) {
                extractKmersFromSequence(
                    sequence,
                    kmerHashOffsetMap,
                    refrangeToBitSet,
                    rangeToHapidIndexMap,
                    rangeToHapIdMap
                )
            }
        }


        /**
         * Function to filter hapIds by the kmer count.  This will return a set of hapIds that have a count greater than
         * or equal to the [minProportionOfMaxCount] of the max count.
         */
        fun filterHapIdsByKmerCount(hapIds: List<String>, minProportionOfMaxCount: Double): Set<String> {
            val hapidCounts = hapIds.groupingBy { it }.eachCount()
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
                val minHash = min(previousHash.first.toULong(), previousHash.second.toULong()).toLong()
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
         * Function that builds a set of minHashes from a sequence.
         * This will work by doing a single pass over the read sequence and building the kmer hash set.
         */
        fun extractMinHash(sequence: String) : Set<Long> {
            var previousHash = Pair(0L, 0L)

            val minHashSet = mutableSetOf<Long>()

            //for first 31 nucleotides just update the hash
            for (nucleotide in sequence.subSequence(0..30)) {
                previousHash = BuildKmerIndex.updateKmerHashAndReverseCompliment(previousHash, nucleotide)
            }

            for (nucleotide in sequence.subSequence(31 until sequence.length)) {
                previousHash = BuildKmerIndex.updateKmerHashAndReverseCompliment(previousHash, nucleotide)
                val minHash = min(previousHash.first, previousHash.second)
                minHashSet.add(minHash)
            }

            return minHashSet
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
                if (hapidIndex == null || hapidBitSet == null) {
                    continue
                }

                val hapidList = hapidIndex.entries.filter { (entry, index) ->
                    hapidBitSet.get(offset + index)
                }
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

        fun filterHapIdsToOneReferenceRange(
            rangeHapIdMap: Map<Int,List<String>>,
            minSameReferenceRange: Double = 0.9
        ): MutableMap<Int,MutableList<String>> {
            val numberOfHapIds = rangeHapIdMap.values.sumOf { it.size }
            val maxHapidCount = rangeHapIdMap.values.maxOf { it.size }

            if(maxHapidCount.toDouble() / numberOfHapIds.toDouble() < minSameReferenceRange) return mutableMapOf()

            rangeHapIdMap.entries.find { it.value.size == maxHapidCount }?.let { entry ->
                return mutableMapOf(entry.key to entry.value.toMutableList())
            }
            return mutableMapOf()
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