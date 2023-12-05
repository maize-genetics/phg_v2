package net.maizegenetics.phgv2.pathing

import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.long
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap
import it.unimi.dsi.fastutil.longs.Long2LongOpenHashMap
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.utils.retrieveAgcContigs
import org.apache.logging.log4j.LogManager
import java.io.File
import java.lang.IllegalStateException
import java.util.*


import kotlin.math.ceil
import kotlin.math.min
import kotlin.time.DurationUnit
import kotlin.time.measureTimedValue

/**
 * Creates a Map of 32-mer hash -> hapid list for the haplotypes in a HaplotypeGraph. Only hashes observed
 * in exactly one reference range will be kept. Also, hashes will be retained only if they map to at
 * most [maxHaplotypeProportion] * the number of haplotyes in a reference range. Only hashes that pass the filter
 * ((hashValue and [hashMask]) == [hashFilterValue]) will be considered. For example, setting [hashMask] = 3u and [hashFilterValue] = 1u
 * only uses hashes from kmers ending in C. To filter on the final two positions set [hashMask] = 16u (0b1111).
 * [hashFilterValue] is based on the two bit encoding of nucleotides: A -> 0, C -> 1, G -> 2, T -> 3.
 * For example, the [hashFilterValue] for CG is 6u (0b0110).
 *
 * To use this class, ...
 */
class BuildKmerIndex: CliktCommand(help="Create a kmer index for a HaplotypeGraph") {

    private val myLogger = LogManager.getLogger(BuildKmerIndex::class.java)

    val tiledbPath by option(help = "Tile DB URI")
        .default("")

    val agcPath by option(help = "AGC fasta DB URI")
        .required()


    val maxHaplotypeProportion by option("-p", "--maxHapProportion", help = "only kmers mapping to less than or " +
            "equal to maxHapProportion of haplotypes in a reference range will be retained.")
        .double()
        .default(0.75)

    val hashMask by option("-m", "--hashMask", help = "with hashFilter, used to mask kmers for filtering. " +
            "Default uses only the last kmer nucleotide. Only change this if you know what you are doing.")
        .long()
        .default(3)

    val hashFilterValue by option("f", "--hashFilter", help = "Only hashes that pass the filter" +
            " ((hashValue and hashMask) == hashFilter) will be considered")
        .long()
        .default(1)

    val hvcfDir by option("--hvcf-dir", help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        .default("")


    override fun run() {

        TODO("Not yet implemented")
    }

    private fun buildHaplotypeGraph(): HaplotypeGraph {
        require(hvcfDir.isNotBlank() or tiledbPath.isNotBlank()) {"Either of --tiledb-path or --agc-path must be provided."}
        val timedValue = measureTimedValue {
            if(hvcfDir != "") {
                val pathList = File(hvcfDir).listFiles { file -> file.extension == "vcf" || file.name.endsWith("vcf.gz") }.map { it.path }
                HaplotypeGraph(pathList)
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
     * Finds the set of kmers meeting the conditions set by [maxHaplotypeProportion], the maximum proportion of
     * haplotypes in a reference range, and [hashMask] and [hashFilterValue], which determine which nucleotide(s)
     * are required to be in the final postion(s) of a kmer.
     *
     * This returns a HashMap of hash -> hapid list for all the kmers in the keep set.
     * Which allows the export to not need to do a second pass over the sequences to get the set of hapIds which contain the unique kmers.
     */
    fun processGraphKmers() : Long2ObjectOpenHashMap<Set<String>> {
        val keepMap = Long2ObjectOpenHashMap<Set<String>>()
        val discardSet = LongOpenHashSet()
        var rangeCount = 0
        val graph = buildHaplotypeGraph()
        val startTime = System.nanoTime()

        for (refrange in graph.ranges()) {
            if (rangeCount++ % 5000 == 0) {
                myLogger.debug("Processing range $rangeCount, keep set size = ${keepMap.size}, discard set size = ${discardSet.size}, elapse time ${(System.nanoTime() - startTime)/1e9} sec")
            }

            val hapidToSampleMap = graph.hapIdToSamples(refrange)
            val maxHaplotypes = ceil(hapidToSampleMap.size * maxHaplotypeProportion)

            //create a map of hash -> count of occurrences for all the haplotypes in this reference range
            //retrieveAgcContigs needs a list of contig@genome:start-end, for all of the regions in the haplotype alt headers
            //subtract 1 from altHeader positions because AGC positions are 0-based.
            //also need a map of source -> hapid in order to associate sequence from agc to the hapid that it came from
            val agcRangeList = mutableListOf<String>()
            val sourceHapidMap = mutableMapOf<String, String>()
            for ((hapid, alt) in graph.altHeaderMap.entries) {
                //mapping source name to hapid assumes there is only one hapid per source in a reference range
                //this seems safe, but it is being checked here just in case
                check(sourceHapidMap.contains(alt.source)) {"Two hapids from ${alt.source} at ${alt.regions[0].first.position}"}
                sourceHapidMap[alt.source] = hapid
                for (range in alt.regions) {
                    if (range.first.position <= range.second.position) {
                        agcRangeList.add("${range.first.contig}@${alt.source}:${range.first.position - 1}-${range.second.position - 1}")
                    } else {
                        agcRangeList.add("${range.first.contig}@${alt.source}:${range.second.position - 1}-${range.first.position - 1}")
                    }
                }
            }

            val agcdataMap = retrieveAgcContigs(agcPath, agcRangeList)
            //TODO method to get source name from the agcdataMap keys
            val hapidSeqMap = agcdataMap.entries.map { (agcid, seq) -> }

            val (kmerHashCounts, longToHapIdMap) = countKmerHashesForHaplotypeSequence(hapidSeqMap)

            for (hashCount in kmerHashCounts.entries) {
                val hashValue = hashCount.key

                //if the hash is in the discard set skip it
                if (discardSet.contains(hashValue)) continue

                when {
                    //if hash count >= numberOfHaplotype add it to the discard set
                    hashCount.value >= maxHaplotypes -> discardSet.add(hashValue)
                    //if the hash is already in the keepSet, it has been seen in a previous reference range
                    // so, remove it from the keep set and add it to the discard set
                    keepMap.containsKey(hashValue) -> {
                        keepMap.remove(hashValue)
                        discardSet.add(hashValue)
                    }
                    else -> {
                        keepMap[hashValue] = longToHapIdMap[hashValue]
                    }
                }
            }

        }
        myLogger.debug("Finished building kmer keep set, keep set size = ${keepMap.size}, discard set size = ${discardSet.size}, elapse time ${(System.nanoTime() - startTime)/1e9} sec")
        return keepMap
    }


    private fun getSequenceForReferenceRange() {

    }

    /**
     * Function to count the kmerHashes for a single reference range's haplotype nodes.
     * This gets put in 2 sets of maps.
     * One for the hash counts and one for a hash to a list of hapIds which contain that hash.
     * The sequenceList input is a map of hapidId -> list of sequences
     */
    private fun countKmerHashesForHaplotypeSequence(sequenceMap: Map<String, List<String>>) : Pair<Map<Long,Int>,Map<Long,Set<String>>> {
        //start by splitting sequence into subsequences without N's
        val mapOfHashes = Long2IntOpenHashMap()
        val mapOfHapIds = Long2ObjectOpenHashMap<MutableSet<String>>()

        for ((hapid, seqList) in sequenceMap) {
            for (sequence in seqList) {
                val splitList = sequence
                    .split("N+".toRegex())
                    .filter{it.length > 31}
                for (sequence in splitList) {
                    var previousHash = Pair(0L, 0L)

                    //for first 31 nucleotides just update the hash
                    for (nucleotide in sequence.subSequence(0..30)) {
                        previousHash = updateKmerHashAndReverseCompliment(previousHash, nucleotide)
                    }

                    //start using kmers starting with the 32nd nucleotide
                    for (nucleotide in sequence.subSequence(31 until sequence.length)) {
                        previousHash = updateKmerHashAndReverseCompliment(previousHash, nucleotide)
                        val minHash = min(previousHash.first, previousHash.second)
                        // Use only kmers with a specific ending nucleotide(s)
                        // hashMask determines how many positions will be used
                        // hashFilterValue determines which nucleotide will be kept at that position
                        //at this point the hash value, which is a ULong, must be converted to a Long in order to make
                        //use of the fastutils Long2IntOpenHashMap
                        if ((minHash and hashMask) == hashFilterValue) {
                            mapOfHashes.addTo(minHash,1)
                            if(mapOfHapIds.containsKey(minHash)) {
                                mapOfHapIds[minHash].add(hapid)
                            }
                            else {
                                mapOfHapIds[minHash] = mutableSetOf(hapid)
                            }
                        }
                    }
                }
            }
            //split sequence on N's then filter on length > 31

        }
        return Pair(mapOfHashes,mapOfHapIds)
    }

    /**
     * Function to save the kmer hash map to the specified file.  If writeToDB is set to true it will also update the haplotype list ID.
     */
    fun saveKmerHashesAndHapids(filePath: String, kmerMapToHapIds: Long2ObjectOpenHashMap<Set<Int>>, writeToDb: Boolean = true) {
        //TODO replace next two lines after graph (interface?) is finished
        val refRangeToHapidMap = getRefRangeToHapidMap(graph)
        val hapIdToRefRangeMap = graph.referenceRanges().flatMap { graph.nodes(it) }.map { Pair(it.id(),it.referenceRange().id()) }.toMap()

        //Walk through the map and associate kmerLongs with specific referenceRanges
        val refRangeToKmerSetMap = mutableMapOf<Int,MutableSet<Long>>()
        for(kmerMap in kmerMapToHapIds.long2ObjectEntrySet()) {
            val kmer = kmerMap.longKey
            val currentRefRange = hapIdToRefRangeMap[kmerMap.value.first()]!!
            if(refRangeToKmerSetMap.containsKey(currentRefRange)) {
                refRangeToKmerSetMap[currentRefRange]?.add(kmer)
            }
            else {
                refRangeToKmerSetMap[currentRefRange] = mutableSetOf(kmer)
            }
        }

//TODO replace haplotypeListId with something that allows the graph to be stored
        val haplotypeListId = if(writeToDb) {
            PHGdbAccess(DBLoadingUtils.connection(false)).use { phgdb ->
                getHaplotypeListIdForGraph(graph, phgdb)
            }
        }
        else {
            -1
        }


        var rangeCount = 0
        val startTime = System.nanoTime()
        var seqscanTime = 0L

        Utils.getBufferedWriter(filePath).use { myWriter ->
            //start with haplistid header
            myWriter.write("#haplotypeListId=$haplotypeListId\n")

            for(refrangeId in refRangeToKmerSetMap.keys) {
                if(rangeCount%1000 == 0) {
                    myLogger.info("Time Spent Processing output: ${(System.nanoTime() - startTime)/1e9} seconds.  Processed ${rangeCount} Ranges.")
                }
                rangeCount++
                val kmers = refRangeToKmerSetMap[refrangeId]?:throw IllegalStateException("Error no kmers associated with this refRange ID: ${refrangeId}")

                val mapForRange = kmers.map { Pair(it, kmerMapToHapIds[it]) }.toMap()

                //generate map of hapidSet -> kmer hash list
                val hapidKmerHashMap = mutableMapOf<MutableSet<Int>, MutableList<Long>>()
                for (entry in mapForRange.entries) {
                    val kmerHashList = hapidKmerHashMap[entry.value]
                    if (kmerHashList == null) hapidKmerHashMap.put(entry.value.toMutableSet(), mutableListOf(entry.key))
                    else kmerHashList.add(entry.key)
                }

                val numberOfHaplotypesInRange = refRangeToHapidMap[refrangeId]?.size?:throw IllegalStateException("Error no haplotypes associated with this refRange.")
                val numberOfHapSets = hapidKmerHashMap.size
                if (numberOfHapSets == 0) continue  //skip to next refrange because there are no kmers here
                val numberOfBitsNeeded = numberOfHaplotypesInRange * numberOfHapSets
                //encodedHapSets is a BitSet containing all of the Hapsets seen in this reference range
                //each hapset is encoded with a series of bits equal to the number of haplotypes in the range
                //with the bits in this hapset set
                //the offset for each kmer is the offset for its associated hapset
                val encodedHapSets = BitSet(numberOfBitsNeeded)
                //make pairs of kmerHash, offset
                val kmerHashOffsets = mutableListOf<Pair<Long,Int>>()

                var offset = 0
                val hapidIndex = refRangeToHapidMap[refrangeId]
                check(hapidIndex != null) {"hapidIndex null for $refrangeId"}

                hapidKmerHashMap.entries.forEachIndexed { index, entry ->
                    //this line encodes a haplotype set (hapset)
                    for (hapid in entry.key) encodedHapSets.set(offset + hapidIndex[hapid]!!)
                    //this line stores a pair of kmerHash, offset for each kmer mapping to this haplotype set
                    for (kmerHash in entry.value) kmerHashOffsets.add(Pair(kmerHash,offset))
                    offset += numberOfHaplotypesInRange
                }

                //write the results of the range to a file
                //line 1: >rangeid
                myWriter.write(">${refrangeId}\n")
                //line 2: long1,long2,...,longn (range has hapid sets encoded into a bitSet, which can be stored as longs)
                myWriter.write("${encodedHapSets.toLongArray().joinToString(",")}\n")
                //line 3: hash1,offset1,hash2,offset2,...,hashn,offsetn
                myWriter.write("${kmerHashOffsets.map { "${it.first}@${it.second}" }.joinToString(",")}\n")

            }
        }

        myLogger.info("Saved kmer mapping to $filePath, elapsed time ${(System.nanoTime() - startTime)/1e9} sec, seqscan time ${seqscanTime/1e9} sec")
    }

    companion object {
        //Coding for nucleotides is A -> 0, C -> 1, G -> 2, T -> 3
        //Generated by nuc_char shr 1 and 3
        val Avalue: Long = 0L
        val Cvalue: Long = 1L
        val Gvalue: Long = 2L
        val Tvalue: Long = 3L
        val revCompA = Tvalue shl 62
        val revCompC = Gvalue shl 62
        val revCompG = Cvalue shl 62
        val revCompT = Avalue shl 62

        /**
         * Updates hashes with the next nucleotide for both the kmer sequence and its reverse compliment.
         */
        fun updateKmerHashAndReverseCompliment(hashes: Pair<Long,Long>, nucleotide: Char): Pair<Long,Long> {
            return when (nucleotide) {
                'A' -> Pair((hashes.first shl 2) or Avalue, (hashes.second ushr 2) or revCompA)
                'T' -> Pair((hashes.first shl 2) or Tvalue, (hashes.second ushr 2) or revCompT)
                'G' -> Pair((hashes.first shl 2) or Gvalue, (hashes.second ushr 2) or revCompG)
                'C' -> Pair((hashes.first shl 2) or Cvalue, (hashes.second ushr 2) or revCompC)
                else -> throw java.lang.IllegalArgumentException("Attempted to update kmer hash with an invalid nucleotide character($nucleotide). Must be one of A,G,C,T")
            }
        }

        /**
         * data class to hold the KmerMap information.
         */
        data class KmerMapData(val haplotypeListId: Int, val rangeToBitSetMap: Map<Int, BitSet>, val kmerHashToLongMap: Long2LongOpenHashMap)

        /**
         * Loads a kmer hash map file to memory for use in mapping reads.
         */
        fun loadKmerMaps(filename: String): KmerMapData {
            //Load the contents of the file into
            //rangeHapidMap: a map of refRangeId to an IntArray of the haplotype ids in the ReferenceRange
            // and the BitSet of all hapid sets
            //kmerHashmap: a map of kmer hash to reference range and offset into its BitSet, encoded as a long
            //These data structures all the reference range and haplotype set to be looked up for a kmer has

            val rangeToBitSetMap = mutableMapOf<Int, BitSet>()
            val kmerHashMap = Long2LongOpenHashMap()
            var lineCount = 0
            var totalLineCount = 0
            var refrangeId = 0
            var haplotypeListId = -1

            Utils.getBufferedReader(filename).useLines {
                it.forEach { inputStr ->
                    totalLineCount++

                    lineCount = when(inputStr.first()) {
                        '>' -> 1
                        '#' -> 4
                        else -> lineCount + 1
                    }

                    when (lineCount) {
                        1 -> {
                            refrangeId = inputStr.substring(1).toInt()                  }
                        2-> {
                            try {
                                val parsedLine = inputStr.split(",")
                                val myBitset = BitSet.valueOf(parsedLine.map { it.toLong() }.toLongArray())
                                rangeToBitSetMap[refrangeId] = myBitset
                            } catch(e: Exception) {
                                println("error at line $totalLineCount for input = $inputStr")
                                throw java.lang.IllegalArgumentException(e)
                            }
                        }
                        3 -> {
                            val refrangeLong = (refrangeId.toLong()) shl 32
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
                        4 -> {
                            if (inputStr.startsWith("#haplotypeListId="))
                                haplotypeListId = inputStr.substringAfter('=', "-1").toInt()

                        }
                    }
                }
            }

            return KmerMapData(haplotypeListId, rangeToBitSetMap, kmerHashMap)
        }
    }
}

}