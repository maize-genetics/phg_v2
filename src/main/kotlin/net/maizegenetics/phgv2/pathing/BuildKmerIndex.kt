package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.long
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.retrieveAgcContigs
import org.apache.logging.log4j.LogManager
import java.io.BufferedWriter
import java.io.File
import java.io.FileOutputStream
import java.lang.IllegalStateException
import java.nio.file.Files
import java.nio.file.Paths
import java.util.*
import java.util.zip.GZIPOutputStream


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

    val dbPath by option(help = "Tile DB URI")
        .required() //Needs to be required now due to the agc archive

    val maxHaplotypeProportion by option("-p", "--maxHapProportion", help = "only kmers mapping to less than or " +
            "equal to maxHapProportion of haplotypes in a reference range will be retained.")
        .double()
        .default(0.75)

    val hashMask by option("-m", "--hashMask", help = "with hashFilter, used to mask kmers for filtering. " +
            "Default uses only the last kmer nucleotide. Only change this if you know what you are doing.")
        .long()
        .default(3)

    val hashFilterValue by option("-f", "--hashFilter", help = "Only hashes that pass the filter" +
            " ((hashValue and hashMask) == hashFilter) will be considered")
        .long()
        .default(1)

    val hvcfDir by option("--hvcf-dir", help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        .required()//Todo: make this optional by adding .default("")



    override fun run() {
        //either tiledbPath or agcPath must be provided


        //build the haplotypeGraph
        val graph = buildHaplotypeGraph()
        val hashToHapidMap = processGraphKmers(graph, dbPath, maxHaplotypeProportion,  hashMask, hashFilterValue)

        //for now, the name of the kmerIndex will be kmerIndex.txt. Later, the file path and name can be set by the user.
        val kmerIndexFilename = "${hvcfDir}kmerIndex.txt"

        //save the kmerIndex
        saveKmerHashesAndHapids(graph, kmerIndexFilename, hashToHapidMap)
    }

    private fun buildHaplotypeGraph(): HaplotypeGraph {
        val timedValue = measureTimedValue {
            if(hvcfDir != "") {
                val pathList = File(hvcfDir).listFiles { file -> file.name.endsWith(".h.vcf") || file.name.endsWith(".h.vcf.gz") }.map { it.path }
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
    fun processGraphKmers(graph: HaplotypeGraph, dbPath: String, maxHaplotypeProportion: Double=.75, hashMask: Long = 3, hashFilterValue:Long = 1) : Long2ObjectOpenHashMap<Set<String>> {
        //keepMap is a map of hash -> Set of haplotype ids
        val keepMap = Long2ObjectOpenHashMap<Set<String>>()
        //discardSet is a Set of hashes
        val discardSet = LongOpenHashSet()
        var rangeCount = 0
        val startTime = System.nanoTime()

        for (refrange in graph.ranges()) {
            if (rangeCount++ % 5000 == 0) {
                myLogger.debug("Processing range $rangeCount, keep set size = ${keepMap.size}, discard set size = ${discardSet.size}, elapse time ${(System.nanoTime() - startTime)/1e9} sec")
            }

            val hapidToSampleMap = graph.hapIdToSampleGametes(refrange)
            val maxHaplotypes = ceil(hapidToSampleMap.size * maxHaplotypeProportion)

            //create a map of hash -> count of occurrences for all the haplotypes in this reference range
            //retrieveAgcContigs needs a list of contig@genome:start-end, for all of the regions in the haplotype alt headers
            //subtract 1 from altHeader positions because AGC positions are 0-based.
            //Also need a map of Pair(sampleName, region) -> hapid in order to associate sequence from agc to the originating hapid.
            //The region has to match the region as requested from agc.
            val agcRangeList = mutableListOf<String>()
            val hapidSampleRangeList = mutableListOf<Array<String>>()

            hapidToSampleMap.keys.forEach { hapid ->
                val alt = graph.altHeader(hapid) ?: throw IllegalStateException("No alt header for $hapid")

                for (range in alt.regions) {
                    if (range.first.position <= range.second.position) {
                        agcRangeList.add("${range.first.contig}@${alt.sampleName()}:${range.first.position - 1}-${range.second.position - 1}")
                        hapidSampleRangeList.add(arrayOf(hapid, alt.sampleName(), "${range.first.contig}:${range.first.position - 1}-${range.second.position - 1}"))
                    } else {
                        agcRangeList.add("${range.first.contig}@${alt.sampleName()}:${range.second.position - 1}-${range.first.position - 1}")
                        hapidSampleRangeList.add(arrayOf(hapid, alt.sampleName(), "${range.first.contig}:${range.second.position - 1}-${range.first.position - 1}"))
                    }
                }

            }

            val sampleRegionToHapidMap = hapidSampleRangeList.groupBy({Pair(it[1], it[2])},{it[0]})

            val agcdataMap = retrieveAgcContigs(dbPath, agcRangeList)
            //retrieveAgcContigs returns a Map<Pair<String,String>,NucSeq> where Pair is
            // sampleName, contig:start-end and NucSeq is the sequence
            // The following gets the hapid that was used to generate each agc entry and creates a map of hapid -> list of sequences
            // from that hapid.
            val hapidSeqMap = agcdataMap.entries.mapNotNull { (sample, seqrec) ->
                val hapidList = sampleRegionToHapidMap[sample]
                hapidList?.map { Pair(it, seqrec) }
            }.flatten().groupBy ({it.first}, {it.second.seq()})

            val (kmerHashCounts, longToHapIdMap) = countKmerHashesForHaplotypeSequence(hapidSeqMap, hashMask, hashFilterValue)

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

    /**
     * Function to count the kmerHashes for a single reference range's haplotype nodes.
     * This gets put in 2 sets of maps.
     * One for the hash counts and one for a hash to a list of hapIds which contain that hash.
     * The sequenceList input is a map of hapidId -> list of sequences
     */
    private fun countKmerHashesForHaplotypeSequence(sequenceMap: Map<String, List<String>>, hashMask: Long, hashFilterValue: Long) : Pair<Map<Long,Int>,Map<Long,Set<String>>> {
        //start by splitting sequence into subsequences without N's
        //mapOfHashes is a map of hash -> count of occurrences
        val mapOfHashes = Long2IntOpenHashMap()
        //mapOfHapIds is a map of hash -> Set of haplotype ids
        val mapOfHapIds = Long2ObjectOpenHashMap<MutableSet<String>>()

        for ((hapid, seqList) in sequenceMap) {
            for (sequenceWithNs in seqList) {
                //split sequence on N's then filter on length > 31 because we are looking for 32-mers
                val splitList = sequenceWithNs
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


        }
        return Pair(mapOfHashes,mapOfHapIds)
    }

    /**
     * Function to save the kmer hash map to the specified file.
     */
    fun saveKmerHashesAndHapids(graph: HaplotypeGraph, filePath: String, kmerMapToHapIds: Long2ObjectOpenHashMap<Set<String>>) {
        //refRangeToHapidMap is a map of refRange -> (map of hapid -> index of hapid)
        val refRangeToHapIndexMap = graph.refRangeToHapIdMap()

        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()


        //Walk through the map and associate kmerLongs (the kmer hashes) with specific referenceRanges
        //refRangeToKmerSetMap is a map of refRange -> Set of kmer hashes
        val refRangeToKmerSetMap = getRefRangeToKmerSetMap(kmerMapToHapIds, hapIdToRefRangeMap)

        // graph persistence not needed for now, will write index to hvcfDir
        //TODO add sufficient information to be able to reproduce the graph from the index file
        //   for example a list of the hvcf files or tiledb sample names

        val startTime = System.nanoTime()

        getBufferedWriter(filePath).use { myWriter ->
            //do NOT start with haplistid header

            for((rangeCount, refrange) in refRangeToKmerSetMap.keys.withIndex()) {

                if(rangeCount % 1000 == 0) {
                    myLogger.info("Time Spent Processing output: ${(System.nanoTime() - startTime)/1e9} seconds.  Processed $rangeCount Ranges.")
                }
                //Do extract out the kmers for the refRange and export to file
                extractKmersAndExportIndexForRefRange(refRangeToKmerSetMap, refrange, kmerMapToHapIds, refRangeToHapIndexMap, myWriter)
            }
        }

        myLogger.info("Saved kmer mapping to $filePath, elapsed time ${(System.nanoTime() - startTime)/1e9} sec")
    }

    fun extractKmersAndExportIndexForRefRange(
        refRangeToKmerSetMap: Map<ReferenceRange, Set<Long>>,
        refrange: ReferenceRange,
        kmerMapToHapIds: Long2ObjectOpenHashMap<Set<String>>,
        refRangeToHapIndexMap: Map<ReferenceRange, Map<String, Int>>,
        myWriter: BufferedWriter
    ) {
        //Check to make sure the refRange maps have the correct keys.
        if(!refRangeToKmerSetMap.containsKey(refrange)) {
            myLogger.warn("No kmers associated with this refRange: $refrange")
            return
        }
        check(refRangeToHapIndexMap.containsKey(refrange)) { "Error no haplotypes associated with this refRange: $refrange" }

        val kmers = refRangeToKmerSetMap[refrange]!!

        //mapForRange is a map of kmer hash -> Set of hapid ids
        val mapForRange = kmers.associateWith { kmerMapToHapIds[it] }

        //generate map of hapidSet -> kmer hash list
        val hapidKmerHashMap = createHapIdToKmerMap(mapForRange)

        val numberOfHaplotypesInRange = refRangeToHapIndexMap[refrange]!!.size

        val numberOfHapSets = hapidKmerHashMap.size
        if (numberOfHapSets == 0) return  //skip to next refrange because there are no kmers here
        val (encodedHapSets, kmerHashOffsets) = buildEncodedHapSetsAndHashOffsets(
            numberOfHaplotypesInRange,
            numberOfHapSets,
            refRangeToHapIndexMap,
            refrange,
            hapidKmerHashMap
        )

        exportRefRangeKmerIndex(myWriter, refrange, encodedHapSets, kmerHashOffsets)
    }

    fun buildEncodedHapSetsAndHashOffsets(
        numberOfHaplotypesInRange: Int,
        numberOfHapSets: Int,
        refRangeToHapIndexMap: Map<ReferenceRange, Map<String, Int>>,
        refrange: ReferenceRange,
        hapidKmerHashMap: Map<Set<String>, List<Long>>
    ): Pair<BitSet, MutableList<Pair<Long, Int>>> {
        val numberOfBitsNeeded = numberOfHaplotypesInRange * numberOfHapSets
        //encodedHapSets is a BitSet containing all of the Hapsets seen in this reference range
        //each hapset is encoded with a series of bits equal to the number of haplotypes in the range
        //with the bits in this hapset set
        //the offset for each kmer is the offset for its associated hapset
        val encodedHapSets = BitSet(numberOfBitsNeeded)
        //make pairs of kmerHash, offset
        val kmerHashOffsets = mutableListOf<Pair<Long, Int>>()

        var offset = 0
        val hapidIndex = refRangeToHapIndexMap[refrange]
        check(hapidIndex != null) { "hapidIndex null for $refrange" }

        hapidKmerHashMap.entries.forEachIndexed { index, entry ->
            //this line encodes a haplotype set (hapset)
            for (hapid in entry.key) encodedHapSets.set(offset + hapidIndex[hapid]!!)
            //this line stores a pair of kmerHash, offset for each kmer mapping to this haplotype set
            for (kmerHash in entry.value) kmerHashOffsets.add(Pair(kmerHash, offset))
            offset += numberOfHaplotypesInRange
        }
        return Pair(encodedHapSets, kmerHashOffsets)
    }

    fun exportRefRangeKmerIndex(
        myWriter: BufferedWriter,
        refrange: ReferenceRange,
        encodedHapSets: BitSet,
        kmerHashOffsets: List<Pair<Long, Int>>
    ) {
        //write the results of the range to a file
        //line 1: >rangeid
        myWriter.write(">$refrange\n")
        //line 2: long1,long2,...,longn (range has hapid sets encoded into a bitSet, which can be stored as longs)
        myWriter.write("${encodedHapSets.toLongArray().joinToString(",")}\n")
        //line 3: hash1,offset1,hash2,offset2,...,hashn,offsetn
        myWriter.write("${kmerHashOffsets.map { "${it.first}@${it.second}" }.joinToString(",")}\n")
    }

    /** Function to reverse the multimap from kmerHash -> Set of hapid ids to hapidSet -> List of kmer hashes
     *  This is needed to be able to encode the haplotype sets into a bitset
     */
    fun createHapIdToKmerMap(kmerMapForRange: Map<Long, Set<String>>): Map<Set<String>, List<Long>> {
        val hapidKmerHashMap = mutableMapOf<Set<String>, MutableList<Long>>()
        for (entry in kmerMapForRange.entries) {
            val kmerHashList = hapidKmerHashMap[entry.value]
            if (kmerHashList == null) hapidKmerHashMap.put(entry.value.toMutableSet(), mutableListOf(entry.key))
            else kmerHashList.add(entry.key)
        }
        return hapidKmerHashMap
    }

    /**
     * Function to create a map of refRange -> Set of kmer hashes based on the kmerMapToHapIds and hapIdToRefRangeMap
     */
    fun getRefRangeToKmerSetMap(
        kmerMapToHapIds: Long2ObjectOpenHashMap<Set<String>>,
        hapIdToRefRangeMap: MutableMap<String, ReferenceRange>
    ): Map<ReferenceRange, Set<Long>> {
        val refRangeToKmerSetMap = mutableMapOf<ReferenceRange, MutableSet<Long>>()
        for (kmerMap in kmerMapToHapIds.long2ObjectEntrySet()) {
            val kmer = kmerMap.longKey
            val currentRefRange = hapIdToRefRangeMap[kmerMap.value.first()]!!
            if (refRangeToKmerSetMap.containsKey(currentRefRange)) {
                refRangeToKmerSetMap[currentRefRange]!!.add(kmer)
            } else {
                refRangeToKmerSetMap[currentRefRange] = mutableSetOf(kmer)
            }
        }
        return refRangeToKmerSetMap
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

        fun getBufferedWriter(filePath: String): BufferedWriter {
            return if (filePath.endsWith(".gz")) {
                GZIPOutputStream(FileOutputStream(filePath)).bufferedWriter()
            } else {
                Files.newBufferedWriter(Paths.get(filePath))
            }
        }
    }
}
