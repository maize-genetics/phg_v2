package net.maizegenetics.phgv2.pathing

import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.int
import com.github.ajalt.clikt.parameters.types.long
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.*
import org.apache.logging.log4j.LogManager
import java.io.BufferedWriter
import java.io.File
import java.nio.file.Files
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
class BuildKmerIndex: CliktCommand(help="Create a kmer index for a HaplotypeGraph. By default the file will be written " +
        "to <hvcfDir>/kmerIndex.txt") {

    private val myLogger = LogManager.getLogger(BuildKmerIndex::class.java)

    val dbPath by option(help = "Tile DB URI")
        .required() //Needs to be required now due to the agc archive

    val indexFile by option(help = "The full path of the kmer index file. Default = <hvcf-dir>/kmerIndex.txt")
        .default("")

    val maxHaplotypeProportion by option("-p", "--maxHapProportion", help = "only kmers mapping to less than or " +
            "equal to maxHapProportion of haplotypes in a reference range will be retained.Default = 0.75")
        .double()
        .default(0.75)

    val hashMask by option("-m", "--hashMask", help = "with hashFilter, used to mask kmers for filtering. " +
            "Default uses only the last kmer nucleotide. Only change this if you know what you are doing. Default = 3")
        .long()
        .default(3)

    val hashFilterValue by option("-f", "--hashFilter", help = "Only hashes that pass the filter" +
            " ((hashValue and hashMask) == hashFilter) will be considered. Do not change this value unless you know " +
            "what you are doing. Default = 1")
        .long()
        .default(1)

    val hvcfDir by option("--hvcf-dir", help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        .required()//Todo: make this optional by adding .default("")

    val maxArgLength by option(help="The maximum argument length for a call to agc. This defaults to 200000. " +
            "If you get an error caused by a call to agc being too long try reducing this value.")
        .int()
        .default(200000)

    val noDiagnostics by option("-n", "--no-diagnostics", help = "Flag that will suppress writing of diagnostics.").flag()

    private val refrangeToAdjacentHashCount = mutableMapOf<ReferenceRange, Int>()
    private var runDiagnostics = true

    override fun run() {
        //build the haplotypeGraph
        myLogger.info("Start of BuildKmerIndex...")
        val graph = buildHaplotypeGraph()

        val hashToHapidMap = processGraphKmers(graph, dbPath, maxHaplotypeProportion,  hashMask, hashFilterValue)

        val kmerIndexFilename = if (indexFile == "") "${hvcfDir}/kmerIndex.txt" else indexFile

        //save the kmerIndex
        saveKmerHashesAndHapids(graph, kmerIndexFilename, hashToHapidMap)

        if (noDiagnostics) {
            runDiagnostics = false
            myLogger.info("BuildKmerIndex: Diagnostic output will not be written because the --no-diagnostic flag was set.")
        }
        else writeDiagnostics(refrangeToAdjacentHashCount)

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
    fun processGraphKmers(graph: HaplotypeGraph, dbPath: String, maxHaplotypeProportion: Double=.75,
                                        hashMask: Long = 3, hashFilterValue:Long = 1) : Long2ObjectOpenHashMap<Set<String>> {
        //keepMap is a map of hash -> Set of haplotype ids
        val keepMap = Long2ObjectOpenHashMap<Set<String>>()

        //discardSet is a Set of hashes
        val discardSet = LongOpenHashSet()
        val startTime = System.nanoTime()
        val sampleGametes = graph.sampleGametesInGraph()

        val contigRangesMap = graph.rangesByContig()

        //needed for diagnostics
        val refRangeToIndexMap = graph.refRangeToIndexMap()
        val hapidToRefrangeMap = graph.hapIdToRefRangeMap()

        val maxHapsToKeep = sampleGametes.size * 2

        for (chr in contigRangesMap.keys) {
            //get all sequence for this chromosome
            val agcRequestLists = rangeListsForAgcCommand(graph, contigRangesMap, chr)
            val agcChromSequence = if (agcRequestLists.sampleContigList.isNotEmpty()) {
                myLogger.info("sampleContigList size = ${agcRequestLists.sampleContigList.size}, first element = ${agcRequestLists.sampleContigList[0]}")
                retrieveAgcContigs(dbPath, agcRequestLists.sampleContigList)
            }
            else emptyMap()
            val agcOtherRegionSequence: Map<Pair<String,String>, NucSeq> = if (agcRequestLists.otherRegionsList.isNotEmpty()) getAgcSequenceForRanges(agcRequestLists.otherRegionsList)
            else emptyMap()

            //iterate through refranges, generate kmers from the sequence
            for (refrange in contigRangesMap[chr]!!) {
                val hapidToSampleMap = graph.hapIdToSampleGametes(refrange)

                //if there are any null haplotypes add 1 to the number of hapids (hapidToSampleMap.size)
                //The reason is that if a read maps to all non-null haplotypes, it is still informative if there are null haplotypes
                val numberOfSampleGametesWithHapids = hapidToSampleMap.values.sumOf { it.size }
                val hasNullHaplotypes = numberOfSampleGametesWithHapids < sampleGametes.size
                val maxHaplotypes = if (hasNullHaplotypes) ceil((hapidToSampleMap.size + 1) * maxHaplotypeProportion)
                else ceil(hapidToSampleMap.size * maxHaplotypeProportion)

                //map haplotype ids to the sequence for that hapid
                val hapidToSequencMap = mutableMapOf<String, List<String>>()
                for (hapid in hapidToSampleMap.keys) {
                    //checked for altHeader existence already
                    val altHeader = graph.altHeader(hapid)!!
                    val sequenceList = altHeader.regions.map { region ->
                        val inThisChrom = region.first.contig == chr
                        if (inThisChrom) {
                            //translate from 1-based Position to 0-based nucseq coordinates
                            val seqRange = if (region.first.position <= region.second.position) (region.first.position - 1..region.second.position - 1)
                            else (region.second.position - 1..region.first.position - 1)
                            agcChromSequence[Pair(altHeader.sampleName(), chr)]?.get(seqRange)?.seq() ?:""
                        } else {
                            agcOtherRegionSequence[Pair(altHeader.sampleName(), regionToString(region))]?.seq() ?:""
                        }
                    }
                    hapidToSequencMap[hapid] = sequenceList
                }

                val (kmerHashCounts, longToHapIdMap) = countKmerHashesForHaplotypeSequence(hapidToSequencMap, hashMask, hashFilterValue)
                for (hashCount in kmerHashCounts.entries) {
                    val hashValue = hashCount.key

                    //if the hash is in the discard set skip it
                    if (discardSet.contains(hashValue)) continue

//                    when {
//                        //if hash count >= numberOfHaplotype add it to the discard set
//                        hashCount.value >= maxHaplotypes -> discardSet.add(hashValue)
//                        //if the hash is already in the keepSet, it has been seen in a previous reference range
//                        //was this hash seen in the range immediately preceeding this one?
//                        // so, remove it from the keep set and add it to the discard set
//                        keepMap.containsKey(hashValue) -> {
//                            val hapidSet = keepMap.remove(hashValue)
//                            if (runDiagnostics) {
//                                val hapidRefrange = hapidToRefrangeMap[hapidSet.first()]
//                                val wasPreviousRange = refRangeToIndexMap[refrange] == ((refRangeToIndexMap[hapidRefrange] ?: -2) + 1)
//                                if (wasPreviousRange) {
//                                    val oldCount = refrangeToAdjacentHashCount.getOrElse(refrange) {0}
//                                    refrangeToAdjacentHashCount[refrange] = oldCount + 1
//                                }
//                            }
//                            discardSet.add(hashValue)
//                        }
//                        else -> {
//                            keepMap[hashValue] = longToHapIdMap[hashValue]
//                        }
//                    }

                    when {
                        //if hash count >= numberOfHaplotype add it to the discard set
                        hashCount.value >= maxHaplotypes -> discardSet.add(hashValue)
                        //if the hash is already in the keepSet, it has been seen in a previous reference range
                        //was this hash seen in the range immediately preceeding this one?
                        // so, remove it from the keep set and add it to the discard set
                        keepMap.containsKey(hashValue) && (keepMap[hashValue].size + longToHapIdMap[hashValue]!!.size) > maxHapsToKeep -> {
                            keepMap.remove(hashValue)
//=======
//                        keepMap.containsKey(hashValue) -> {
//                            val hapidSet = keepMap.remove(hashValue)
//                            if (runDiagnostics) {
//                                //check whether kmer was seen in previous refrange
//                                //hapidRefrangeIndexSet is the set of refrange indices of the refranges contain the hapids in hapidSet
//                                //that is the refranges that contain this kmer (mostly but not always a single refrange)
//                                val hapidRefrangeIndexSet = hapidSet.mapNotNull { hapidToRefrangeMap[it] }.flatten()
//                                    .map { refRangeToIndexMap[it] }.toSet()
//                                val wasPreviousRange = hapidRefrangeIndexSet.contains((refRangeToIndexMap[refrange] ?: 0)  - 1)
//                                if (wasPreviousRange) {
//                                    val oldCount = refrangeToAdjacentHashCount.getOrElse(refrange) {0}
//                                    refrangeToAdjacentHashCount[refrange] = oldCount + 1
//                                }
//                            }
//>>>>>>> main
                            discardSet.add(hashValue)
                        }
                        keepMap.containsKey(hashValue) && (keepMap[hashValue].size + longToHapIdMap[hashValue]!!.size) <= maxHapsToKeep -> {
                            keepMap[hashValue] = keepMap[hashValue]!!.union(longToHapIdMap[hashValue]!!)
                        }
                        else -> {
                            keepMap[hashValue] = longToHapIdMap[hashValue]
                        }
                    }
                }

            }

        }
        myLogger.debug("Finished building kmer keep set, keep set size = ${keepMap.size}, discard set size = ${discardSet.size}, elapse time ${(System.nanoTime() - startTime)/1e9} sec")

        return keepMap
    }

    private fun writeDiagnostics(adjacentHashCounts: Map<ReferenceRange, Int>) {
        val diagnosticFileName = "kmerIndexStatistics.txt"
        val diagnosticFilePath = if (indexFile.isBlank()) {
            File(hvcfDir).resolve(diagnosticFileName).absolutePath
        } else {
            File(indexFile).parentFile.resolve(diagnosticFileName).absolutePath
        }

        val tmpStats = Files.createTempFile("stats", ".txt").toString()

        //start by running KmerIndexStatistics
        val argList = listOf("--hvcf-dir", hvcfDir, "--index-file", indexFile, "--output-file", tmpStats)
        KmerIndexStatistics().main(argList)

        //Add the counts from adjacentHashCounts
        getBufferedReader(tmpStats).use { myReader ->
            getBufferedWriter(diagnosticFilePath).use { myWriter ->
                val header = myReader.readLine()
                myWriter.write("$header\tadjacentCount\n")
                var inputLine = myReader.readLine()
                while (inputLine  != null) {
                    val data = inputLine.split("\t")
                    val refrange = ReferenceRange.parse("${data[0]}:${data[1]}-${data[2]}")
                    val adjacentCount = adjacentHashCounts.getOrElse(refrange){0}
                    myWriter.write("$inputLine\t${adjacentCount}\n")
                    inputLine = myReader.readLine()
                }
            }
        }

    }

    private fun regionToString(region: Pair<Position,Position>): String {
        return if (region.first.position <= region.second.position) {
            "${region.first.contig}:${region.first.position - 1}-${region.second.position - 1}"
        } else {
            "${region.first.contig}:${region.second.position - 1}-${region.first.position - 1}"
        }
    }

    private fun getAgcSequenceForRanges(ranges: List<String>): Map<Pair<String,String>, NucSeq> {
        val argLength = ranges.sumOf { it.length } / ranges.size + 1
        val maxArgs = maxArgLength / argLength
        val sequenceMap = mutableMapOf<Pair<String,String>, NucSeq>()
        ranges.windowed(maxArgs, maxArgs, true).forEach {
            myLogger.debug("getting sequence for $it")
            sequenceMap.putAll(retrieveAgcContigs(dbPath, it))
        }
        return sequenceMap
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

        val startTime = System.nanoTime()

        getBufferedWriter(filePath).use { myWriter ->

            for((rangeCount, refrange) in refRangeToKmerSetMap.keys.withIndex()) {

                if(rangeCount % 1000 == 0) {
                    myLogger.info("Time Spent Processing output: ${(System.nanoTime() - startTime)/1e9} seconds.  Processed $rangeCount Ranges.")
                }
                //extract out the kmers for the refRange and export to file
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
            for (hapid in entry.key) {
                if(hapidIndex.containsKey(hapid)) encodedHapSets.set(offset + hapidIndex[hapid]!!)
                else myLogger.warn("BuildKmerIndex.buildEncodedHapSetsAndHashOffsets: ndx = null for hapid = $hapid.")
            }

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
        hapIdToRefRangeMap: Map<String, List<ReferenceRange>>
    ): Map<ReferenceRange, Set<Long>> {
        val refRangeToKmerSetMap = mutableMapOf<ReferenceRange, MutableSet<Long>>()
        for (kmerMap in kmerMapToHapIds.long2ObjectEntrySet()) {
            val kmer = kmerMap.longKey
            val hapidSet = kmerMap.value

            //use the most frequent reference range
            //all the haplotype ids map to the same refrange
            //rarely some kmers will map to additional ref ranges but not all
            //the most frequent ReferenceRange will be the one used to create this hapid set
            val referenceRangeList = hapidSet.mapNotNull { hapIdToRefRangeMap[it] }.flatten()
            val referenceRangeCounts = referenceRangeList.groupingBy { it }.eachCount()
            val currentRefRange = referenceRangeCounts.maxBy { it.value }.key
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

        /**
         * sampleContigList is a list of Strings of the form contig@sample. otherRegionsList is a list
         * of Strings of the form contig@sample:xxx-yyy.
         */
        data class AgcLists(val sampleContigList: List<String>, val otherRegionsList: List<String>)

        /**
         * @param graph a [HaplotypeGraph]
         * @param contigRangesMap   a map of contigs -> list of [ReferenceRange] for that contig generated by graph.rangesByContig()
         * @param chr a chromosome (contig)
         *
         * @return  an [AgcLists] containing lists that can be used as part of an agc command
         *
         * This function uses the alt header regions to make lists of contigs or ranges that can be used as part of
         * an agc command to get all the sequence for a chromosome from a [HaplotypeGraph].
         */
        fun rangeListsForAgcCommand(graph: HaplotypeGraph, contigRangesMap: Map<String, List<ReferenceRange>>, chr: String): AgcLists {
            require(contigRangesMap.containsKey(chr)) {"In BuildKmerIndex.rangeListsForAgcCommand() chr is not in contigRangesMap."}
            //If at least half the ranges in a chromosome come from the sample contig, add the contig to the list of
            //  contig@sample to be retrieved. Otherwise add the individual ranges to otherRegions.
            val minRangeCount = 0.5 * contigRangesMap[chr]!!.size

            //sampleNameToRegionsMap is a map of sample -> (map of contig -> set of ranges in that contig)
            val sampleNameToRegionsMap = mutableMapOf<String, MutableMap<String,MutableSet<String>>>()
            for (refrange in contigRangesMap[chr]!!) {
                for (hapid in graph.hapIdToSampleGametes(refrange).keys) {
                    val altheader = graph.altHeader(hapid)
                    check(altheader != null) {"altheader null for $hapid in $refrange"}
                    val sampleName = altheader.sampleName()
                    //regionsMap is a map of contig -> set of ranges in that contig
                    val regionsMap = sampleNameToRegionsMap.getOrPut(sampleName){mutableMapOf()}
                    val regions = altheader.regions

                    //for each range in region, that range is added to regionsMap for the sample
                    for (region in regions) {
                        val rangeSet = regionsMap.getOrPut(region.first.contig) { mutableSetOf() }
                        val rangeStr = if (region.first.position <= region.second.position) "${region.first.position}-${region.second.position}"
                        else "${region.second.position}-${region.first.position}"
                        rangeSet.add(rangeStr)
                    }
                }
            }

            //For each sample name, remove from its regionsMap the contig with the highest range count if the count exceeds minRangeCount
            //For each sample-contig combination removed, add the string "contig@sampleName" to the sampleContigList
            val sampleContigList = mutableListOf<String>()
            val otherRegions = mutableListOf<String>()
            for ((sampleName, regionsMap) in sampleNameToRegionsMap.entries) {
                val biggestEntry = regionsMap.entries.maxBy { (_, regionSet) -> regionSet.size }
                if (biggestEntry.value.size > minRangeCount) {
                    regionsMap.remove(biggestEntry.key)
                    sampleContigList.add("${biggestEntry.key}@$sampleName")
                }
                //add remaining entries to otherRegions
                for ((contig, rangeSet) in regionsMap) {
                    for (range in rangeSet) otherRegions.add("$contig@$sampleName:$range")
                }
            }

            return AgcLists(sampleContigList, otherRegions)
        }

    }
}
