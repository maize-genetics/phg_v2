package net.maizegenetics.phgv2.pathing

import biokotlin.seq.NucSeq
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.int
import com.github.ajalt.clikt.parameters.types.long
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.longs.LongOpenHashBigSet
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.utils.*
import org.apache.logging.log4j.LogManager
import java.io.BufferedWriter
import java.io.File
import java.nio.file.Files
import java.util.*


import kotlin.math.min
import kotlin.time.DurationUnit
import kotlin.time.measureTimedValue

/**
 * Creates a Map of 32-mer hash -> hapid list for the haplotypes in a HaplotypeGraph. Hashes will be retained only if they map to at
 * most [maxHaplotypeProportion] * the number of haplotyes in a reference range.  If a kmer is too repetitive(>2 * numSamples), the kmer will be purged.
 * Only hashes that pass the filter ((hashValue and [hashMask]) == [hashFilterValue]) will be considered.
 * For example, setting [hashMask] = 3u and [hashFilterValue] = 1u
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

    val discardFile by option(help = "The full path of the discard set file. If Left blank no file will be written out.")
        .default("")

    val maxHaplotypeProportion by option("-p", "--maxHapProportion", help = "only kmers mapping to less than or " +
            "equal to maxHapProportion of haplotypes in a reference range will be retained.")
        .double()
        .default(0.75)

    val hashMask by option("-m", "--hashMask", help = "with hashFilter, used to mask kmers for filtering. " +
            "Default uses only the last kmer nucleotide. Only change this if you know what you are doing.")
        .long()
        .default(3)

    val hashFilterValue by option("-f", "--hashFilter", help = "Only hashes that pass the filter" +
            " ((hashValue and hashMask) == hashFilter) will be considered. Do not change this value unless you know " +
            "what you are doing.")
        .long()
        .default(1)

    val hvcfDir by option("--hvcf-dir", help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        .required()//Todo: make this optional by adding .default("")

    val maxArgLength by option(help="The maximum argument length for a call to agc." +
            "If you get an error caused by a call to agc being too long try reducing this value.")
        .int()
        .default(200000)

    val noDiagnostics by option("-n", "--no-diagnostics", help = "Flag that will suppress writing of diagnostics.").flag()

    val initialKeepSize by option(help="Initial size of the keep map.  Default is 500,000,000")
        .int()
        .default(500_000_000)

    private val refrangeToAdjacentHashCount = mutableMapOf<ReferenceRange, Int>()
    private var runDiagnostics = true

    override fun run() {
        //build the haplotypeGraph
        myLogger.info("Start of BuildKmerIndex...")
        val graph = buildHaplotypeGraph()


        val (hashToHapidMap, discardSet) = processGraphKmers(graph, dbPath, maxHaplotypeProportion,  hashMask, hashFilterValue, initialKeepSize)

        val kmerIndexFilename = if (indexFile == "") "${hvcfDir}/kmerIndex.txt" else indexFile


        //save the kmerIndex
        saveKmerHashesAndHapids(graph, kmerIndexFilename, hashToHapidMap)

        if(discardFile.isNotEmpty()) {
            bufferedWriter(discardFile).use { writer ->
                for (element in discardSet) {
                    writer.write("$element\n")
                }
            }
        }

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
                          hashMask: Long = 3, hashFilterValue:Long = 1, initialKeepSize: Int = 500_000_000) : Pair<Long2ObjectOpenHashMap<Set<String>>,LongOpenHashSet> {
        //keepMap is a map of hash -> Set of haplotype ids
        //TODO Change this to be more RAM efficient.  Probably need something like this:
        // Map<Long, Set<Pair<refRangeIdx, hapIdBitset>>
        val keepMap = Long2ObjectOpenHashMap<MutableSet<String>>(initialKeepSize)

        //discardSet is a Set of hashes
        val discardSet = LongOpenHashSet(initialKeepSize)
        val startTime = System.nanoTime()
        val sampleGametes = graph.sampleGametesInGraph()

        val contigRangesMap = graph.rangesByContig()

        //We set this to be numSampleGametes * .5 so we only keep track of informative kmers
        //We originally tried 2 * numSampleGametes but that was too high as we had too much repetitive mappings
        val maxHapsToKeep = sampleGametes.size * .5

        for (chr in contigRangesMap.keys) {
            //get all sequence for this chromosome
            val agcRequestLists = rangeListsForAgcCommand(graph, contigRangesMap, chr)
            val agcChromSequence = if (agcRequestLists.sampleContigList.isNotEmpty()) {
                myLogger.info("sampleContigList size = ${agcRequestLists.sampleContigList.size}, first element = ${agcRequestLists.sampleContigList[0]}")
                retrieveAgcContigs(dbPath, agcRequestLists.sampleContigList,"")
            }
            else emptyMap()
            val agcOtherRegionSequence: Map<Pair<String,String>, NucSeq> = if (agcRequestLists.otherRegionsList.isNotEmpty()) getAgcSequenceForRanges(agcRequestLists.otherRegionsList)
            else emptyMap()

            //iterate through refranges, generate kmers from the sequence
            for (refrange in contigRangesMap[chr]!!) {
                val hapidToSampleMap = graph.hapIdToSampleGametes(refrange)

                //map haplotype ids to the sequence for that hapid
                val hapidToSequenceMap = extractSequenceForCurrentHapIds(hapidToSampleMap, graph, agcChromSequence, agcOtherRegionSequence)

                //Count Kmer Hashes for the Haplotype Sequences
                countKmerHashesForHaplotypeSequenceSimplified(hapidToSequenceMap, hashMask, hashFilterValue, keepMap, discardSet, maxHapsToKeep)
            }
        }
        myLogger.debug("Finished building kmer keep set, keep set size = ${keepMap.size}, discard set size = ${discardSet.size}, elapse time ${(System.nanoTime() - startTime)/1e9} sec")
        //make sure no hashes in the keepMap are also in discard
        var numberInDiscardSet = 0
        for (element in keepMap.long2ObjectEntrySet()) {
            if (discardSet.contains(element.longKey)) numberInDiscardSet++
        }
        myLogger.info("$numberInDiscardSet kmers in the keepMap are also in the discardSet.")

        //Save the discard Set here as well:


        return Pair(keepMap as Long2ObjectOpenHashMap<Set<String>>, discardSet)
    }

    /**
     * Function to extract out the sequence for a current list of haplotype ids
     */
    private fun extractSequenceForCurrentHapIds(
        hapidToSampleMap: Map<String, List<SampleGamete>>,
        graph: HaplotypeGraph,
        agcChromSequence: Map<Pair<String, String>, NucSeq>,
        agcOtherRegionSequence: Map<Pair<String, String>, NucSeq>
    ): Map<String, List<String>> {

        val hapIdToSequenceMap = hapidToSampleMap.keys.associateWith {
            val altHeader = graph.altHeader(it)!!
            getSeqListForHapId(altHeader, agcChromSequence, agcOtherRegionSequence, it, hapidToSampleMap)
        }

        return hapIdToSequenceMap
    }

    /**
     * Function to pull the needed sequences from AGC
     */
    private fun getSeqListForHapId(
        altHeader: AltHeaderMetaData,
        agcChromSequence: Map<Pair<String, String>, NucSeq>,
        agcOtherRegionSequence: Map<Pair<String, String>, NucSeq>,
        hapid: String,
        hapidToSampleMap: Map<String, List<SampleGamete>>
    ): List<String> {
        val sequenceList = altHeader.regions.map { region ->
            //Regions need to be corrected for any inversions. That is, when requesting sequence
            // make sure that  range.start < range.end.
            // Note that agc positions are 1-based but nucseq positions are 0-based.
            // Because chromosome names can vary between assemblies, for any given assembly
            // the chromosome name in agcChromSequence may be specific to that assembly.
            // To deal with that, first check whether a needed sampleName, contig is in agcChromSequence.
            // If so, use that. If not, get the sequence from afcOtherRegionSequence

            //when requesting from agc, use this range
            val seqRangeStr =
                if (region.first.position <= region.second.position) "${region.first.position}-${region.second.position}"
                else "${region.second.position}-${region.first.position}"
            //when requesting from a NucSeq use this range
            val nucseqRange =
                if (region.first.position <= region.second.position) (region.first.position - 1..<region.second.position)
                else (region.second.position - 1..<region.first.position)

            val contigNuqseq = agcChromSequence[Pair(altHeader.sampleName(), region.first.contig)]

            val regionNucSeq = if (contigNuqseq != null) contigNuqseq[nucseqRange] else {
                agcOtherRegionSequence[Pair(altHeader.sampleName(), "${region.first.contig}:$seqRangeStr")]
            }

            //regionNucSeq should always be found. If it is not, throw an error because something has gone wrong somewhere.
            check(regionNucSeq != null) { "No sequence for ${altHeader.sampleName()} at ${region.first.contig}:$seqRangeStr; hapid $hapid, sample ${hapidToSampleMap[hapid]}" }
            regionNucSeq.seq()
        }
        return sequenceList
    }

    private fun writeDiagnostics(adjacentHashCounts: Map<ReferenceRange, Int>) {
        val diagnosticFileName = "kmerIndexStatistics.txt"
        val diagnosticFilePath = if (indexFile.isBlank() || File(indexFile).parentFile == null) {
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

    private fun getAgcSequenceForRanges(ranges: List<String>): Map<Pair<String,String>, NucSeq> {
        val argLength = ranges.sumOf { it.length } / ranges.size + 1
        val maxArgs = maxArgLength / argLength
        val sequenceMap = mutableMapOf<Pair<String,String>, NucSeq>()
        ranges.windowed(maxArgs, maxArgs, true).forEach {
            sequenceMap.putAll(retrieveAgcContigs(dbPath, it,""))
        }
        return sequenceMap
    }

    /**
     * Function to count the kmerHashes for a single reference range's haplotype nodes.
     * This gets put in 2 sets of maps.
     * One for the hash counts and one for a hash to a list of hapIds which contain that hash.
     * The sequenceList input is a map of hapidId -> list of sequences
     * Note this will likely need to be heavily updated to keep track of relative position in the hapId due to the splitting on N's code
     */
    private fun countKmerHashesForHaplotypeSequenceSimplified(sequenceMap: Map<String, List<String>>, hashMask: Long, hashFilterValue: Long, keepMap:Long2ObjectOpenHashMap<MutableSet<String>>, discardSet: LongOpenHashSet, maxHapsToKeep: Double) {
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
                        //Need to convert these to ULong otherwise the min function will not work correctly
                        //If we leave these as signed Longs, G and T are sorted before A and C due to negative values
                        val minHash = min(previousHash.first.toULong(), previousHash.second.toULong()).toLong()

                        // Use only kmers with a specific ending nucleotide(s)
                        // hashMask determines how many positions will be used
                        // hashFilterValue determines which nucleotide will be kept at that position
                        when {
                            ((minHash and hashMask) != hashFilterValue) -> continue //Doesn't pass hash Mask filter so skip
                            discardSet.contains(minHash) -> continue // Already in discard set so skip
                            keepMap.containsKey(minHash) && keepMap[minHash]!!.size + 1 > maxHapsToKeep -> {
                                //This means that adding a hapID will push us over the maxHapsToKeep parameter so we
                                //should remove from the keepMap and add to the discardSet
                                keepMap.remove(minHash)
                                discardSet.add(minHash)
                            }
                            keepMap.containsKey(minHash) -> keepMap[minHash]!!.add(hapid) // Just add the hapID to the current set for this hash
                            else -> keepMap[minHash] = mutableSetOf(hapid) //Newly seen hash so we need to add a new set in
                        }
                    }
                }
            }
        }
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

            for(refRange in referenceRangeList) {
                if (refRangeToKmerSetMap.containsKey(refRange)) {
                    refRangeToKmerSetMap[refRange]!!.add(kmer)
                } else {
                    refRangeToKmerSetMap[refRange] = mutableSetOf(kmer)
                }
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
                    //agc requests are 1-based positions, so the region positions are correct
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
