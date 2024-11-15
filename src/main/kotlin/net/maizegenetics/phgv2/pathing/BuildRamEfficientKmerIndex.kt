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
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.pathing.AlignmentUtils.Companion.buildHaplotypeGraph
import net.maizegenetics.phgv2.utils.AltHeaderMetaData
import net.maizegenetics.phgv2.utils.retrieveAgcContigs
import org.apache.logging.log4j.LogManager
import java.io.File
import kotlin.math.min
import kotlin.time.DurationUnit
import kotlin.time.measureTimedValue

data class IndexedKmerData(val chr: String, val pos: UInt, var gameteLong0 : ULong, var gameteLong1 : ULong)
class BuildRamEfficientKmerIndex: CliktCommand(help="Create a kmer index for a HaplotypeGraph. By default the file will be written " +
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
        myLogger.info("Start of BuildRamEfficientKmerIndex...")
        val graph = buildHaplotypeGraph(hvcfDir)

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
        else {
            writeDiagnostics(refrangeToAdjacentHashCount)
        }

    }


    // This is derived from the function processGraphKmers in BuildKmerIndex.kt
    // This implementation uses a mutableMap, which should be able to
    // grow as needed.  The bitmaps that contain the gamete information will take up
    // less room and need less storage.
    fun processGraphKmers(graph: HaplotypeGraph, dbPath: String, maxHaplotypeProportion: Double=.75,
                          hashMask: Long = 3, hashFilterValue:Long = 1, initialKeepSize: Int = 500_000_000) : Pair<Map<Long,IndexedKmerData>, LongOpenHashSet> {

        // keepMap is a map of hapidHash -> IndexedKmerData
        val keepMap = mutableMapOf<Long,IndexedKmerData>()
        //discardSet is a Set of hashes
        val discardSet = LongOpenHashSet(initialKeepSize)
        val startTime = System.nanoTime()

        // This returns a SortedSEt of <SampleGamete> - should it be turned into an ordered
        // list of just sample names?  This list is needed to index into the bit arrays
        val sampleGametes = graph.sampleGametesInGraph()
        val samples = sampleGametes.map { it.name }.sorted()

        val contigRangesMap = graph.rangesByContig()

        //We set this to be numSampleGametes * .5 so we only keep track of informative kmers
        //We originally tried 2 * numSampleGametes but that was too high as we had too much repetitive mappings
        val maxHapsToKeep = sampleGametes.size * .5

        for (chr in contigRangesMap.keys) {
            //get all sequence for this chromosome
            val agcRequestLists = BuildKmerIndex.rangeListsForAgcCommand(graph, contigRangesMap, chr)
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

                //  TODO - implement this function - it will have most of our changes based on keepMap changes
                // We  will need to add the hapidToSampleMap as a parameter - this will be needed to get the sample
                // for the hapid to index into the bit array
                //Count Kmer Hashes for the Haplotype Sequences
                countKmerHashesForHaplotypeSequenceSimplified(hapidToSequenceMap, hashMask, hashFilterValue, keepMap,
                    discardSet, maxHapsToKeep,hapidToSampleMap,samples,refrange)
            }
        }
        myLogger.debug("Finished building kmer keep set, keep set size = ${keepMap.size}, discard set size = ${discardSet.size}, elapse time ${(System.nanoTime() - startTime)/1e9} sec")

        //make sure no hashes in the keepMap are also in discard
        var numberInDiscardSet = 0
        for (element in keepMap.keys) {
            if (discardSet.contains(element)) numberInDiscardSet++
        }
        myLogger.info("$numberInDiscardSet kmers in the keepMap are also in the discardSet.")

        return Pair(keepMap, discardSet)

    }

    /**
     * Function to count the kmerHashes for a single reference range's haplotype nodes.
     * This gets put in 2 sets of maps.
     * One for the hash counts and one for a hash to a list of hapIds which contain that hash.
     * The sequenceList input is a map of hapidId -> list of sequences
     * Note this will likely need to be heavily updated to keep track of relative position in the hapId due to the splitting on N's code
     * "samples" is an ordered list of samples from the graph
     * "hapidToSampleMap" maps the hapid to a list of samples which contain that hapid.
     */
    private fun countKmerHashesForHaplotypeSequenceSimplified(sequenceMap: Map<String, List<String>>, hashMask: Long, hashFilterValue: Long,
                                                              keepMap:MutableMap<Long,IndexedKmerData>, discardSet: LongOpenHashSet, maxHapsToKeep: Double,
                                                              hapidToSampleMap: Map<String, List<SampleGamete>>, samplesList:List<String>, refRange: ReferenceRange) {
        for ((hapid, seqList) in sequenceMap) {
            val sample = hapidToSampleMap[hapid]!!.first().name
            for (sequenceWithNs in seqList) {
                //split sequence on N's then filter on length > 31 because we are looking for 32-mers
                val splitList = sequenceWithNs
                    .split("N+".toRegex())
                    .filter{it.length > 31}
                for (sequence in splitList) {
                    var previousHash = Pair(0L, 0L)

                    //for first 31 nucleotides just update the hash
                    for (nucleotide in sequence.subSequence(0..30)) {
                        previousHash = BuildKmerIndex.updateKmerHashAndReverseCompliment(previousHash, nucleotide)
                    }

                    //start using kmers starting with the 32nd nucleotide
                    for (nucleotide in sequence.subSequence(31 until sequence.length)) {
                        previousHash = BuildKmerIndex.updateKmerHashAndReverseCompliment(previousHash, nucleotide)
                        //Need to convert these to ULong otherwise the min function will not work correctly
                        //If we leave these as signed Longs, G and T are sorted before A and C due to negative values
                        val minHash = min(previousHash.first.toULong(), previousHash.second.toULong()).toLong()

                        val keyInMap = keepMap.containsKey(minHash)
                        var bitCount1 = 0
                        var bitCount2 = 0
                        if (keyInMap) {
                            bitCount1 = keepMap[minHash]!!.gameteLong0.countOneBits()
                            bitCount2 = keepMap[minHash]!!.gameteLong1.countOneBits()
                        }
                        // Use only kmers with a specific ending nucleotide(s)
                        // hashMask determines how many positions will be used
                        // hashFilterValue determines which nucleotide will be kept at that position

                        when {
                            ((minHash and hashMask) != hashFilterValue) -> continue //Doesn't pass hash Mask filter so skip
                            discardSet.contains(minHash) -> continue // Already in discard set so skip
                            keyInMap && (bitCount1 + bitCount2 + 1 > maxHapsToKeep) -> {
                                //This means that adding a hapID will push us over the maxHapsToKeep parameter so we
                                //should remove from the keepMap and add to the discardSet
                                keepMap.remove(minHash)
                                discardSet.add(minHash)
                            }
                            keyInMap -> {
                                var kmerDetails = keepMap[minHash]!!

                                val index = samplesList.indexOf(sample)
                                if (index == -1) {
                                    throw IllegalArgumentException("Sample ${sample} not found in the list.")
                                }
                                when (index) {
                                    in 0..63 -> {
                                        kmerDetails.gameteLong0 = kmerDetails.gameteLong0 or (1uL shl index)
                                        keepMap[minHash] = kmerDetails
                                    } // Update gamete0
                                    in 64..127 -> {
                                        kmerDetails.gameteLong1 = kmerDetails.gameteLong1 or (1uL shl (index - 64))
                                        keepMap[minHash] = kmerDetails
                                    } // Update gamete1
                                    else -> {
                                        // no changes?  Are we ignoring when have more than 128 samples?
                                    }
                                }

                            } // New Entry: add the hapID to the current set for this hash
                            else  -> { // create new entry for map
                                var gameteLong0 = 0uL
                                var gameteLong1 = 0uL
                                val index = samplesList.indexOf(sample)
                                if (index == -1) {
                                    throw IllegalArgumentException("Sample ${sample} not found in the list.")
                                }
                                when (index) {
                                    in 0..63 -> {
                                        gameteLong0 = gameteLong0 or (1uL shl index)
                                    } // Update gamete0
                                    in 64..127 -> {
                                        gameteLong1 = gameteLong1 or (1uL shl (index - 64))
                                    } // Update gamete1
                                    else -> {
                                        // no changes?  Are we ignoring when have more than 128 samples?
                                    }
                                }
                                keepMap[minHash] = IndexedKmerData(refRange.contig, refRange.start.toUInt(),gameteLong0,gameteLong1)
                            } //Newly seen hash so we need to add a new object
                        }
                    }
                }
            }
        }
    }
    private fun saveKmerHashesAndHapids(
        graph: HaplotypeGraph,
        kmerIndexFilename: String,
        hashToHapidMap: Map<Long, IndexedKmerData>
    ) {
        TODO("Not yet implemented")
    }

    private fun writeDiagnostics(refrangeToAdjacentHashCount: MutableMap<ReferenceRange, Int>) {
        TODO("Not yet implemented")

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

}