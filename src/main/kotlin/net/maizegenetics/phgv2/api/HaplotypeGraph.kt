package net.maizegenetics.phgv2.api

import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import net.maizegenetics.phgv2.utils.AltHeaderMetaData
import net.maizegenetics.phgv2.utils.parseALTHeader
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.*

/**
 * Class to create a HaplotypeGraph from a list of hvcf files.
 * The HaplotypeGraph is a data structure that allows for fast
 * lookup of haplotype sequences for a given ReferenceRange.
 */
class HaplotypeGraph(hvcfFiles: List<String>, dbPath: String? = null) {

    private val myLogger = LogManager.getLogger(HaplotypeGraph::class.java)

    // Map<sampleName, sampleId>
    private lateinit var sampleNameToIdMap: Map<String, Int>
    private lateinit var sampleNames: Array<String>

    // rangeByGameteIdToHapid[refRangeId][sampleId][gameteID] -> checksum / hapid
    // jagged array because different number of haplotypes for each refRange
    private lateinit var rangeByGameteIdToHapid: Array<Array<Array<String>>>

    // Map<ReferenceRange, refRangeId>
    private lateinit var refRangeMap: SortedMap<ReferenceRange, Int>

    // A list of contigs in this graph
    val contigs: List<String>

    // Map<ID (checksum), AltHeaderMetaData>
    private lateinit var altHeaderMap: Map<String, AltHeaderMetaData>

    init {

        getSampleNamesALTHeaders(hvcfFiles)

        if (dbPath == null) {
            processFiles(hvcfFiles)
        } else {
            processFiles(hvcfFiles, dbPath)
        }

        if (rangeByGameteIdToHapid.isNotEmpty()) {
            myLogger.info("rangeToSampleToChecksum: ${rangeByGameteIdToHapid.size} x ${rangeByGameteIdToHapid[0].size}")
            myLogger.info("numOfSamples: ${numberOfSamples()}")
            myLogger.info("numOfRanges: ${numberOfRanges()}")
        }

        contigs = refRangeMap.keys.map { it.contig }.toSortedSet().toList()
    }

    /**
     * Returns the number of samples for this graph.
     */
    fun numberOfSamples() = sampleNameToIdMap.size

    /**
     * Returns a list of samples for this graph.
     */
    fun samples() = sampleNames.toList()

    /**
     * Returns the number of ReferenceRanges for this graph.
     */
    fun numberOfRanges(): Int = refRangeMap.size

    /**
     * Returns a list of ReferenceRanges for this graph.
     */
    fun ranges(): List<ReferenceRange> = refRangeMap.keys.sorted()

    /**
     * Returns a map of contig to a sorted list of ReferenceRanges in each contig
     */
    fun rangesByContig(): Map<String, List<ReferenceRange>> {
        return ranges().groupBy { refrange -> refrange.contig }.mapValues { (_, rangeList) -> rangeList.sorted() }
    }

    /**
     * A sorted set of all SampleGametes in the graph
     */
    fun sampleGametesInGraph(): SortedSet<SampleGamete> {
        val gameteSet = sortedSetOf<SampleGamete>()
        for (refrange in ranges()) {
            gameteSet.addAll(hapIdToSampleGametes(refrange).values.flatten())
        }
        return gameteSet
    }

    fun numSampleGametes(): Int {
        return sampleGametesInGraph().size
    }

    /**
     * Returns a map of hapid -> SampleGamete(s) for the specified ReferenceRange.
     */
    fun hapIdToSampleGametes(range: ReferenceRange): Map<String, List<SampleGamete>> {
        val rangeId = refRangeMap[range]
        require(rangeId != null) { "hapIdToSampleGametes: range: $range not found" }

        val result = mutableMapOf<String, MutableList<SampleGamete>>()
        for (sampleId in rangeByGameteIdToHapid[rangeId].indices) {
            for (gameteId in rangeByGameteIdToHapid[rangeId][sampleId].indices) {
                val hapid = rangeByGameteIdToHapid[rangeId][sampleId][gameteId]
                if (hapid != null && hapid.isNotEmpty()) {
                    result.getOrPut(hapid) { mutableListOf() }.add(SampleGamete(sampleNames[sampleId], gameteId))
                }
            }
        }
        return result
    }

    /**
     * Simple function to make a map of all the haplotypeIds and all the SampleGametes which have that haplotypeId.
     */
    fun hapIdsToSampleGametes(): Map<String, List<SampleGamete>> {
        return ranges().map { hapIdToSampleGametes(it) }
            .flatMap { it.toList() }.toMap()
    }

    /**
     * Creates a map of each SampleGamete in range to its haplotype id.
     * @param range a reference range in this graph
     *
     * If range is not in this graph throws [IllegalArgumentException].
     */
    fun sampleGameteToHaplotypeId(range: ReferenceRange): Map<SampleGamete, String> {
        val rangeId = refRangeMap[range]
        require(rangeId != null) { "hapIdToSamples: range: $range not found" }

        val result = mutableMapOf<SampleGamete, String>()
        for (sampleId in rangeByGameteIdToHapid[rangeId].indices) {
            for (gameteId in rangeByGameteIdToHapid[rangeId][sampleId].indices) {
                val hapid = rangeByGameteIdToHapid[rangeId][sampleId][gameteId]
                if (hapid != "") {
                    result[SampleGamete(sampleNames[sampleId], gameteId)] = hapid
                }
            }
        }
        return result

    }

    /**
     * Returns the hapId for the sample in the specified ReferenceRange.
     */
    fun sampleToHapId(range: ReferenceRange, sample: SampleGamete): String? {
        val rangeId = refRangeMap[range]
        require(rangeId != null) { "sampleToHapId: range: $range not found" }
        val sampleId = sampleNameToIdMap[sample.name]
        require(sampleId != null) { "sampleToHapId: sample: $sample not found" }
        val sampleHapids = rangeByGameteIdToHapid[rangeId][sampleId]
        val result = runCatching { sampleHapids[sample.gameteId] }.getOrNull()
        return if (result == "") {
            null
        } else {
            result
        }
    }

    /**
     * Returns the AltHeaderMetaData for the specified hapId / checksum.
     */
    fun altHeader(hapid: String): AltHeaderMetaData? {
        return altHeaderMap[hapid]
    }

    /**
     * Returns all the AltHeaderMetaData for this graph.
     * Map<ID (checksum), AltHeaderMetaData>
     */
    fun altHeaders() = altHeaderMap

    private fun getSampleNamesALTHeaders(hvcfFiles: List<String>) {

        // Step 1: Get sample names and ALT headers
        val sampleChannel = Channel<Deferred<Pair<List<String>, Map<String, AltHeaderMetaData>>>>(15)
        CoroutineScope(Dispatchers.IO).launch {

            hvcfFiles.forEach { hvcfFile ->

                myLogger.info("processFiles: get samples and ALT headers: $hvcfFile")

                sampleChannel.send(async {
                    VCFFileReader(File(hvcfFile), false).use { reader ->
                        Pair(
                            reader.header.sampleNamesInOrder,
                            parseALTHeader(reader.header)
                        )
                    }
                })

            }

            sampleChannel.close()

        }

        val sampleNamesSet = mutableSetOf<String>()
        val sampleNamesList = mutableListOf<String>()
        val mutableAltHeaderMap: MutableMap<String, AltHeaderMetaData> = mutableMapOf()

        // Single thread to process the sampleChannel
        // This collects the sample names and ALT headers
        // from the HVCF files
        runBlocking {

            for (deferred in sampleChannel) {
                val sampleHeader = deferred.await()

                // sampleNames are being added to both a Set and List so that they can be
                // compared to detect and report duplicate sample names
                sampleNamesSet.addAll(sampleHeader.first)
                sampleNamesList.addAll(sampleHeader.first)

                mutableAltHeaderMap.putAll(sampleHeader.second)
            }

        }

        // make the alt header map immutable
        altHeaderMap = mutableAltHeaderMap.toMap()

        // Step 2: check for duplicate sample names
        myLogger.info("Parsing ALT headers, now checking for duplicate sample names.")
        sampleNamesSet.forEach { sampleNamesList.remove(it) }
        if (sampleNamesList.isNotEmpty()) {
            throw IllegalArgumentException("processFiles: duplicate sample names: $sampleNamesList")
        }

        // Step 3: Initialize sampleNameToIdMap and sampleNames after validation
        myLogger.info("processFiles: initializing sampleNameToIdMap")
        sampleNames = sampleNamesSet.sorted().toTypedArray()
        sampleNameToIdMap = sampleNames.mapIndexed { index, sampleName ->
            Pair(sampleName, index)
        }.toMap()

    }

    private fun processFiles(hvcfFiles: List<String>) {

        myLogger.info("processFiles: ${hvcfFiles.size} hvcf files")

        // rangeIdToSampleToChecksum[rangeId][gameteId][sampleId] -> checksum
        val rangeIdToSampleToChecksum = mutableListOf<MutableList<Array<String?>>>()

        // Step 4: Process the context variants of each HVCF file
        val rangeInfoChannel = Channel<Deferred<List<RangeInfo>>>(5)
        CoroutineScope(Dispatchers.IO).launch {
            hvcfFiles.forEach { hvcfFile ->
                myLogger.info("processFile: $hvcfFile")
                rangeInfoChannel.send(async {
                    VCFFileReader(File(hvcfFile), false).use { reader ->
                        processRanges(reader)
                    }
                })
            }
            rangeInfoChannel.close()
        }

        val rangeMap = mutableMapOf<ReferenceRange, Int>()

        // Single thread to process the rangeInfoChannel
        // This is needed because the rangeMap and rangeIdToSampleToChecksum
        // are shared data structures
        runBlocking {

            for (deferred in rangeInfoChannel) {
                val rangeInfoList = deferred.await()
                for (rangeInfo in rangeInfoList) {
                    val rangeId = rangeMap.getOrPut(rangeInfo.range) { rangeMap.size }

                    if (rangeId < rangeIdToSampleToChecksum.size) {
                        rangeIdToSampleToChecksum[rangeId] =
                            mergeStringArrays(
                                rangeIdToSampleToChecksum.getOrNull(rangeId),
                                rangeInfo.rangeSampleToChecksum
                            )
                    } else {
                        rangeIdToSampleToChecksum.add(
                            rangeId,
                            mergeStringArrays(
                                rangeIdToSampleToChecksum.getOrNull(rangeId),
                                rangeInfo.rangeSampleToChecksum
                            )
                        )
                    }
                }
            }

        }

        refRangeMap = rangeMap.toSortedMap()

        rangeByGameteIdToHapid = convert(rangeIdToSampleToChecksum)

    }

    private fun processFiles(hvcfFiles: List<String>, dbPath: String) {

        myLogger.info("processFiles: ${hvcfFiles.size} hvcf files")

        // rangeIdToSampleToChecksum[rangeId][gameteId][sampleId] -> checksum
        val rangeIdToSampleToChecksum = mutableListOf<MutableList<Array<String?>>>()

        // Step 4: Process the context variants of each HVCF file
        val rangeInfoChannel = Channel<Deferred<List<RangeInfo>>>(5)
        CoroutineScope(Dispatchers.IO).launch {
            hvcfFiles.forEach { hvcfFile ->
                myLogger.info("processFile: $hvcfFile")
                rangeInfoChannel.send(async {
                    VCFFileReader(File(hvcfFile), false).use { reader ->
                        processRanges(reader)
                    }
                })
            }
            rangeInfoChannel.close()
        }

        val rangeMap = mutableMapOf<ReferenceRange, Int>()

        // Single thread to process the rangeInfoChannel
        // This is needed because the rangeMap and rangeIdToSampleToChecksum
        // are shared data structures
        runBlocking {

            for (deferred in rangeInfoChannel) {
                val rangeInfoList = deferred.await()
                for (rangeInfo in rangeInfoList) {
                    val rangeId = rangeMap.getOrPut(rangeInfo.range) { rangeMap.size }

                    if (rangeId < rangeIdToSampleToChecksum.size) {
                        rangeIdToSampleToChecksum[rangeId] =
                            mergeStringArrays(
                                rangeIdToSampleToChecksum.getOrNull(rangeId),
                                rangeInfo.rangeSampleToChecksum
                            )
                    } else {
                        rangeIdToSampleToChecksum.add(
                            rangeId,
                            mergeStringArrays(
                                rangeIdToSampleToChecksum.getOrNull(rangeId),
                                rangeInfo.rangeSampleToChecksum
                            )
                        )
                    }
                }
            }

        }

        refRangeMap = rangeMap.toSortedMap()

        rangeByGameteIdToHapid = convert(rangeIdToSampleToChecksum)

    }

    /**
     * Converts
     * rangeIdToSampleToChecksum[rangeId][gameteId][sampleId] -> checksum
     * to
     * rangeByGameteIdToHapid[refRangeId][sampleId][gameteID] -> checksum
     *
     * First data structure is for faster processing, and
     * the second data structure is for use of the built graph.
     */
    private fun convert(
        rangeIdToSampleToChecksum: MutableList<MutableList<Array<String?>>>
    ): Array<Array<Array<String>>> {

        // Determine the dimensions of the arrays based on the input list
        val rangeSize = rangeIdToSampleToChecksum.size
        val sampleSize = rangeIdToSampleToChecksum[0][0].size
        val gameteSize = rangeIdToSampleToChecksum[0].size

        // Initialize the array for rangeByGameteIdToHapid
        val rangeByGameteIdToHapid = Array(rangeSize) {
            Array(sampleSize) {
                Array(gameteSize) { "" }
            }
        }

        // Loop through the data and reorganize it into the new structure
        for (rangeId in rangeIdToSampleToChecksum.indices) {
            for (gameteId in rangeIdToSampleToChecksum[rangeId].indices) {
                for (sampleId in rangeIdToSampleToChecksum[rangeId][gameteId].indices) {
                    rangeIdToSampleToChecksum[rangeId][gameteId][sampleId]?.let { checksum ->
                        // Move checksum from [rangeId][gameteId][sampleId] -> [rangeId][sampleId][gameteId]
                        rangeByGameteIdToHapid[rangeId][sampleId][gameteId] = checksum
                    }
                }
            }
        }

        return rangeByGameteIdToHapid

    }

    private fun processRanges(reader: VCFFileReader): List<RangeInfo> {
        return reader.map { context ->
            // NOTE: This can throw an exception if there are blank lines in the VCF file
            // The problem seems to occur in the VCFReader processing - checking for blank
            // at this stage does not help.
            contextToRange(context)
        }
    }

    /**
     * ReferenceRange Information
     *
     * @param rangeSampleToChecksum This is the sequence checksums for
     * one Reference Range indexed by gameteId and sampleId
     * @param range This is the ReferenceRange
     */
    data class RangeInfo(
        val rangeSampleToChecksum: MutableList<Array<String?>>,
        val range: ReferenceRange
    )

    /**
     * Convert a VariantContext to the ReferenceRange Information
     */
    private fun contextToRange(
        context: VariantContext
    ): RangeInfo {

        val range = ReferenceRange(context.contig, context.start, context.end)

        // rangeSampleToChecksum[gameteID][sampleID] -> Checksum
        val rangeSampleToChecksum: MutableList<Array<String?>> = mutableListOf()

        context.genotypes.forEach { genotype ->

            val sampleId = sampleNameToIdMap[genotype.sampleName]!!

            genotype.alleles.forEachIndexed { gameteId, allele ->
                var checksumsForSamples = rangeSampleToChecksum.getOrNull(gameteId)
                if (checksumsForSamples == null) {
                    checksumsForSamples = arrayOfNulls(numberOfSamples())
                    rangeSampleToChecksum.add(gameteId, checksumsForSamples)
                }
                checksumsForSamples[sampleId] = allele.displayString.substringAfter("<").substringBefore(">")
            }

        }

        return RangeInfo(rangeSampleToChecksum, range)

    }

    /**
     * Merges two jagged arrays of strings for each sample.
     * This is needed since the input files are processed in
     * multiple threads. This aggregates the results.
     *
     * sampleGameteHapid[gameteId][sampleId] -> Checksum / hapid
     */
    private fun mergeStringArrays(
        gameteSampleHapid1: MutableList<Array<String?>>?,
        gameteSampleHapid2: MutableList<Array<String?>>
    ): MutableList<Array<String?>> {

        return if (gameteSampleHapid1 == null) {
            gameteSampleHapid2
        } else {
            gameteSampleHapid2.forEachIndexed { gameteId, sampleList ->
                var sampleHapid1 = gameteSampleHapid1.getOrNull(gameteId)
                if (sampleHapid1 == null) {
                    sampleHapid1 = arrayOfNulls(sampleList.size)
                    gameteSampleHapid1.add(gameteId, sampleHapid1)
                }
                sampleList.forEachIndexed { sampleId, checksum ->
                    if (sampleHapid1[sampleId] == null) {
                        sampleHapid1[sampleId] = checksum
                    }
                }
            }

            gameteSampleHapid1
        }

    }

    /**
     * Returns a map of hapid -> ReferenceRange
     */
    fun hapIdToRefRangeMap(): Map<String, List<ReferenceRange>> {

        // hapIdToRefRangeMap is a map of hapid -> list of ReferenceRange
        val hapIdToRefRangeMap = mutableMapOf<String, MutableList<ReferenceRange>>()
        for (range in ranges()) {
            for (hapid in hapIdToSampleGametes(range).keys) {
                val refrangeList = hapIdToRefRangeMap.getOrPut(hapid) { mutableListOf() }
                refrangeList.add(range)
            }
        }
        return hapIdToRefRangeMap
    }

    /**
     * Creates a map of ReferenceRange -> (map of hapid -> index)
     */
    fun refRangeToHapIdMap(): Map<ReferenceRange, Map<String, Int>> {
        // This creates a map of ReferenceRangeId -> (map of hapid -> index)
        return ranges().associateWith { range ->
            hapIdToSampleGametes(range).keys.toSortedSet()
                .mapIndexed { index, hapid -> hapid to index }.toMap()
        }
    }

    /**
     * Returns a map of ReferenceRange -> list of all haplotype ids in that range
     */
    fun refRangeToHapIdList(): Map<ReferenceRange, List<String>> {
        return ranges().associateWith { range ->
            hapIdToSampleGametes(range).keys.toList()
        }
    }

    /**
     * Creates a map of ReferenceRangeId -> (map of hapid -> index)
     */
    fun refRangeIdToHapIdMap(): Map<Int, Map<String, Int>> {
        return ranges().mapIndexed { rangeIndex, range ->
            rangeIndex to hapIdToSampleGametes(range).keys.toSortedSet()
                .mapIndexed { hapIndex, hapid -> hapid to hapIndex }.toMap()
        }.toMap()
    }

    /**
     * Creates a map of ReferenceRange -> index for this [HaplotypeGraph].
     * Because the ranges are sorted when calling the ranges() method,
     * a graph always returns the same map in the same order.
     */
    fun refRangeToIndexMap(): Map<ReferenceRange, Int> {
        return ranges().mapIndexed { index, range -> range to index }.toMap()
    }

    /**
     * Creates a map of ReferenceRange Strings -> index for this [HaplotypeGraph].
     */
    fun refRangeStrToIndexMap(): Map<String, Int> {
        return ranges().mapIndexed { index, range -> range.toString() to index }.toMap()
    }

    /**
     * Returns an array of haplotype IDs for every reference range
     * based on a given sample
     */
    fun sampleGameteToHaplotypeId(sample: SampleGamete): Array<String> {
        val result = Array(numberOfRanges()) { "" }

        ranges().forEachIndexed { index, referenceRange ->
            result[index] = sampleToHapId(referenceRange, sample).toString()
        }

        return result
    }

}
