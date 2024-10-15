package net.maizegenetics.phgv2.api

import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.utils.AltHeaderMetaData
import net.maizegenetics.phgv2.utils.parseALTHeader
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.*

/**
 * Class to create a HaplotypeGraph from a list of hvcf files.
 * The HaplotypeGraph is a data structure that allows for fast
 * lookup of haplotype sequences for a given ReferenceRange.
 *
 * Current version is re-worked to removed co-routines and use
 * single-threaded.  The computation for convering variant context
 * to refRange is very quick.  Co-routine overhead caused major
 * delays in processing - this is much faster.
 */
class HaplotypeGraph(hvcfFiles: List<String>) {

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
        processFiles(hvcfFiles)

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
            .reduce { acc, map -> acc + map }
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
                result[SampleGamete(sampleNames[sampleId], gameteId)] =
                    rangeByGameteIdToHapid[rangeId][sampleId][gameteId]
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
        return if (sampleHapids.size > sample.gameteId) sampleHapids[sample.gameteId] else null
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

    fun processFiles(hvcfFiles: List<String>) {

        val rangeIdToSampleToChecksum = mutableListOf<Array<MutableList<String?>>>()

        val rangeMap = mutableMapOf<ReferenceRange, Int>()

        val sampleNamesSet = mutableSetOf<String>()
        val sampleNamesList = mutableListOf<String>()
        val mutableAltHeaderMap: MutableMap<String, AltHeaderMetaData> = mutableMapOf()

        myLogger.info("processFiles: ${hvcfFiles.size} hvcf files")

        // Step 1: Get sample names and parse ALT headers
        hvcfFiles.forEach { hvcfFile ->

            myLogger.info("processFiles: $hvcfFile")

            VCFFileReader(File(hvcfFile), false).use { reader ->
                // sampleNames are being added to both a Set and List so that they can be compared to detect and report
                // duplicate sample names
                sampleNamesSet.addAll(reader.header.sampleNamesInOrder)
                sampleNamesList.addAll(reader.header.sampleNamesInOrder)

                parseALTHeader(reader.header, mutableAltHeaderMap)

                reader.close()
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

        // Initialize sampleNameToIdMap and sampleNames after validation
        myLogger.info("processFiles: initializing sampleNameToIdMap")
        sampleNames = sampleNamesSet.sorted().toTypedArray()
        sampleNameToIdMap = sampleNames.mapIndexed { index, sampleName ->
            Pair(sampleName, index)
        }.toMap()

        // Step 3: Process the files sequentially
        hvcfFiles.forEach { hvcfFile ->

            myLogger.info("processFiles: $hvcfFile")

            VCFFileReader(File(hvcfFile), false).use { reader ->

                val rangeInfoList = processRanges(reader)
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

                reader.close()

            }

        }

        myLogger.info("Finished with readers loop to process ranges.")

        refRangeMap = rangeMap.toSortedMap()

        rangeByGameteIdToHapid = rangeIdToSampleToChecksum.map { sampleGameteHapid ->
            sampleGameteHapid.map { gameteHapid ->
                gameteHapid.map { it ?: "" }.toTypedArray()
            }.toTypedArray()
        }.toTypedArray()

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
     * one Reference Range indexed by sampleId and gameteId
     * @param range This is the ReferenceRange
     */
    data class RangeInfo(
        val rangeSampleToChecksum: Array<MutableList<String?>>,
        val range: ReferenceRange
    )

    /**
     * Convert a VariantContext to the ReferenceRange Information
     */
    private fun contextToRange(
        context: VariantContext
    ): RangeInfo {

        val range = ReferenceRange(context.contig, context.start, context.end)

        val rangeSampleToChecksum = Array<MutableList<String?>>(numberOfSamples()) { mutableListOf() }
        context.genotypes.forEach { genotype ->

            val sampleId = sampleNameToIdMap[genotype.sampleName]!!

            genotype.alleles.forEach { allele ->
                rangeSampleToChecksum[sampleId].add(allele.displayString.substringAfter("<").substringBefore(">"))
            }

        }

        return RangeInfo(rangeSampleToChecksum, range)

    }

    /**
     * Merges two jagged arrays of strings for each sample.
     * This is needed since the input files are processed in
     * multiple threads. This aggregates the results.
     *
     * sampleGameteHapid[sampleId][gameteId] -> Checksum / hapid
     */
    private fun mergeStringArrays(
        sampleGameteHapid1: Array<MutableList<String?>>?,
        sampleGameteHapid2: Array<MutableList<String?>>
    ): Array<MutableList<String?>> {

        return if (sampleGameteHapid1 == null) {
            sampleGameteHapid2
        } else {
            sampleGameteHapid2.forEachIndexed { sampleId, gameteList ->
                gameteList.forEachIndexed { gameteId, gamete ->
                    if (gameteId < sampleGameteHapid1[sampleId].size) {
                        if (sampleGameteHapid1[sampleId][gameteId] == null) {
                            sampleGameteHapid1[sampleId][gameteId] = gamete
                        }
                    } else {
                        sampleGameteHapid1[sampleId].add(gameteId, gamete)
                    }
                }
            }
            sampleGameteHapid1
        }

    }


    /**
     * Returns a map of hapid -> ReferenceRange
     */
    fun hapIdToRefRangeMap(): Map<String, List<ReferenceRange>> {

        //hapIdToRefRangeMap is a map of hapid -> list of ReferenceRange
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
        //This creates a map of ReferenceRangeId -> (map of hapid -> index)
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

