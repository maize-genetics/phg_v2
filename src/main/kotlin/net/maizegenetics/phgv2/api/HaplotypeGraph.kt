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
class HaplotypeGraph(hvcfFiles: List<String>) {

    private val myLogger = LogManager.getLogger(HaplotypeGraph::class.java)

    // Map<sampleName, sampleId>
    private lateinit var sampleNameToIdMap: Map<String, Int>
    private lateinit var sampleNames: Array<String>

    // rangeByGameteIdToHapid[refRangeId][sampleId][gameteID] -> checksum / hapid
    // jagged array because different number of haplotypes for each refRange
    private lateinit var rangeByGameteIdToHapid: Array<Array<Array<String?>>>

    // Map<ReferenceRange, refRangeId>
    private lateinit var refRangeMap: SortedMap<ReferenceRange, Int>

    // Map<ID (checksum), AltHeaderMetaData>
    private lateinit var altHeaderMap: Map<String, AltHeaderMetaData>

    private val processingFiles = Channel<Deferred<Job>>(100)
    private val processingChannel = Channel<RangeInfo>(10)

    init {

        CoroutineScope(Dispatchers.IO).launch {
            processFiles(hvcfFiles)
            closeChannel()
        }

        runBlocking { addSites() }

        if (rangeByGameteIdToHapid.isNotEmpty()) {
            myLogger.info("rangeToSampleToChecksum: ${rangeByGameteIdToHapid.size} x ${rangeByGameteIdToHapid[0].size}")
            myLogger.info("numOfSamples: ${numberOfSamples()}")
            myLogger.info("numOfRanges: ${numberOfRanges()}")
        }

    }

    /**
     * Returns the number of samples for this graph.
     */
    fun numberOfSamples() = sampleNameToIdMap.size

    /**
     * Returns the number of ReferenceRanges for this graph.
     */
    fun numberOfRanges(): Int = refRangeMap.size

    /**
     * Returns a list of ReferenceRanges for this graph.
     */
    fun ranges(): List<ReferenceRange> = refRangeMap.keys.sorted()

    /**
     * Returns a hapId -> sample list map for the given ReferenceRange.
     * Returned Map<hapId, List<sampleName>>
     */
    fun hapIdToSamples(range: ReferenceRange, gameteId: Int = 0): Map<String, List<String>> {
        val rangeId = refRangeMap[range]
        require(rangeId != null) { "hapIdToSamples: range: $range not found" }
        val result = mutableMapOf<String, MutableList<String>>()
        rangeByGameteIdToHapid[rangeId].forEachIndexed { sampleId, hapIdList ->
            result.getOrPut(hapIdList[gameteId]!!) { mutableListOf() }.add(sampleNames[sampleId])
        }
        return result
    }

    /**
     * Returns the hapId for the sample in the specified ReferenceRange.
     */
    fun sampleToHapId(range: ReferenceRange, sample: String, gameteId: Int = 0): String? {
        val rangeId = refRangeMap[range]
        require(rangeId != null) { "sampleToHapId: range: $range not found" }
        val sampleId = sampleNameToIdMap[sample]
        require(sampleId != null) { "sampleToHapId: sample: $sample not found" }
        return rangeByGameteIdToHapid[rangeId][sampleId][gameteId]
    }

    /**
     * Returns the AltHeaderMetaData for the specified hapId / checksum.
     */
    fun altHeader(hapid: String): AltHeaderMetaData? {
        return altHeaderMap[hapid]
    }

    private suspend fun processFiles(hvcfFiles: List<String>) {

        val readers = mutableListOf<VCFFileReader>()
        val sampleNamesSet = mutableSetOf<String>()
        val sampleNamesList = mutableListOf<String>()

        val mutableAltHeaderMap: MutableMap<String, AltHeaderMetaData> = mutableMapOf()
        hvcfFiles.forEach { hvcfFile ->

            myLogger.info("processFiles: $hvcfFile")

            val reader = VCFFileReader(File(hvcfFile), false)
            readers.add(reader)
            sampleNamesSet.addAll(reader.header.sampleNamesInOrder)
            sampleNamesList.addAll(reader.header.sampleNamesInOrder)

            try {
                // extract out the haplotype sequence boundaries for each haplotype from the hvcf
                mutableAltHeaderMap.putAll(parseALTHeader(reader.header))
            } catch (exc: Exception) {
                myLogger.error("processFiles: $hvcfFile: ${exc.message}")
                processingFiles.close()
                processingChannel.close()
                readers.forEach { it.close() }
                throw IllegalArgumentException("processFiles: $hvcfFile: ${exc.message}")
            }
        }
        altHeaderMap = mutableAltHeaderMap.toMap()

        sampleNamesSet.forEach { sampleNamesList.remove(it) }
        if (sampleNamesList.isNotEmpty()) {
            throw IllegalArgumentException("processFiles: duplicate sample names: $sampleNamesList")
        }

        sampleNames = sampleNamesSet.sorted().toTypedArray()
        sampleNameToIdMap = sampleNames.mapIndexed { index, sampleName ->
            Pair(sampleName, index)
        }.toMap()

        readers.forEach { reader ->

            processingFiles.send(CoroutineScope(Dispatchers.IO).async {

                reader.use { reader ->

                    CoroutineScope(Dispatchers.IO).launch {
                        processRanges(reader)
                    }

                }

            })

        }

        processingFiles.close()

    }

    /**
     * Wait for all file processing to complete.
     */
    private suspend fun closeChannel() {
        for (deferred in processingFiles) {
            deferred.await().join()
        }
        processingChannel.close()
    }

    private suspend fun processRanges(reader: VCFFileReader) =
        withContext(Dispatchers.IO) {

            reader.forEach { context ->
                processingChannel.send(contextToRange(context))
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
     * Add reference ranges to data structures, as
     * made available on the processingChannel.
     */
    private suspend fun addSites() {

        // rangeIdToSampleToChecksum[refRangeId][sampleId][gameteId] -> Checksum / hapid
        val rangeIdToSampleToChecksum = mutableListOf<Array<MutableList<String?>>>()

        val rangeMap = mutableMapOf<ReferenceRange, Int>()

        for (rangeInfo in processingChannel) {

            val rangeId = rangeMap.getOrPut(rangeInfo.range) { rangeMap.size }

            if (rangeId < rangeIdToSampleToChecksum.size) {
                rangeIdToSampleToChecksum[rangeId] =
                    mergeStringArrays(rangeIdToSampleToChecksum.getOrNull(rangeId), rangeInfo.rangeSampleToChecksum)
            } else {
                rangeIdToSampleToChecksum.add(
                    rangeId,
                    mergeStringArrays(rangeIdToSampleToChecksum.getOrNull(rangeId), rangeInfo.rangeSampleToChecksum)
                )
            }

        }

        refRangeMap = rangeMap.toSortedMap()

        //rangeByGameteIdToHapid = rangeIdToSampleToChecksum.toTypedArray()
        rangeByGameteIdToHapid = rangeIdToSampleToChecksum.map { sampleGameteHapid ->
            sampleGameteHapid.map { gameteHapid ->
                gameteHapid.toTypedArray()
            }.toTypedArray()
        }.toTypedArray()

    }

    // sampleGameteHapid[sampleId][gameteId] -> Checksum / hapid
    private fun mergeStringArrays(
        sampleGameteHapid1: Array<MutableList<String?>>?,
        sampleGameteHapid2: Array<MutableList<String?>>
    ): Array<MutableList<String?>> {

        return if (sampleGameteHapid1 == null) {
            sampleGameteHapid2
        } else {
            sampleGameteHapid2.forEachIndexed { sampleId, gameteList ->
                gameteList.forEachIndexed { gameteId, gamete ->
                    if (sampleGameteHapid1[sampleId] == null) {
                        sampleGameteHapid1[sampleId] = mutableListOf()
                    }
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

}

