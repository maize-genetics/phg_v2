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

    // seqHash[refRangeId][sampleID] -> checksum / hapid
    // jagged array because different number of haplotypes for each refRange
    private lateinit var rangeToSampleToChecksum: Array<Array<String?>>

    // Map<ReferenceRange, refRangeId>
    private lateinit var refRangeMap: SortedMap<ReferenceRange, Int>

    // Map<ID (checksum), AltHeaderMetaData>
    val altHeaderMap: MutableMap<String, AltHeaderMetaData> = mutableMapOf()

    private val processingFiles = Channel<Deferred<Job>>(100)
    private val processingChannel = Channel<RangeInfo>(10)

    init {

        CoroutineScope(Dispatchers.IO).launch {
            processFiles(hvcfFiles)
            closeChannel()
        }

        runBlocking { addSites() }

        myLogger.info("rangeToSampleToChecksum: ${rangeToSampleToChecksum.size} x ${rangeToSampleToChecksum[0].size}")
        myLogger.info("numOfSamples: ${numberOfSamples()}")
        myLogger.info("numOfRanges: ${numberOfRanges()}")

    }

    /**
     * Returns the number of nodes for this graph.
     */
    fun numberOfNodes(): Int {
        TODO()
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
    fun hapIdToSamples(range: ReferenceRange): Map<String, List<String>> {
        val rangeId = refRangeMap[range]
        require(rangeId != null) { "hapIdToSamples: range: $range not found" }
        val result = mutableMapOf<String, MutableList<String>>()
        rangeToSampleToChecksum[rangeId].forEachIndexed { sampleId, hapId ->
            result.getOrPut(hapId!!) { mutableListOf() }.add(sampleNames[sampleId])
        }
        return result
    }

    /**
     * Returns the hapId for the sample in the specified ReferenceRange.
     */
    fun sampleToHapId(range: ReferenceRange, sample: String): String? {
        val rangeId = refRangeMap[range]
        require(rangeId != null) { "sampleToHapId: range: $range not found" }
        val sampleId = sampleNameToIdMap[sample]
        require(sampleId != null) { "sampleToHapId: sample: $sample not found" }
        return rangeToSampleToChecksum[rangeId][sampleId]
    }

    private suspend fun processFiles(hvcfFiles: List<String>) {

        val readers = mutableListOf<VCFFileReader>()
        val sampleNamesSet = mutableSetOf<String>()
        val sampleNamesList = mutableListOf<String>()

        hvcfFiles.forEach { hvcfFile ->
            val reader = VCFFileReader(File(hvcfFile), false)
            readers.add(reader)
            sampleNamesSet.addAll(reader.header.sampleNamesInOrder)
            sampleNamesList.addAll(reader.header.sampleNamesInOrder)

            // extract out the haplotype sequence boundaries for each haplotype from the hvcf
            altHeaderMap.putAll(parseALTHeader(reader.header))
        }

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
     * one Reference Range indexed by sampleId
     * @param range This is the ReferenceRange
     */
    data class RangeInfo(
        val rangeSampleToChecksum: Array<String?>,
        val range: ReferenceRange
    )

    /**
     * Convert a VariantContext to the ReferenceRange Information
     */
    private fun contextToRange(
        context: VariantContext
    ): RangeInfo {

        val range = ReferenceRange(context.contig, context.start, context.end)

        val rangeSampleToChecksum = Array<String?>(numberOfSamples()) { null }
        context.genotypes.forEach { genotype ->

            val sampleId = sampleNameToIdMap[genotype.sampleName]!!

            // Number of alleles should always be 1 for a haplotype
            // Diploids need to be handled by sample names
            check(genotype.alleles.size == 1) { "genotype.alleles.size != 1" }

            rangeSampleToChecksum[sampleId] = genotype.alleles[0].displayString.substringAfter("<").substringBefore(">")
        }

        return RangeInfo(rangeSampleToChecksum, range)

    }

    /**
     * Add reference ranges to data structures, as
     * made available on the processingChannel.
     */
    private suspend fun addSites() {

        // rangeIdToSampleToChecksum[refRangeId][sampleId] -> Checksum / hapid
        val rangeIdToSampleToChecksum = mutableListOf<Array<String?>>()

        val rangeMap = mutableMapOf<ReferenceRange, Int>()

        for (rangeInfo in processingChannel) {

            val rangeId = rangeMap.getOrPut(rangeInfo.range) { rangeMap.size }
            myLogger.info("processingChannel rangeId: $rangeId")

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

        rangeToSampleToChecksum = rangeIdToSampleToChecksum.toTypedArray()

    }

    private fun mergeStringArrays(seqHash1: Array<String?>?, seqHash2: Array<String?>): Array<String?> {
        return seqHash1?.let {
            Array(seqHash2.size) { index ->
                if (seqHash1[index] != null) seqHash1[index] else seqHash2[index]
            }
        } ?: seqHash2
    }

}

