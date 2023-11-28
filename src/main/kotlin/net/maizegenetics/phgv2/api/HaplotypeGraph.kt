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

    // lookup the sequence checksum index
    // which is the index into seqHash
    // lookup[refRangeId][sampleId]
    private lateinit var lookup: Array<Array<UByte>>

    // seqHash[refRangeId][lookup: UByte]
    // jagged array because different number of haplotypes for each refRange
    private lateinit var seqHash: Array<Array<String>>

    // Map<ReferenceRange, refRangeId>
    private lateinit var refRangeMap: SortedMap<ReferenceRange, Int>

    // Map<ID (checksum), AltHeaderMetaData>
    private val altHeaderMap: MutableMap<String, AltHeaderMetaData> = mutableMapOf()

    private val processingFiles = Channel<Deferred<Job>>(100)
    private val processingChannel = Channel<RangeInfo>(10)

    init {

        CoroutineScope(Dispatchers.IO).launch {
            processFiles(hvcfFiles)
            closeChannel()
        }

        runBlocking { addSites() }

        myLogger.info("lookup: ${lookup.size} x ${lookup[0].size}")
        myLogger.info("seqHash: ${seqHash.size} x ${seqHash[0].size}")
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
        lookup[rangeId].forEachIndexed { index, seqIndex ->
            val hapId = seqHash[rangeId][seqIndex.toInt()]
            val sampleName = sampleNames[index]
            result.getOrPut(hapId) { mutableListOf() }.add(sampleName)
        }
        return result
    }

    /**
     * Returns the hapId for the sample in the specified ReferenceRange.
     */
    fun sampleToHapId(range: ReferenceRange, sample: String): String {
        val rangeId = refRangeMap[range]
            ?: throw IllegalArgumentException("sampleToHapId: range: $range not found")
        val sampleId = sampleNameToIdMap[sample]
            ?: throw IllegalArgumentException("sampleToHapId: sample: $sample not found")
        return seqHash[rangeId][lookup[rangeId][sampleId].toInt()]
    }

    private suspend fun processFiles(hvcfFiles: List<String>) {

        val readers = mutableListOf<VCFFileReader>()
        val sampleNamesSet = mutableSetOf<String>()

        hvcfFiles.forEach { hvcfFile ->
            val reader = VCFFileReader(File(hvcfFile), false)
            readers.add(reader)
            sampleNamesSet.addAll(reader.header.sampleNamesInOrder)

            // extract out the haplotype sequence boundaries for each haplotype from the hvcf
            altHeaderMap.putAll(parseALTHeader(reader.header))
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
     * @param rangeLookup This is the lookup sequence checksum
     * indices for one Reference Range rangeLookup[sampleId]
     * @param rangeSeqHash This is the sequence checksums for
     * one Reference Range rangeSeqHash[lookup: UByte]
     * @param range This is the ReferenceRange
     */
    data class RangeInfo(
        val rangeLookup: Array<UByte>,
        val rangeSeqHash: Array<String>,
        val range: ReferenceRange
    )

    /**
     * Convert a VariantContext to the ReferenceRange Information
     */
    private fun contextToRange(
        context: VariantContext
    ): RangeInfo {

        val range = ReferenceRange(context.contig, context.start, context.end)

        val symToID = context.alleles.mapIndexed { index, allele ->
            val symbolicAllele = allele.displayString.substringAfter("<").substringBefore(">")
            Pair(symbolicAllele, index)
        }.toMap()

        val rangeSeqHash = context.alleles.map { allele ->
            allele.displayString.substringAfter("<").substringBefore(">")
        }.toTypedArray()

        val rangeLookup = Array(numberOfSamples()) { UByte.MAX_VALUE }
        context.genotypes.forEach { genotype ->
            val sampleId = sampleNameToIdMap[genotype.sampleName]!!
            check(genotype.alleles.size == 1) { "genotype.alleles.size != 1" }
            val alleleId = symToID[genotype.alleles[0].displayString.substringAfter("<").substringBefore(">")]!!
            rangeLookup[sampleId] = alleleId.toUByte()
        }

        return RangeInfo(rangeLookup, rangeSeqHash, range)

    }

    /**
     * Add reference ranges to data structures, as
     * made available on the processingChannel.
     */
    private suspend fun addSites() {

        val lookupList = mutableListOf<Array<UByte>>()
        val seqHashList = mutableListOf<Array<String>>()
        val rangeMap = mutableMapOf<ReferenceRange, Int>()

        for (rangeInfo in processingChannel) {
            val rangeId = rangeMap.getOrPut(rangeInfo.range) { rangeMap.size }
            lookupList.add(rangeId, rangeInfo.rangeLookup)
            seqHashList.add(rangeId, rangeInfo.rangeSeqHash)
        }

        lookup = lookupList.toTypedArray()
        seqHash = seqHashList.toTypedArray()
        refRangeMap = rangeMap.toSortedMap()

    }

}

