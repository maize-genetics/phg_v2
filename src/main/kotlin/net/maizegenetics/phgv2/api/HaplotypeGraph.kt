package net.maizegenetics.phgv2.api

import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import net.maizegenetics.phgv2.utils.AltHeaderMetaData
import net.maizegenetics.phgv2.utils.parseALTHeader
import java.io.File
import java.util.*

class HaplotypeGraph(hvcfFile: String) {

    // Map<sampleName, sampleId>
    private val sampleNameToIdMap: Map<String, Int>

    private val numOfSamples: Int

    // lookup[refRangeId][ploidy][sampleId]
    private lateinit var lookup: Array<Array<Array<UByte>>>

    // seqHash[refRangeId][lookup: UByte]
    // jagged array because different number of haplotypes for each refRange
    private lateinit var seqHash: Array<Array<String>>

    private lateinit var refRangeMap: SortedMap<Int, ReferenceRange>

    private val altHeaderMap: Map<String, AltHeaderMetaData>

    fun numOfRanges(): Int = refRangeMap.size

    private val processingChannel = Channel<RangeInfo>(5)

    init {

        VCFFileReader(File(hvcfFile), false).use { reader ->

            // create a map of sampleName to sampleId
            sampleNameToIdMap = reader.header.sampleNamesInOrder.mapIndexed { index, sampleName ->
                Pair(sampleName, index)
            }.toMap()

            numOfSamples = sampleNameToIdMap.size

            // extract out the haplotype sequence boundaries for each haplotype from the hvcf
            altHeaderMap = parseALTHeader(reader.header)

            CoroutineScope(Dispatchers.IO).launch {
                processRanges(reader)
            }

            runBlocking { addSites() }

        }

    }

    private suspend fun processRanges(reader: VCFFileReader) =
        withContext(Dispatchers.IO) {

            reader.forEachIndexed { index, context ->
                processingChannel.send(contextToRange(context, index))
            }

            processingChannel.close()

        }

    data class RangeInfo(
        val rangeLookup: Array<Array<UByte>>,
        val rangeSeqHash: Array<String>,
        val rangeId: Int,
        val range: ReferenceRange
    )

    private fun contextToRange(
        context: VariantContext,
        rangeId: Int
    ): RangeInfo {

        val range = ReferenceRange(rangeId, context.contig, context.start, context.end)

        val ploidy = context.getMaxPloidy(2)

        val symToID = context.alleles.mapIndexed { index, allele ->
            val symbolicAllele = allele.displayString.substringAfter("<").substringBefore(">")
            Pair(symbolicAllele, index)
        }.toMap()

        val rangeSeqHash = context.alleles.map { allele ->
            allele.displayString.substringAfter("<").substringBefore(">")
        }.toTypedArray()

        val rangeLookup = Array(ploidy) { Array(numOfSamples) { UByte.MAX_VALUE } }
        context.genotypes.forEach { genotype ->
            val sampleId = sampleNameToIdMap[genotype.sampleName]!!
            genotype.alleles.forEachIndexed { index, allele ->
                val alleleId = symToID[allele.displayString.substringAfter("<").substringBefore(">")]!!
                rangeLookup[index][sampleId] = alleleId.toUByte()
            }
        }

        return RangeInfo(rangeLookup, rangeSeqHash, rangeId, range)

    }

    private suspend fun addSites() {

        val lookupList = mutableListOf<Array<Array<UByte>>>()
        val seqHashList = mutableListOf<Array<String>>()
        val rangeMap = mutableMapOf<Int, ReferenceRange>()

        for (rangeInfo in processingChannel) {
            lookupList.add(rangeInfo.rangeLookup)
            seqHashList.add(rangeInfo.rangeSeqHash)
            rangeMap[rangeInfo.rangeId] = rangeInfo.range
        }

        lookup = lookupList.toTypedArray()
        seqHash = seqHashList.toTypedArray()
        refRangeMap = rangeMap.toSortedMap()

    }

}

