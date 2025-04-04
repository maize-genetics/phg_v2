package net.maizegenetics.phgv2.api

import biokotlin.seq.NucSeqRecord
import biokotlin.seqIO.NucSeqIO
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFHeader
import htsjdk.variant.vcf.VCFHeaderLine
import net.maizegenetics.phgv2.utils.AltHeaderMetaData
import net.maizegenetics.phgv2.utils.altHeaderMetadataToVCFHeaderLine
import net.maizegenetics.phgv2.utils.createGenericHeaderLineSet
import net.maizegenetics.phgv2.utils.loadRanges
import org.apache.logging.log4j.LogManager
import java.io.File

private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.api.ExportHaplotypeGraph")

enum class SymbolicAllele {
    CHECKSUM, RANGE_SAMPLE_GAMETE
}

/**
 * Export a HaplotypeGraph to a multi-sample h.vcf file.
 * The h.vcf is the format designed by the PHGv2.
 * The alternate alleles are symbolic alleles that are check sums of the haplotype sequences.
 *
 * The rangeBedfile is a bedfile that contains the master list of reference ranges.
 * If no rangeBedfile is supplied, it will use the reference ranges from the graph.
 *
 * @param graph The HaplotypeGraph to export.
 * @param filename The name of the file to export to.
 * @param referenceGenome The filename of the reference genome to use.
 * If null, the reference genome will not be used.
 * @param symbolicAllele The type of symbolic allele to use.
 * @param rangeBedfile The filename of the bedfile to use for the ranges.
 */
fun exportMultiSampleHVCF(
    graph: HaplotypeGraph,
    filename: String,
    referenceGenome: String? = null,
    symbolicAllele: SymbolicAllele = SymbolicAllele.CHECKSUM,
    rangeBedfile: String? = null
) {

    // Load the reference genome into memory if filename is supplied.
    val referenceSequence = referenceGenome?.let {
        NucSeqIO(it).readAll()
    }

    VariantContextWriterBuilder()
        .unsetOption(Options.INDEX_ON_THE_FLY)
        .setOutputFile(File(filename))
        .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
        .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
        .build()
        .use { writer ->

            // Map of <rangeStr: String, index: Int>
            // chr1:1-1000 -> 0
            val rangeStrToIndex = if (rangeBedfile != null) {
                loadRanges(rangeBedfile)
                    .mapIndexed { index, range ->
                        "${range.first.contig}:${range.first.position}-${range.second.position}" to index
                    }.toMap()
            } else {
                graph.refRangeStrToIndexMap()
            }

            val altHeaders = graph.altHeaders()

            // fun altHeaderMetadataToVCFHeaderLine(altHeaderData: AltHeaderMetaData, altHeaderId: String? = null): VCFAltHeaderLine {
            // map of hapid to VCFAltHeaderLine
            val hapidToHeaderLine = mutableMapOf<String, Pair<AltHeaderMetaData, VCFHeaderLine>>()

            graph.ranges()
                .forEach { range ->
                    val rangeStr = range.toString()
                    rangeStrToIndex[rangeStr]?.let { index ->
                        graph.hapIds(range)
                            .forEach { hapid ->
                                hapidToHeaderLine.computeIfAbsent(hapid) {
                                    val metaData = createAltHeaderMetaData(
                                        altHeaders[hapid]!!,
                                        rangeStr,
                                        index,
                                        symbolicAllele
                                    )
                                    Pair(metaData, altHeaderMetadataToVCFHeaderLine(metaData))
                                }
                            }
                    }
                }

            val headerLines = createGenericHeaderLineSet().toMutableSet()
            headerLines.addAll(hapidToHeaderLine.values.map { it.second })
            headerLines.addAll(graph.contigHeaderLines())
            val header = VCFHeader(headerLines, graph.samples())

            writer.writeHeader(header)

            graph.ranges()
                .forEach { range ->
                    val rangeStr = range.toString()
                    rangeStrToIndex[rangeStr]?.let { index ->
                        val variantContext =
                            createVariantContext(
                                range,
                                index,
                                graph.hapIdToSampleGametes(range),
                                referenceSequence,
                                symbolicAllele,
                                hapidToHeaderLine
                            )
                        writer.add(variantContext)
                    }
                }

        }

}

/**
 * Checks if the input is a valid MD5 checksum.
 */
private fun isMd5Checksum(input: String): Boolean {
    // Regex to match a valid MD5 checksum (32 hexadecimal characters)
    val md5Regex = Regex("^[a-fA-F0-9]{32}$")
    return md5Regex.matches(input)
}

private fun createAltHeaderMetaData(
    origAltHeader: AltHeaderMetaData,
    rangeStr: String,
    rangeStrIndex: Int,
    symbolicAllele: SymbolicAllele
): AltHeaderMetaData {

    val sample = origAltHeader.sampleGamete

    // Find what is the checksum
    val checksum = when {
        isMd5Checksum(origAltHeader.checksum) -> origAltHeader.checksum
        isMd5Checksum(origAltHeader.id) -> origAltHeader.id
        else -> throw IllegalArgumentException("Invalid checksum: ${origAltHeader.id} ${origAltHeader.checksum}")
    }

    val refChecksum = when {
        origAltHeader.refChecksum.isNotBlank() -> origAltHeader.refChecksum
        isMd5Checksum(origAltHeader.refRange) -> origAltHeader.refRange
        else -> ""
    }

    val id = when (symbolicAllele) {

        SymbolicAllele.RANGE_SAMPLE_GAMETE -> {
            symbolicAlleleRangeSampleGameteStr(
                rangeStrIndex,
                sample.name,
                sample.gameteId
            )
        }

        SymbolicAllele.CHECKSUM -> checksum

    }

    return AltHeaderMetaData(
        id,
        origAltHeader.description,
        origAltHeader.source,
        sample,
        origAltHeader.regions,
        checksum,
        rangeStr,
        refChecksum
    )

}

private fun createVariantContext(
    range: ReferenceRange,
    rangeIndex: Int,
    hapIdToSampleGametes: Map<String, List<SampleGamete>>,
    referenceSequence: Map<String, NucSeqRecord>?,
    symbolicAllele: SymbolicAllele,
    hapidToHeaderLine: Map<String, Pair<AltHeaderMetaData, VCFHeaderLine>>
): VariantContext {

    // alleles: Map<hapid: String, Allele>
    val alleles = hapIdToSampleGametes.keys
        .asSequence()
        .map { hapid ->

            when (symbolicAllele) {
                SymbolicAllele.CHECKSUM -> Pair(hapid, symbolicAlleleAlt(hapidToHeaderLine[hapid]!!.first.checksum))
                SymbolicAllele.RANGE_SAMPLE_GAMETE -> {
                    val sample = hapIdToSampleGametes[hapid]!!.firstNotNullOf { it }
                    Pair(hapid, symbolicAlleleRangeSampleGamete(rangeIndex, sample.name, sample.gameteId))
                }
            }

        }
        .toMap()

    val taxaToHapids = mutableMapOf<String, MutableList<String>>()
        .apply {
            hapIdToSampleGametes.forEach { (hapid, gametes) ->
                gametes
                    .sorted()
                    .forEach { gamete ->
                        val hapids = getOrPut(gamete.name) { mutableListOf() }
                        hapids.add(hapid)
                    }
            }
        }

    val genotypes = taxaToHapids
        .map { (taxon, hapids) ->
            if (hapids.isEmpty()) {
                GenotypeBuilder(taxon, listOf(Allele.NO_CALL)).phased(false).make()
            } else {
                val allelesForHapids = hapids.map { alleles[it] }
                GenotypeBuilder(taxon, allelesForHapids).phased(true).make()
            }
        }
        .toList()

    val resultAlleles = mutableListOf<Allele>()
    resultAlleles.add(alleleRef(range, referenceSequence))
    alleles.values
        .filter { it != Allele.NO_CALL }
        .forEach { resultAlleles.add(it) }

    return VariantContextBuilder()
        .source(".")
        .alleles(resultAlleles)
        .chr(range.contig)
        .start(range.start.toLong())
        .stop(range.end.toLong())
        .attribute("END", range.end.toString())
        .genotypes(genotypes)
        .make()

}

private fun alleleRef(range: ReferenceRange, referenceSequence: Map<String, NucSeqRecord>?): Allele {
    val allele = if (referenceSequence != null) {
        referenceSequence[range.contig]!!.sequence[range.start - 1].name
    } else {
        "A"
    }
    return Allele.create(allele, true)
}

private fun symbolicAlleleAlt(hapid: String): Allele {
    if (hapid.isBlank()) return Allele.NO_CALL
    return Allele.create("<${hapid}>", false)
}

// <R000100_B73_G1>
private fun symbolicAlleleRangeSampleGamete(rangeIndex: Int, sample: String, gamete: Int): Allele {
    return Allele.create("<${symbolicAlleleRangeSampleGameteStr(rangeIndex, sample, gamete)}>", false)
}

private fun symbolicAlleleRangeSampleGameteStr(rangeIndex: Int, sample: String, gamete: Int): String {
    val rangeStr = rangeIndex.toString().padStart(6, '0')
    return "R${rangeStr}_${sample}_G${gamete}"
}
