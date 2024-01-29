package net.maizegenetics.phgv2.api

import biokotlin.seq.NucSeqRecord
import biokotlin.seqIO.NucSeqIO
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFAltHeaderLine
import htsjdk.variant.vcf.VCFHeader
import htsjdk.variant.vcf.VCFHeaderLine
import net.maizegenetics.phgv2.utils.altHeaderMetadataToVCFHeaderLine
import net.maizegenetics.phgv2.utils.createGenericHeaderLineSet
import org.apache.logging.log4j.LogManager
import java.io.File

private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.api.ExportHaplotypeGraph")

fun exportMultiSampleHVCF(
    graph: HaplotypeGraph,
    filename: String,
    referenceGenome: String? = null
) {

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

            val headerLines = graph.altHeaders().values
                .map { altHeaderMetadataToVCFHeaderLine(it) }
                .toMutableSet()
            headerLines.addAll(createGenericHeaderLineSet() as Set<VCFAltHeaderLine>)
            val header = VCFHeader(headerLines as Set<VCFHeaderLine>, graph.samples())

            writer.writeHeader(header)

            graph.ranges().forEach { range ->
                val variantContext = createVariantContext(range, graph.hapIdToSampleGametes(range), referenceSequence)
                writer.add(variantContext)
            }

        }

}

private fun createVariantContext(
    range: ReferenceRange,
    hapIdToSampleGametes: Map<String, List<SampleGamete>>,
    referenceSequence: Map<String, NucSeqRecord>?
): VariantContext {

    val alleles = hapIdToSampleGametes.keys
        .asSequence()
        .map { hapid -> Pair(hapid, symbolicAlleleAlt(hapid)) }
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
                val alleles = hapids.map { alleles[it] }
                GenotypeBuilder(taxon, alleles).phased(true).make()
            }
        }
        .toList()

    val resultAlleles = mutableListOf<Allele>()
    resultAlleles.add(alleleRef(range, referenceSequence))
    resultAlleles.addAll(alleles.values)

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
    return Allele.create("<${hapid}>", false)
}