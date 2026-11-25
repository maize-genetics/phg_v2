package net.maizegenetics.phgv2.pathing.ropebwt

import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFFormatHeaderLine
import htsjdk.variant.vcf.VCFHeaderLineType
import java.io.File

/**
 * Writes a sample's own genotype call over the composed one, wherever the sample has a full call at a
 * site of the panel being composed from. This is what `--keep-sample-calls` does.
 *
 * The path reproduces the sample's calls at most of its own sites, but not all of them: a call the
 * path cannot explain -- a genotyping error, a new mutation, a founder allele the panel has wrong -- is
 * outvoted by its neighbours. Where the sample was genotyped, its call is direct evidence, so this
 * keeps it. Every genotype then records where it came from in `FORMAT/SRC`: `O` for observed, the
 * sample's own call, and `I` for imputed.
 *
 * ## Which calls are used
 *
 * A sample record is matched to a panel record on contig, position and REF, as the path pairing
 * matches them; see [pairVcfSites]. Its call is used only when
 *
 *  - it is a **full call**: both alleles of a diploid genotype, or the one allele of a haploid one,
 *    which is written as homozygous, as the pairing reads it. A half call such as `0/.` is not used,
 *    since the missing half could be either allele.
 *  - every allele it calls is one the panel record carries, matched by its bases. Otherwise the call
 *    is counted in [SampleCallCounts.unmatchedAlleles] and the composed genotype is kept: the record
 *    has no way to write an allele it does not list.
 *
 * Where several sample records share the key -- a multi-allelic site split across records -- the
 * first whose alleles the panel record carries is used.
 *
 * ## Ordering
 *
 * The sample VCF is read once, in step with the panel, so both must be sorted the same way: contigs in
 * [contigRank]'s order and positions ascending within each. A contig the panel does not have is
 * skipped, since none of its records can match. A record out of order fails, since a merge-join over
 * it would quietly match almost nothing.
 *
 * @param sampleVcf the VCF of the samples being imputed
 * @param contigRank the order of the panel's contigs, as from [panelContigOrder]
 */
class SampleCallOverlay(sampleVcf: File, private val contigRank: Map<String, Int>) : AutoCloseable {

    private val reader = VCFFileReader(sampleVcf, false)
    private val iterator = reader.iterator()
    private val name = sampleVcf.name
    private var lookahead: VariantContext? = null
    private var lastRank = -1
    private var lastPosition = -1

    /** Sample records at [bufferRank]:[bufferPosition], kept while the panel has records there too. */
    private val buffer = mutableListOf<VariantContext>()
    private var bufferRank = -1
    private var bufferPosition = -1

    /** What was done for each sample, by name. */
    val counts = mutableMapOf<String, SampleCallCounts>()

    init { advance() }

    /**
     * The sample record to take calls from at [variant]'s site, or null when the sample VCF has none
     * there. Must be called with panel records in sorted order.
     */
    fun recordAt(variant: VariantContext): VariantContext? {
        val rank = contigRank[variant.contig] ?: return null
        val position = variant.start
        if (rank != bufferRank || position != bufferPosition) {
            buffer.clear()
            bufferRank = rank; bufferPosition = position
            while (true) {
                val next = lookahead ?: break
                val nextRank = contigRank.getValue(next.contig)
                if (nextRank > rank || (nextRank == rank && next.start > position)) break
                if (nextRank == rank && next.start == position) buffer.add(next)
                advance()
            }
        }
        val sameRef = buffer.filter { it.reference.baseString == variant.reference.baseString }
        return sameRef.firstOrNull { record -> record.alleles.all { variant.hasAlleleBases(it) } }
            ?: sameRef.firstOrNull()
    }

    /**
     * The genotype to write for [sampleName] at [variant]: the sample's own call from [observed] when
     * it can be used, otherwise [composed], which is null where the path does not cover the site. The
     * result carries `FORMAT/SRC`, and [counts] records what happened.
     */
    fun overlay(variant: VariantContext, observed: VariantContext?, sampleName: String,
                composed: List<Allele>?): Genotype? {
        val count = counts.getOrPut(sampleName) { SampleCallCounts() }
        val called = observed?.getGenotype(sampleName)?.alleles ?: emptyList()
        val observedAlleles: List<Allele>? = when {
            called.isEmpty() || called.all { it.isNoCall } -> null
            called.any { it.isNoCall } -> { count.halfCalls++; null }
            else -> {
                val mapped = called.map { allele -> variant.alleles.firstOrNull { it.baseString == allele.baseString } }
                if (mapped.any { it == null }) { count.unmatchedAlleles++; null }
                else if (mapped.size == 1) listOf(mapped[0]!!, mapped[0]!!) else mapped.map { it!! }
            }
        }

        if (observedAlleles == null) {
            return composed?.let { GenotypeBuilder(sampleName).alleles(it).attribute(SOURCE_KEY, IMPUTED).make() }
        }

        count.sampleCalls++
        when {
            composed == null || composed.any { it.isNoCall } -> count.filled++
            sameGenotype(composed, observedAlleles) -> count.agreed++
            else -> count.overwritten++
        }
        return GenotypeBuilder(sampleName).alleles(observedAlleles).attribute(SOURCE_KEY, OBSERVED).make()
    }

    private fun sameGenotype(a: List<Allele>, b: List<Allele>): Boolean =
        a.map { it.baseString }.sorted() == b.map { it.baseString }.sorted()

    private fun VariantContext.hasAlleleBases(allele: Allele) =
        allele.isNoCall || alleles.any { it.baseString == allele.baseString }

    /** Moves the lookahead to the next record on a contig the panel has, checking the order. */
    private fun advance() {
        while (iterator.hasNext()) {
            val record = iterator.next()
            val rank = contigRank[record.contig] ?: continue
            if (rank == lastRank) {
                require(record.start >= lastPosition) {
                    "$name is not coordinate sorted: ${record.contig}:${record.start} follows position " +
                            "$lastPosition, so its calls could not be matched to the panel's sites."
                }
            } else {
                require(rank > lastRank) {
                    "$name lists contig ${record.contig} out of the panel's contig order, so its calls " +
                            "could not be matched to the panel's sites."
                }
                lastRank = rank
            }
            lastPosition = record.start
            lookahead = record
            return
        }
        lookahead = null
    }

    override fun close() = reader.close()

    companion object {
        const val SOURCE_KEY = "SRC"
        const val OBSERVED = "O"
        const val IMPUTED = "I"

        /** The `FORMAT/SRC` header line. */
        val sourceHeaderLine = VCFFormatHeaderLine(SOURCE_KEY, 1, VCFHeaderLineType.String,
            "Source of the genotype: O, the sample's own call (--keep-sample-calls); I, imputed from the founder path")
    }
}

/**
 * What [SampleCallOverlay] did for one sample, over every site of the composed VCF.
 *
 * @param sampleCalls sites where the sample's own call was written
 * @param agreed of those, sites where the composed genotype was the same
 * @param overwritten of those, sites where the composed genotype was a different full call
 * @param filled of those, sites where there was no composed genotype or only a partial one
 * @param unmatchedAlleles full calls not used, because they call an allele the panel record lacks
 * @param halfCalls half calls not used
 */
data class SampleCallCounts(
    var sampleCalls: Long = 0,
    var agreed: Long = 0,
    var overwritten: Long = 0,
    var filled: Long = 0,
    var unmatchedAlleles: Long = 0,
    var halfCalls: Long = 0
) {
    /**
     * Overwritten as a fraction of the sites where both a sample call and a full composed call
     * existed: how often the path disagrees with what was genotyped. An optimistic estimate of
     * imputation error, since those sites also informed the path. NaN where there were none.
     */
    val discordance: Double get() = overwritten.toDouble() / (agreed + overwritten)
}

/**
 * Writes [counts] as a tab-delimited table, one row per sample in name order.
 */
fun writeSampleCallSummary(counts: Map<String, SampleCallCounts>, file: File) {
    file.bufferedWriter().use { out ->
        out.write("sample\tsample_calls\tagreed\toverwritten\tfilled\tunmatched_alleles\thalf_calls\tdiscordance\n")
        for ((sample, c) in counts.toSortedMap()) {
            val discordance = if (c.agreed + c.overwritten == 0L) "NA" else "%.6f".format(c.discordance)
            out.write("$sample\t${c.sampleCalls}\t${c.agreed}\t${c.overwritten}\t${c.filled}\t" +
                    "${c.unmatchedAlleles}\t${c.halfCalls}\t$discordance\n")
        }
    }
}
