package net.maizegenetics.phgv2.pathing.ropebwt

import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import org.apache.logging.log4j.LogManager
import net.maizegenetics.phgv2.pathing.ropebwt.ContigSites.Companion.MAX_ALLELES_PER_SITE
import net.maizegenetics.phgv2.pathing.ropebwt.ContigSites.Companion.MISSING
import java.io.File

/**
 * What the pairing saw, for the log. A path is only as good as the sites behind it, and these are the
 * numbers that say whether the two files actually correspond.
 *
 * @param sitesUsed sites present in both files and carried into the model
 * @param sampleSitesNotInPanel sites the sample VCF carries that the panel does not. The panel cannot
 *   say which founder any allele belongs to there, so they are dropped. A large count is the signature
 *   of a mismatch -- wrong reference build, wrong contig naming, or a panel that does not cover the
 *   assay -- and is otherwise invisible, since the sites are simply skipped.
 * @param panelSitesNotInSample sites only the panel carries. Expected and harmless when the sample
 *   VCF is the sparser of the two, which is the normal case.
 * @param duplicateKeysSkipped extra records sharing a contig, position and REF allele with one
 *   already taken; see [pairVcfSites] on why the first is used.
 */
data class PairingCounts(
    var sitesUsed: Long = 0,
    var sampleSitesNotInPanel: Long = 0,
    var panelSitesNotInSample: Long = 0,
    var duplicateKeysSkipped: Long = 0
)

/**
 * Walks a sample VCF and a reference panel VCF together and hands [perContig] the paired, allele-index
 * encoded sites for each contig in turn.
 *
 * Both files are coordinate sorted, so one pass over each suffices and the panel is never held beyond
 * the contig in hand. That matters because the "low density" panel is only conceptually low density:
 * nothing stops a user passing the same dense panel they use for composition.
 *
 * ## Which sites pair
 *
 * A site is carried only when **both** files have it. A site the panel lacks cannot say which founder
 * any allele belongs to; a site the sample lacks has nothing to explain.
 *
 * Records are matched on **contig, position and REF allele**, not position alone. A SNP and an
 * insertion can sit at one position, and they are different variants that must not be conflated; the
 * REF allele separates them. Where several records still share that key -- a multi-allelic site split
 * across records, which is the case REF cannot distinguish -- the first is used and the rest counted
 * in [PairingCounts.duplicateKeysSkipped]. A multi-allelic site expressed the usual way, as one record
 * with several ALTs, needs none of this and is handled in full.
 *
 * ## Allele encoding
 *
 * The emission model only asks whether two alleles are the *same*, so alleles become small integers.
 * Indices are assigned per site from the union of the base strings both records carry, which makes the
 * comparison independent of the order either file happens to list its ALTs in -- the two VCFs are
 * unrelated files and need not agree on that.
 *
 * A genotype of one allele, whether a haploid call or a panel founder written that way, is encoded as
 * homozygous. That is what makes a haploid founder contribute a single gamete downstream.
 *
 * ## Ordering
 *
 * The contig order comes from the panel's own header. Every contig the sample VCF mentions must appear
 * there, and positions must ascend within a contig in both files; a merge-join over files that
 * disagree would silently pair almost nothing rather than fail, so both are checked.
 */
fun pairVcfSites(
    sampleVcf: File,
    panelVcf: File,
    contigsToUse: Set<String> = emptySet(),
    perContig: (ContigSites) -> Unit
): PairingCounts {
    val counts = PairingCounts()

    VCFFileReader(panelVcf, false).use { panelReader ->
        VCFFileReader(sampleVcf, false).use { sampleReader ->
            val founderNames = panelReader.fileHeader.sampleNamesInOrder.toList()
            val sampleNames = sampleReader.fileHeader.sampleNamesInOrder.toList()
            require(founderNames.isNotEmpty()) { "$panelVcf has no samples to use as founders" }
            require(sampleNames.isNotEmpty()) { "$sampleVcf has no samples to impute" }

            val contigRank = panelContigOrder(panelReader, panelVcf)

            val accumulator = SiteAccumulator(founderNames, sampleNames, perContig)
            val panel = PeekingVariants(panelReader.iterator(), panelVcf.name, contigRank)
            val sample = PeekingVariants(sampleReader.iterator(), sampleVcf.name, contigRank)

            while (panel.peek() != null && sample.peek() != null) {
                val panelRecord = panel.peek()!!
                val sampleRecord = sample.peek()!!
                val comparison = compareValuesBy(panelRecord, sampleRecord,
                    { contigRank.getValue(it.contig) }, { it.start })
                when {
                    comparison < 0 -> { counts.panelSitesNotInSample++; panel.next() }
                    comparison > 0 -> { counts.sampleSitesNotInPanel++; sample.next() }
                    // Same contig and position. REF still has to agree, or these are different
                    // variants that merely share a coordinate.
                    panelRecord.reference.baseString != sampleRecord.reference.baseString -> {
                        // Advance the one whose REF sorts first, so a run of records at one position
                        // is worked through rather than deadlocking.
                        if (panelRecord.reference.baseString < sampleRecord.reference.baseString) {
                            counts.panelSitesNotInSample++; panel.next()
                        } else {
                            counts.sampleSitesNotInPanel++; sample.next()
                        }
                    }
                    else -> {
                        val useContig = contigsToUse.isEmpty() || contigsToUse.contains(panelRecord.contig)
                        if (useContig) {
                            accumulator.add(panelRecord, sampleRecord)
                            counts.sitesUsed++
                        }
                        panel.next()
                        sample.next()
                        // Any further records sharing this key cannot be told apart from the one just
                        // taken, so they are dropped and counted.
                        counts.duplicateKeysSkipped += panel.skipSameKeyAs(panelRecord)
                        counts.duplicateKeysSkipped += sample.skipSameKeyAs(sampleRecord)
                    }
                }
            }
            while (panel.peek() != null) { counts.panelSitesNotInSample++; panel.next() }
            while (sample.peek() != null) { counts.sampleSitesNotInPanel++; sample.next() }
            accumulator.flush()
        }
    }
    return counts
}

/**
 * The contig order the merge-join walks in, taken from the panel.
 *
 * `##contig` header lines are preferred, but they are optional in the VCF specification and real panels
 * do without them -- the 25-founder maize panel this pipeline targets has none. The fallback reads the
 * panel once to learn the order its contigs first appear in, which for a coordinate-sorted file is the
 * same thing. That costs a pass over the file, so it is logged rather than done silently.
 */
private fun panelContigOrder(panelReader: VCFFileReader, panelVcf: File): Map<String, Int> {
    val fromHeader = panelReader.fileHeader.contigLines
        .withIndex().associate { (index, line) -> line.id to index }
    if (fromHeader.isNotEmpty()) return fromHeader

    LogManager.getLogger("net.maizegenetics.phgv2.pathing.ropebwt.VcfSitePairing").info(
        "${panelVcf.name} has no ##contig header lines; reading it once to learn its contig order"
    )
    val order = LinkedHashMap<String, Int>()
    VCFFileReader(panelVcf, false).use { scan ->
        for (record in scan) order.getOrPut(record.contig) { order.size }
    }
    require(order.isNotEmpty()) { "$panelVcf holds no variant records" }
    return order
}

/**
 * One-record lookahead over a VCF, checking as it goes that the file is sorted the way the merge-join
 * assumes: contigs in the panel header's order, positions ascending within a contig.
 */
private class PeekingVariants(
    private val iterator: Iterator<VariantContext>,
    private val name: String,
    private val contigRank: Map<String, Int>
) {
    private var current: VariantContext? = if (iterator.hasNext()) iterator.next() else null
    private var lastRank = -1
    private var lastPosition = -1

    init { current?.let { check(it) } }

    fun peek(): VariantContext? = current

    fun next(): VariantContext? {
        current = if (iterator.hasNext()) iterator.next() else null
        current?.let { check(it) }
        return current
    }

    /** Advances past any record sharing [record]'s contig, position and REF, returning how many. */
    fun skipSameKeyAs(record: VariantContext): Long {
        var skipped = 0L
        while (current?.let {
                it.contig == record.contig && it.start == record.start &&
                        it.reference.baseString == record.reference.baseString
            } == true) {
            skipped++
            next()
        }
        return skipped
    }

    private fun check(record: VariantContext) {
        val rank = contigRank[record.contig] ?: throw IllegalArgumentException(
            "$name has contig ${record.contig}, which is not in the panel's ##contig header lines. " +
                    "The two files must use the same contig names."
        )
        if (rank == lastRank) {
            require(record.start >= lastPosition) {
                "$name is not coordinate sorted: ${record.contig}:${record.start} follows " +
                        "position $lastPosition. A merge-join over unsorted input would pair almost " +
                        "nothing rather than fail."
            }
        } else {
            require(rank > lastRank) {
                "$name lists contig ${record.contig} out of the panel's header order. The two files " +
                        "must order their contigs the same way."
            }
            lastRank = rank
        }
        lastPosition = record.start
    }
}

/**
 * Accumulates paired sites for the contig in hand, emitting a [ContigSites] when the contig changes.
 */
private class SiteAccumulator(
    private val founderNames: List<String>,
    private val sampleNames: List<String>,
    private val emit: (ContigSites) -> Unit
) {
    private var contig: String? = null
    private var positions = IntArray(1024)
    private var siteCount = 0
    private var founderAlleles = ByteArray(1024 * founderNames.size * 2)
    private var sampleAlleles = ByteArray(1024 * sampleNames.size * 2)

    fun add(panelRecord: VariantContext, sampleRecord: VariantContext) {
        if (contig != panelRecord.contig) {
            flush()
            contig = panelRecord.contig
        }
        grow()

        // One index per distinct allele string, from the union of the two records. REF is shared, so
        // it always lands on 0.
        val alleleIndex = HashMap<String, Byte>(8)
        fun indexOf(baseString: String): Byte = alleleIndex.getOrPut(baseString) {
            require(alleleIndex.size < MAX_ALLELES_PER_SITE) {
                "${panelRecord.contig}:${panelRecord.start} has more than " +
                        "$MAX_ALLELES_PER_SITE distinct alleles, which an allele index cannot hold"
            }
            alleleIndex.size.toByte()
        }
        indexOf(panelRecord.reference.baseString)

        encode(panelRecord, founderNames, founderAlleles,
            (siteCount * founderNames.size) * 2, ::indexOf)
        encode(sampleRecord, sampleNames, sampleAlleles,
            (siteCount * sampleNames.size) * 2, ::indexOf)
        positions[siteCount] = panelRecord.start
        siteCount++
    }

    /**
     * Writes each named sample's genotype as two allele halves. A genotype of one allele is stored as
     * homozygous, so a haploid call contributes a single gamete; no call, or a sample the record does
     * not mention, is [MISSING] in both halves.
     */
    private fun encode(
        record: VariantContext,
        names: List<String>,
        into: ByteArray,
        offset: Int,
        indexOf: (String) -> Byte
    ) {
        names.forEachIndexed { index, name ->
            val slot = offset + index * 2
            val genotype = record.getGenotype(name)
            val called = genotype?.alleles?.filter { !it.isNoCall } ?: emptyList()
            when {
                called.isEmpty() -> { into[slot] = MISSING; into[slot + 1] = MISSING }
                called.size == 1 -> {
                    val allele = indexOf(called[0].baseString)
                    into[slot] = allele; into[slot + 1] = allele
                }
                else -> {
                    into[slot] = indexOf(called[0].baseString)
                    into[slot + 1] = indexOf(called[1].baseString)
                }
            }
        }
    }

    private fun grow() {
        if (siteCount < positions.size) return
        positions = positions.copyOf(positions.size * 2)
        founderAlleles = founderAlleles.copyOf(founderAlleles.size * 2)
        sampleAlleles = sampleAlleles.copyOf(sampleAlleles.size * 2)
    }

    fun flush() {
        val name = contig ?: return
        if (siteCount > 0) {
            emit(
                ContigSites(
                    name, founderNames, sampleNames,
                    positions.copyOf(siteCount),
                    founderAlleles.copyOf(siteCount * founderNames.size * 2),
                    sampleAlleles.copyOf(siteCount * sampleNames.size * 2)
                )
            )
        }
        siteCount = 0
        contig = null
    }
}
