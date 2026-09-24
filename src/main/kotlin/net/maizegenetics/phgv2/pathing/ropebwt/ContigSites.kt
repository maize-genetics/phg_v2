package net.maizegenetics.phgv2.pathing.ropebwt

/**
 * Genotypes for one contig, from a reference panel and the samples being imputed, paired site by
 * site and encoded as allele indices.
 *
 * Both VCFs are reduced to small integers here because the emission model only ever asks whether two
 * alleles are the *same*, never what they are. The panel and the samples must therefore share one
 * index per site, which is why they are encoded together: the index is assigned from the union of the
 * allele strings both carry at that site, so comparing a founder's allele to a sample's allele is an
 * integer comparison and is independent of the order the two VCFs happen to list their ALTs in.
 *
 * Only sites present in **both** the panel and the sample VCF appear. A site the panel does not carry
 * cannot say anything about founders, and a site the sample does not carry has nothing to explain.
 *
 * ## Layout
 *
 * Two halves are stored per genotype, so a heterozygous call is preserved rather than collapsed. A
 * haploid call, or a homozygous diploid one, has both halves equal. [MISSING] marks a no-call.
 *
 *     founderAlleles[(site * nFounders + founder) * 2 + half]
 *     sampleAlleles [(site * nSamples  + sample ) * 2 + half]
 *
 * ## Size
 *
 * `2 * nSites * (nFounders + nSamples)` bytes for a contig. A low-density panel -- the intended
 * input, a few tens of thousands of sites -- is a few megabytes. A dense panel used for path
 * inference as well as for composition is far larger: 16 M sites x 25 founders x 2 is 800 MB for the
 * founders alone, plus as much again per hundred samples. If that becomes the normal case, the two
 * halves can collapse to one byte per founder wherever the panel is homozygous, with a side table for
 * the heterozygous minority, which halves the founder term.
 *
 * @param founderNames panel sample names, in the order the founder axis is indexed
 * @param sampleNames names of the samples being imputed, in the order the sample axis is indexed
 * @param positions reference positions of the sites, ascending; used for distance-scaled transitions
 *   and to convert a path back to reference coordinates
 */
class ContigSites(
    val founderNames: List<String>,
    val sampleNames: List<String>,
    val positions: IntArray,
    val founderAlleles: ByteArray,
    val sampleAlleles: ByteArray
) {
    val nFounders = founderNames.size
    val nSamples = sampleNames.size
    val nSites = positions.size

    init {
        require(founderAlleles.size == nSites * nFounders * 2) {
            "founderAlleles must hold 2 halves for each of $nSites sites x $nFounders founders, " +
                    "but holds ${founderAlleles.size}"
        }
        require(sampleAlleles.size == nSites * nSamples * 2) {
            "sampleAlleles must hold 2 halves for each of $nSites sites x $nSamples samples, " +
                    "but holds ${sampleAlleles.size}"
        }
    }

    /** First allele index of [founder] at [site], or [MISSING]. */
    fun founderAllele1(site: Int, founder: Int) = founderAlleles[(site * nFounders + founder) * 2]

    /** Second allele index of [founder] at [site]; equal to the first unless the founder is heterozygous. */
    fun founderAllele2(site: Int, founder: Int) = founderAlleles[(site * nFounders + founder) * 2 + 1]

    /** First allele index of [sample] at [site], or [MISSING]. */
    fun sampleAllele1(site: Int, sample: Int) = sampleAlleles[(site * nSamples + sample) * 2]

    /** Second allele index of [sample] at [site]; equal to the first when the call is homozygous. */
    fun sampleAllele2(site: Int, sample: Int) = sampleAlleles[(site * nSamples + sample) * 2 + 1]

    companion object {
        /** No call, for either a founder or a sample. */
        const val MISSING: Byte = -1

        /**
         * The most distinct alleles one site may carry. Allele indices are stored in a byte, and
         * [MISSING] takes -1, so the usable range is 0..126. No real site approaches this; the limit
         * exists so that a pathological one fails loudly rather than wrapping round to [MISSING].
         */
        const val MAX_ALLELES_PER_SITE = 127
    }
}
