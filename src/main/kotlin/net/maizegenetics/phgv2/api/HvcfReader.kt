package net.maizegenetics.phgv2.api

/**
 * An Interface defining methods to return genotypes from a data source for requested samples and ranges
 * or, if sample names and ranges are not set, for all the data from the source.
 * Samples are identified by name. Ranges are defined by the PositionRange.
 * Other methods return either a single header or a sample -> header map.
 *
 * For now header is returned as a String, but it would make sense to implement a Header interface with methods to
 * return specific information. The methods would probably look like a subset of htsjdk header methods.
 *
 * Typical usage:
 * val myReader = HvcfReader().dataSource("data uri")
 * val someData = myReader.samples(sampleNameList).range(positionRangeList).data()
 * Subsequent calls for more data chunks would not need to set samples again.
 *
 * The interface could support tiledb-vcf, a single hvcf file or a set of hvcf files as data sources.
 *
 * Are there other attributes needed in addition to genotype?
 *
 * Note: there is no method in the interface for setting a datasource as that is left to the implementation to handle.
 */
interface HvcfReader {

    /**
     * Sets the range to be returned. Range is a list of PositionRanges
     */
    fun range(range: List<PositionRange>): HvcfReader

    /**
     * Returns the header for a sample as a String
     */
    fun header(sampleName: String): String

    /**
     * Returns the headers for all samples in this source as a
     */
    fun headerMap(): Map<String, String>

    /**
     * Returns the list of sample names that will be returned by this HvcfReader
     */
    fun sampleNames(): List<String>

    /**
     * Returns data for the requested source, sample names, and ranges
     */
    fun data(): List<SampleData>

    /**
     * A position range and genotype. Genotype is the allele/haplotype id not the integer GT code.
     * genotype and AD are lists to accommodate different ploidy levels
     */
    data class SampleData(val sampleName: String, val contig: String, val startPos: Int, val endPos: Int, val genotype: List<String>, val AD: List<Int>? = null, val DP: Int? = null)


    /**
     * A range of genomic positions. If start and end positions are null then it represents an entire chromosome.
     */
    data class PositionRange(val contig: String, val startPos: Int? = null, val endPos: Int? = null) {
        override fun toString(): String {
            return if (startPos == null && endPos == null) "$contig"
            else if (endPos == null) "$contig:$startPos-$startPos"
            else if (startPos == null) throw IllegalArgumentException("startPos is null and endPos = $endPos")
            else "$contig:$startPos-$endPos"
        }
    }

}