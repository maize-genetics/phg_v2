package net.maizegenetics.phgv2.api

/**
 * A ReferenceRange is a contiguous region of a reference genome.
 * The region is defined by a contig, start and end position.
 */
data class ReferenceRange(val contig: String, val start: Int, val end: Int) :
    Comparable<ReferenceRange> {

    override fun compareTo(other: ReferenceRange): Int {
        return when {
            contig != other.contig -> contig.compareTo(other.contig)
            start != other.start -> start.compareTo(other.start)
            else -> end.compareTo(other.end)
        }
    }

    override fun toString(): String {
        return "$contig:$start-$end"
    }

    companion object {
        /**
         * Parse a string of the form "contig:start-end" into a ReferenceRange.
         */
        fun parse(rangeString: String): ReferenceRange {
            val contig = rangeString.substringBefore(":")
            val start = rangeString.substringAfter(":").substringBefore("-").toInt()
            val end = rangeString.substringAfter("-").toInt()
            return ReferenceRange(contig, start, end)
        }
    }
}