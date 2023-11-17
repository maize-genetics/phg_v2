package net.maizegenetics.phgv2.api

data class ReferenceRange(val contig: String, val start: Int, val end: Int) :
    Comparable<ReferenceRange> {
    override fun compareTo(other: ReferenceRange): Int {
        return when {
            contig != other.contig -> contig.compareTo(other.contig)
            start != other.start -> start.compareTo(other.start)
            else -> end.compareTo(other.end)
        }
    }
}