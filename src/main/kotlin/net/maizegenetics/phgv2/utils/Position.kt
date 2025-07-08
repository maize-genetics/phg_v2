package net.maizegenetics.phgv2.utils

/**
 * Data class to represent a position on a contig.
 * Position is 1-based
 */
data class Position (val contig: String, val position: Int) : Comparable<Position> {
    override fun compareTo(other: Position): Int {

        if (this.contig == other.contig) {
            return this.position.compareTo(other.position)
        }

        val thisContig = this.contig.replace("chr", "", ignoreCase = true).trim()
        val otherContig = other.contig.replace("chr", "", ignoreCase = true).trim()

        return try {
            thisContig.toInt() - otherContig.toInt()
        } catch (e: NumberFormatException) {
            // If we can't convert contigs to an int, then compare the strings
            contig.compareTo(other.contig)
        }

    }

    override fun toString(): String {
        return "$contig:$position"
    }
}