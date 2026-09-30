package net.maizegenetics.phgv2.utils

/**
 * Data class to represent a position on a contig.
 * Position is 1-based
 *
 * Positions are ordered by contig, then by position. Contigs are ordered naturally: a leading "chr"
 * (any case) is ignored, contigs whose name then starts with a number come first, ordered by that
 * number and then by the rest of the name, and all other contigs follow, ordered by name. So
 * "1" and "chr1" are the same contig, and the order runs 1, 1A, 1B, 2, 10, and then X, Pt, scaffolds.
 *
 * The comparison is a consistent total order, which sorting and sorted collections such as TreeMap and
 * TreeRangeMap depend on, and it allocates nothing, since compareTo runs on very large data sets.
 *
 * Note that compareTo is not consistent with equals when two positions differ only in the "chr"
 * prefix: Position("chr1", 5) and Position("1", 5) compare as 0 but are not equal, so a sorted
 * collection would treat them as the same key.
 */
data class Position (val contig: String, val position: Int) : Comparable<Position> {
    override fun compareTo(other: Position): Int {

        if (this.contig == other.contig) {
            return this.position.compareTo(other.position)
        }

        val contigComparison = compareContigs(this.contig, other.contig)
        return if (contigComparison != 0) contigComparison else this.position.compareTo(other.position)

    }

    override fun toString(): String {
        return "$contig:$position"
    }

    companion object {

        /**
         * Compares two contig names in the natural order described on the class.
         *
         * Works on the strings in place, without substrings or parsing: this is on the path of every
         * comparison between positions on different contigs. The leading number is compared by its
         * digits rather than converted, so a number of any length compares correctly and cannot
         * overflow, and leading zeros are ignored, so "01" and "1" are the same contig.
         */
        fun compareContigs(a: String, b: String): Int {
            val aStart = if (a.startsWith("chr", ignoreCase = true)) 3 else 0
            val bStart = if (b.startsWith("chr", ignoreCase = true)) 3 else 0

            var aDigitsEnd = aStart
            while (aDigitsEnd < a.length && a[aDigitsEnd] in '0'..'9') aDigitsEnd++
            var bDigitsEnd = bStart
            while (bDigitsEnd < b.length && b[bDigitsEnd] in '0'..'9') bDigitsEnd++

            val aHasNumber = aDigitsEnd > aStart
            val bHasNumber = bDigitsEnd > bStart
            if (aHasNumber != bHasNumber) return if (aHasNumber) -1 else 1

            var ai = aStart
            var bi = bStart
            if (aHasNumber) {
                // skip leading zeros, then the longer run of digits is the larger number, and runs of
                // equal length compare digit by digit
                while (ai < aDigitsEnd - 1 && a[ai] == '0') ai++
                while (bi < bDigitsEnd - 1 && b[bi] == '0') bi++
                val byLength = (aDigitsEnd - ai).compareTo(bDigitsEnd - bi)
                if (byLength != 0) return byLength
                while (ai < aDigitsEnd) {
                    val byDigit = a[ai].compareTo(b[bi])
                    if (byDigit != 0) return byDigit
                    ai++
                    bi++
                }
            }

            // the rest of the name, or the whole name after "chr" when there is no leading number
            while (ai < a.length && bi < b.length) {
                val byChar = a[ai].compareTo(b[bi])
                if (byChar != 0) return byChar
                ai++
                bi++
            }
            return (a.length - ai).compareTo(b.length - bi)
        }
    }
}
