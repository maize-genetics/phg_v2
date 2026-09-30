package net.maizegenetics.phgv2.utils

import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class PositionTest {

    @Test
    fun testSortingIntChromosomes() {

        val contig1 = "1"
        val contig2 = "2"
        val contig10 = "10"
        val position1_1 = Position(contig1, 1)
        val position1_100 = Position(contig1, 100)
        val position2_1 = Position(contig2, 1)
        val position2_100 = Position(contig2, 100)
        val position10_1 = Position(contig10, 1)
        val position10_100 = Position(contig10, 100)

        val expectedSortOrder = listOf(
            position1_1,
            position1_100,
            position2_1,
            position2_100,
            position10_1,
            position10_100
        )

        val positions = listOf(position1_1, position1_100, position10_1, position10_100, position2_1, position2_100)

        val sortedPositions = positions.sorted()

        assertEquals(expectedSortOrder, sortedPositions, "Positions are not sorted correctly")

    }

    @Test
    fun testSortingChrIntChromosomes() {

        val contig1 = "chr1"
        val contig2 = "chr2"
        val contig10 = "chr10"
        val position1_1 = Position(contig1, 1)
        val position1_100 = Position(contig1, 100)
        val position2_1 = Position(contig2, 1)
        val position2_100 = Position(contig2, 100)
        val position10_1 = Position(contig10, 1)
        val position10_100 = Position(contig10, 100)

        val expectedSortOrder = listOf(
            position1_1,
            position1_100,
            position2_1,
            position2_100,
            position10_1,
            position10_100
        )

        val positions = listOf(position1_1, position1_100, position10_1, position10_100, position2_1, position2_100)

        val sortedPositions = positions.sorted()

        assertEquals(expectedSortOrder, sortedPositions, "Positions are not sorted correctly")

    }

    @Test
    fun testSortingWheatChromosomes() {

        val contig1A = "1A"
        val contig1B = "1B"
        val contig1D = "1D"
        val position1A_1 = Position(contig1A, 1)
        val position1A_100 = Position(contig1A, 100)
        val position1B_1 = Position(contig1B, 1)
        val position1B_100 = Position(contig1B, 100)
        val position1D_1 = Position(contig1D, 1)
        val position1D_100 = Position(contig1D, 100)

        val expectedSortOrder = listOf(
            position1A_1,
            position1A_100,
            position1B_1,
            position1B_100,
            position1D_1,
            position1D_100
        )

        val positions = listOf(position1A_1, position1A_100, position1D_1, position1D_100, position1B_1, position1B_100)

        val sortedPositions = positions.sorted()

        assertEquals(expectedSortOrder, sortedPositions, "Positions are not sorted correctly")

    }

    @Test
    fun testSortingMixedChromosomes() {

        val contig1 = "1"
        val contig2 = "2"
        val contig10 = "10"
        val position1_1 = Position(contig1, 1)
        val position1_100 = Position(contig1, 100)
        val position2_1 = Position(contig2, 1)
        val position2_100 = Position(contig2, 100)
        val position10_1 = Position(contig10, 1)
        val position10_100 = Position(contig10, 100)

        val contigChr1 = "chr1"
        val contigChr2 = "chr2"
        val contigChr10 = "chr10"
        val positionChr1_1 = Position(contigChr1, 1)
        val positionChr1_100 = Position(contigChr1, 100)
        val positionChr2_1 = Position(contigChr2, 1)
        val positionChr2_100 = Position(contigChr2, 100)
        val positionChr10_1 = Position(contigChr10, 1)
        val positionChr10_100 = Position(contigChr10, 100)

        val contig1A = "1A"
        val contig1B = "1B"
        val contig1D = "1D"
        val position1A_1 = Position(contig1A, 1)
        val position1A_100 = Position(contig1A, 100)
        val position1B_1 = Position(contig1B, 1)
        val position1B_100 = Position(contig1B, 100)
        val position1D_1 = Position(contig1D, 1)
        val position1D_100 = Position(contig1D, 100)

        // "1" and "chr1" are the same contig, so their positions interleave. Positions that differ only
        // in the prefix, such as 1:1 and chr1:1, compare as equal, and the stable sort keeps them in
        // input order. 1A, 1B and 1D follow contig 1 and precede contig 2.
        val expectedSortOrder = listOf(
            position1_1,
            positionChr1_1,
            position1_100,
            positionChr1_100,
            position1A_1,
            position1A_100,
            position1B_1,
            position1B_100,
            position1D_1,
            position1D_100,
            position2_1,
            positionChr2_1,
            position2_100,
            positionChr2_100,
            position10_1,
            positionChr10_1,
            position10_100,
            positionChr10_100
        )

        val positions = listOf(
            position1_1,
            position1_100,
            position10_1,
            position10_100,
            position2_1,
            position2_100,
            positionChr1_1,
            positionChr1_100,
            positionChr2_1,
            positionChr2_100,
            positionChr10_1,
            positionChr10_100,
            position1A_1,
            position1A_100,
            position1B_1,
            position1B_100,
            position1D_1,
            position1D_100
        )

        val sortedPositions = positions.sorted()

        assertEquals(expectedSortOrder, sortedPositions, "Positions are not sorted correctly")

    }

    @Test
    fun testContigsMatchingAfterStrippingChrCompareByPosition() {
        // "chr1" and "1" are the same contig once "chr" is stripped, so the positions must decide.
        // compareTo used to return 0 for any two such positions, whatever their positions.
        assertTrue(Position("chr1", 5) < Position("1", 900), "chr1:5 is before 1:900")
        assertTrue(Position("1", 900) > Position("chr1", 5), "and the reverse agrees")
        assertTrue(Position("1", 5) < Position("chr1", 900), "whichever side carries the prefix")
        assertEquals(0, Position("chr1", 5).compareTo(Position("1", 5)), "same contig, same position")

        // The prefix is stripped case-insensitively, and a non-numeric contig works the same way.
        assertTrue(Position("Chr1", 5) < Position("chr1", 900), "Chr1 and chr1")
        assertTrue(Position("chrX", 5) < Position("X", 900), "chrX and X")
        assertTrue(Position("X", 900) > Position("chrX", 5), "X and chrX, reversed")

        // Different contigs are still ordered by contig first, whatever the positions.
        assertTrue(Position("chr1", 900) < Position("2", 5), "contig 1 sorts before contig 2")

        // The case that matters downstream: a sorted map no longer collapses the two keys into one.
        val map = java.util.TreeMap<Position, String>()
        map[Position("chr1", 5)] = "first"
        map[Position("1", 900)] = "second"
        assertEquals(2, map.size, "two different positions are two keys")
    }

    @Test
    fun testNaturalContigOrder() {
        fun sortedContigs(vararg contigs: String) =
            contigs.map { Position(it, 1) }.sorted().map { it.contig }

        // numbered contigs by number, then the rest of the name; unnumbered contigs after, by name
        assertEquals(listOf("1", "1A", "2", "10", "10A", "Mt", "Pt", "X"),
            sortedContigs("X", "10A", "Pt", "2", "1A", "10", "Mt", "1"))
        // "chr" is ignored in any case, including for unnumbered contigs
        assertEquals(listOf("chr2", "CHR10", "Chr11", "chrUn"), sortedContigs("chrUn", "Chr11", "chr2", "CHR10"))
        assertEquals(listOf("chr1", "chr2", "scaffold_1", "scaffold_10", "scaffold_2"),
            sortedContigs("scaffold_2", "chr2", "scaffold_10", "chr1", "scaffold_1"))

        // leading zeros are ignored, so "01" is contig 1
        assertEquals(0, Position("01", 5).compareTo(Position("chr1", 5)))
        assertTrue(Position("chr01", 900) > Position("1", 5))
        // a number too long for an Int or a Long still compares by value, without overflowing
        assertTrue(Position("99999999999999999999", 1) < Position("100000000000000000000", 1))
        assertTrue(Position("2", 1) < Position("99999999999999999999", 1))
    }

    @Test
    fun testMixedContigNamesSortConsistently() {
        // The old compareTo was not a consistent order for a mix of "1", "chr1" and "1A": it called
        // 1 and chr1 equal but ordered them differently against 1A. Java's sort can detect that and
        // throw "Comparison method violates its general contract!", which it did on two million
        // positions drawn from this set. Distinct positions per contig, so there are no ties and every
        // shuffle must sort to the same order.
        val contigs = (1..10).map { "$it" } + (1..10).map { "chr$it" } + listOf("1A", "1B", "1D", "X", "chrX")
        val positions = contigs.flatMapIndexed { index, contig ->
            (1..20).map { Position(contig, it * 100 + index) }
        }
        val expected = positions.sorted()
        val random = kotlin.random.Random(42)
        repeat(50) {
            assertEquals(expected, positions.shuffled(random).sorted(), "every input order sorts the same")
        }
        // and the order is the natural one
        for ((earlier, later) in expected.zipWithNext()) {
            assertTrue(earlier < later, "$earlier sorts before $later")
        }
    }

}