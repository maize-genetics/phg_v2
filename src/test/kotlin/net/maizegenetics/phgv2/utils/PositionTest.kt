package net.maizegenetics.phgv2.utils

import org.junit.jupiter.api.Test
import kotlin.test.assertEquals

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

        val expectedSortOrder = listOf(
            position1_1,
            position1_100,
            position1A_1,
            position1A_100,
            position1B_1,
            position1B_100,
            position1D_1,
            position1D_100,
            positionChr1_1,
            positionChr1_100,
            position2_1,
            position2_100,
            positionChr2_1,
            positionChr2_100,
            position10_1,
            position10_100,
            positionChr10_1,
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

}