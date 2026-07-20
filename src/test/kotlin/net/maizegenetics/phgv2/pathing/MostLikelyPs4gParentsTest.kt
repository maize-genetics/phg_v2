package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader
import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.random.Random
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNull
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
class MostLikelyPs4gParentsTest {

    companion object {
        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.testOutputDir).deleteRecursively()
            File(TestExtension.testOutputDir).mkdirs()
        }

        @JvmStatic
        @AfterAll
        fun tearDown() {
            //comment out the following line to inspect the test results after the tests have been run
            File(TestExtension.testOutputDir).deleteRecursively()
        }
    }

    /** One body line of a PS4G file: the gametes hit, where the hit was, and how many reads hit. */
    private data class Ps4gLine(val gameteIndices: List<Int>, val contig: String, val pos: Int, val count: Int)

    /**
     * Writes a PS4G (version 2.0) file holding [lines], and returns its path. The per-gamete counts in the gamete
     * index block are summed from [lines] so that the file is internally consistent, though [MostLikelyPs4gParents]
     * recounts from the body rather than reading those counts.
     */
    private fun writePs4gFile(name: String, numberOfGametes: Int, lines: List<Ps4gLine>): String {
        val filename = "${TestExtension.testOutputDir}$name.ps4g"
        File(filename).bufferedWriter(Charsets.UTF_8).use { writer ->
            writer.write("#PS4G\n")
            writer.write("#version=2.0\n")
            writer.write("#Command: PHGv2 Command:\n")
            writer.write("#gamete\tgameteIndex\tcount\n")
            for (ndx in 0 until numberOfGametes) {
                val gameteCount = lines.filter { it.gameteIndices.contains(ndx) }.sumOf { it.count }
                writer.write("#gamete$ndx\t$ndx\t$gameteCount\n")
            }
            writer.write("gameteSet\trefContig\trefPosBinned\tcount\n")
            for (line in lines) {
                writer.write("${line.gameteIndices.joinToString(",")}\t${line.contig}\t${line.pos}\t${line.count}\n")
            }
        }
        return filename
    }

    /** A [MostLikelyPs4gParents] over a trivial PS4G file, for testing methods that take their data as arguments. */
    private fun parentFinderForHelpers(name: String): MostLikelyPs4gParents {
        val filename = writePs4gFile(name, 1, listOf(Ps4gLine(listOf(0), "chr1", 0, 1)))
        return MostLikelyPs4gParents(Ps4gFileReader(filename), setOf("chr1"))
    }

    /**
     * Read mappings where gamete 1 has the second highest read count but only ever appears alongside gamete 0,
     * while gamete 2 has a lower count but is hit by reads that gamete 0 misses. Raw counts are
     * 0 -> 15, 1 -> 10, 2 -> 6, 3 -> 1.
     */
    private fun complementaryParentLines(contig: String = "chr1") = listOf(
        Ps4gLine(listOf(0, 1), contig, 100, 10),
        Ps4gLine(listOf(0), contig, 200, 5),
        Ps4gLine(listOf(2), contig, 300, 6),
        Ps4gLine(listOf(3), contig, 400, 1)
    )

    @Test
    fun testSecondParentIsComplementaryRatherThanSecondHighestCount() {
        val filename = writePs4gFile("complementary", 4, complementaryParentLines())
        val parentFinder = MostLikelyPs4gParents(Ps4gFileReader(filename), setOf("chr1"))

        //gamete 0 has the highest count, so it is chosen first. Gamete 1 has the next highest count, but every read
        //hitting it also hits gamete 0, so gamete 2 is the better complement and should be chosen second.
        assertEquals(setOf(0, 2), parentFinder.bestParents(2))
    }

    @Test
    fun testAdditionalParentsAddedByDecreasingCount() {
        val filename = writePs4gFile("threeParents", 4, complementaryParentLines())
        val parentFinder = MostLikelyPs4gParents(Ps4gFileReader(filename), setOf("chr1"))

        //after 0 and 2 are chosen, gamete 1 has the highest remaining count, then gamete 3
        assertEquals(setOf(0, 1, 2), parentFinder.bestParents(3))
        assertEquals(setOf(0, 1, 2, 3), parentFinder.bestParents(4))
    }

    @Test
    fun testOnlyContigsInChromosomeSetAreCounted() {
        //chr1 on its own picks gametes 0 and 1; chr2 has the larger counts and picks gametes 2 and 3
        val lines = listOf(
            Ps4gLine(listOf(0), "chr1", 100, 10),
            Ps4gLine(listOf(1), "chr1", 200, 1),
            Ps4gLine(listOf(2), "chr2", 100, 20),
            Ps4gLine(listOf(3), "chr2", 200, 15)
        )
        val filename = writePs4gFile("twoContigs", 4, lines)
        val ps4gReader = Ps4gFileReader(filename)

        assertEquals(setOf(0, 1), MostLikelyPs4gParents(ps4gReader, setOf("chr1")).bestParents(2))
        assertEquals(setOf(2, 3), MostLikelyPs4gParents(ps4gReader, setOf("chr2")).bestParents(2))

        //with both contigs the chr2 gametes outcount the chr1 gametes
        assertEquals(setOf(2, 3), MostLikelyPs4gParents(ps4gReader, setOf("chr1", "chr2")).bestParents(2))
    }

    @Test
    fun testCountsAreSummedAcrossPositionsAndContigs() {
        //no single position makes gamete 0 the winner, but its counts sum to 9 against 6 and 5 for gametes 1 and 2
        val lines = listOf(
            Ps4gLine(listOf(0), "chr1", 100, 3),
            Ps4gLine(listOf(0), "chr1", 200, 3),
            Ps4gLine(listOf(0), "chr2", 100, 3),
            Ps4gLine(listOf(1), "chr1", 300, 6),
            Ps4gLine(listOf(2), "chr2", 200, 5)
        )
        val filename = writePs4gFile("summedCounts", 3, lines)
        val parentFinder = MostLikelyPs4gParents(Ps4gFileReader(filename), setOf("chr1", "chr2"))

        assertEquals(setOf(0, 1), parentFinder.bestParents(2))
    }

    @Test
    fun testMultipleGameteSetsAtSamePosition() {
        //a single binned position can hold more than one gamete set; all of them must be counted
        val lines = listOf(
            Ps4gLine(listOf(0), "chr1", 100, 4),
            Ps4gLine(listOf(1), "chr1", 100, 3),
            Ps4gLine(listOf(2), "chr1", 100, 1)
        )
        val filename = writePs4gFile("samePosition", 3, lines)
        val parentFinder = MostLikelyPs4gParents(Ps4gFileReader(filename), setOf("chr1"))

        assertEquals(setOf(0, 1), parentFinder.bestParents(2))
    }

    @Test
    fun testNumberOfParentsMustBeGreaterThanOne() {
        val filename = writePs4gFile("tooFewParents", 4, complementaryParentLines())
        val parentFinder = MostLikelyPs4gParents(Ps4gFileReader(filename), setOf("chr1"))

        for (numberOfParents in listOf(1, 0, -1)) {
            val exception = assertFailsWith<IllegalArgumentException> { parentFinder.bestParents(numberOfParents) }
            assertEquals("Number of best parents must be greater than 1.", exception.message)
        }
    }

    @Test
    fun testContigNotInPs4gFileThrows() {
        val filename = writePs4gFile("missingContig", 4, complementaryParentLines())
        val parentFinder = MostLikelyPs4gParents(Ps4gFileReader(filename), setOf("chr1", "chr99"))

        val exception = assertFailsWith<IllegalStateException> { parentFinder.bestParents(2) }
        assertTrue(exception.message!!.contains("chr99"), "Exception message should name the missing contig: ${exception.message}")
    }

    @Test
    fun testTooFewGametesInFileThrows() {
        //only 4 gametes are hit by reads, so a 5th parent cannot be chosen
        val filename = writePs4gFile("notEnoughGametes", 4, complementaryParentLines())
        val parentFinder = MostLikelyPs4gParents(Ps4gFileReader(filename), setOf("chr1"))

        val exception = assertFailsWith<IllegalStateException> { parentFinder.bestParents(5) }
        assertEquals("Could only choose 4 best parents.", exception.message)
    }

    @Test
    fun testSingleGameteInFileThrows() {
        //the first parent can be chosen but has no complement, since every read hits it
        val lines = listOf(Ps4gLine(listOf(0), "chr1", 100, 5), Ps4gLine(listOf(0), "chr1", 200, 3))
        val filename = writePs4gFile("oneGamete", 1, lines)
        val parentFinder = MostLikelyPs4gParents(Ps4gFileReader(filename), setOf("chr1"))

        assertFailsWith<IllegalStateException> { parentFinder.bestParents(2) }
    }

    @Test
    fun testGameteCountsFromGameteSets() {
        val parentFinder = parentFinderForHelpers("gameteCounts")
        val gameteSets = listOf(
            Ps4gGameteSet(intArrayOf(0, 1), 5),
            Ps4gGameteSet(intArrayOf(1, 2), 3),
            Ps4gGameteSet(intArrayOf(0), 2)
        )

        //each gamete gets the summed count of every set it appears in, not the number of sets
        assertEquals(mapOf(0 to 7, 1 to 8, 2 to 3), parentFinder.gameteCountsFromGameteSets(gameteSets))
    }

    @Test
    fun testGameteCountsFromEmptyGameteSets() {
        val parentFinder = parentFinderForHelpers("emptyCounts")
        assertEquals(emptyMap(), parentFinder.gameteCountsFromGameteSets(listOf()))
    }

    @Test
    fun testBestParentFromFilteredGameteSets() {
        val parentFinder = parentFinderForHelpers("filteredSets")
        val gameteSets = listOf(
            Ps4gGameteSet(intArrayOf(0, 1), 5),
            Ps4gGameteSet(intArrayOf(2), 3),
            Ps4gGameteSet(intArrayOf(0), 2)
        )

        //excluding gamete 0 discards the sets it appears in, leaving only gamete 2
        assertEquals(2, parentFinder.bestParentFromFilteredGameteSets(gameteSets, 0))

        //excluding gamete 2 leaves gamete 0 with 7 and gamete 1 with 5
        assertEquals(0, parentFinder.bestParentFromFilteredGameteSets(gameteSets, 2))

        //excluding a gamete that appears in no set filters nothing out
        assertEquals(0, parentFinder.bestParentFromFilteredGameteSets(gameteSets, 9))
    }

    @Test
    fun testBestParentFromFilteredGameteSetsReturnsNullWhenAllSetsExcluded() {
        val parentFinder = parentFinderForHelpers("allFiltered")
        val gameteSets = listOf(
            Ps4gGameteSet(intArrayOf(0), 5),
            Ps4gGameteSet(intArrayOf(0, 1), 3)
        )

        assertNull(parentFinder.bestParentFromFilteredGameteSets(gameteSets, 0))
        assertNull(parentFinder.bestParentFromFilteredGameteSets(listOf(), 0))
    }

    @Test
    fun testLikelyParents() {
        //a randomly generated file, seeded so that the test is reproducible, in which gametes 0 and 1 are the parents
        //and gametes 2, 3, and 4 are hit by only some of the reads
        val random = Random(12345)
        val numberOfPositions = 25
        val numberOfReadsAtAPosition = 2
        val probNonParent = 0.5

        val lines = mutableListOf<Ps4gLine>()
        for (pos in 0 until numberOfPositions) {
            for (readNumber in 0 until numberOfReadsAtAPosition) {
                val gameteIndices = mutableListOf<Int>()
                if (random.nextDouble() < 0.5) {
                    gameteIndices.add(0)
                    if (random.nextDouble() < probNonParent) gameteIndices.add(1)
                } else {
                    if (random.nextDouble() < probNonParent) gameteIndices.add(0)
                    gameteIndices.add(1)
                }
                for (nonParent in 2..4) {
                    if (random.nextDouble() < probNonParent) gameteIndices.add(nonParent)
                }
                lines.add(Ps4gLine(gameteIndices, "chr1", pos, 1))
            }
        }

        val filename = writePs4gFile("likely-parent", 5, lines)
        val parentFinder = MostLikelyPs4gParents(Ps4gFileReader(filename), setOf("chr1"))

        assertContentEquals(listOf(0, 1), parentFinder.bestParents(2).sorted())
    }
}
