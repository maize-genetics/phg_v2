package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.Position
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertFalse
import org.junit.jupiter.api.Assertions.assertNotEquals
import org.junit.jupiter.api.Assertions.assertNull
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.assertContentEquals

class Ps4gFileReaderTest {

    companion object {
        val tempTestDir = "${TestExtension.tempDir}PS4GUtilsTestDir/"

        @JvmStatic
        @BeforeAll
        fun setup() {
            resetDirs()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            resetDirs()
        }

        private fun resetDirs() {
            File(TestExtension.tempDir).deleteRecursively()
            File(tempTestDir).deleteRecursively()

            File(TestExtension.tempDir).mkdirs()
            File(tempTestDir).mkdirs()
        }

        /** The gamete index map shared by the multi-contig tests. */
        private val multiContigGametes = mapOf(
            SampleGamete("sampleA", 0) to 0,
            SampleGamete("sampleA", 1) to 1,
            SampleGamete("sampleB", 0) to 2,
            SampleGamete("sampleB", 1) to 3
        )

        /**
         * Writes a PS4G file with mappings on chr1, chr2, and scaf1, and returns its path.
         * Each contig gets two binned positions.
         */
        private fun createMultiContigPs4gFile(fileName: String): String {
            val ps4gData = listOf(
                PS4GData(listOf(0, 2), Position("chr1", 100), 5),
                PS4GData(listOf(1, 3), Position("chr1", 200), 3),
                PS4GData(listOf(0), Position("chr2", 100), 7),
                PS4GData(listOf(2), Position("chr2", 300), 2),
                PS4GData(listOf(3), Position("scaf1", 100), 1),
                PS4GData(listOf(3), Position("scaf1", 200), 4)
            )

            val gameteCounts = multiContigGametes.keys.associateWith { sampleGamete ->
                val index = multiContigGametes[sampleGamete]!!
                ps4gData.filter { it.gameteList.contains(index) }.sumOf { it.count }
            }

            val outputFile = "$tempTestDir/$fileName"
            PS4GUtils.writeOutPS4GFile(
                ps4gData,
                gameteCounts,
                multiContigGametes,
                outputFile,
                listOf("test header"),
                "test command"
            )
            return outputFile
        }
    }

    @Test
    fun testReadPS4GFile() {
        // Create test data
        val ps4gData = listOf(
            PS4GData(listOf(0, 2), Position("1", 100), 5),
            PS4GData(listOf(1, 3), Position("1", 200), 3)
        )

        // Create gameteToIdxMap with 5 gametes (0-4)
        val gameteToIdxMap = mapOf(
            SampleGamete("sampleA", 0) to 0,
            SampleGamete("sampleA", 1) to 1,
            SampleGamete("sampleB", 0) to 2,
            SampleGamete("sampleB", 1) to 3,
            SampleGamete("sampleC", 0) to 4  // This gamete has no data
        )

        // Only gametes 0-3 have counts
        val sampleGameteCount = mapOf(
            SampleGamete("sampleA", 0) to 5,
            SampleGamete("sampleA", 1) to 3,
            SampleGamete("sampleB", 0) to 5,
            SampleGamete("sampleB", 1) to 3
        )

        val outputFile = "${PS4GUtilsTest.tempTestDir}/test_all_gametes.txt"

        PS4GUtils.writeOutPS4GFile(
            ps4gData,
            sampleGameteCount,
            gameteToIdxMap,
            outputFile,
            listOf("test header"),
            "test command"
        )

        // Read and verify the output
        val ps4gFileReader = Ps4gFileReader(outputFile)
        val indexToGameteMap = ps4gFileReader.gameteIndexMap()
        assertEquals(5, indexToGameteMap.size, "Test of indexToGameteMap size failed.")
        assertEquals("sampleA:0", indexToGameteMap[0])
        assertEquals("sampleA:1", indexToGameteMap[1])
        assertEquals("sampleB:0", indexToGameteMap[2])

        val contigs = ps4gFileReader.contigSet()
        assertEquals(1, contigs.size)
        assertEquals("1", contigs.first())

        val reads = ps4gFileReader.readMapForContig("1")
        assertEquals(2, reads?.size)
        assertContentEquals(intArrayOf(0,2), reads!![100]!!.first().gameteIndices)
        assertContentEquals(intArrayOf(1,3), reads!![200]!!.first().gameteIndices)
        assertEquals(5, reads!![100]!!.first().count)
        assertEquals(3, reads!![200]!!.first().count)

    }

    @Test
    fun testNoContigFilter() {
        // The contigSet parameter defaults to an empty set, which means "keep every contig".
        // Passing an empty set explicitly must behave the same way.
        val outputFile = createMultiContigPs4gFile("noFilter.ps4g")

        val defaultReader = Ps4gFileReader(outputFile)
        assertEquals(setOf("chr1", "chr2", "scaf1"), defaultReader.contigSet())

        val emptySetReader = Ps4gFileReader(outputFile, emptySet())
        assertEquals(setOf("chr1", "chr2", "scaf1"), emptySetReader.contigSet())

        // Contigs that used to be dropped by the "scaf" name check are now retained.
        assertEquals(2, emptySetReader.readMapForContig("scaf1")?.size)
    }

    @Test
    fun testContigFilterKeepsOnlyRequestedContigs() {
        val outputFile = createMultiContigPs4gFile("filtered.ps4g")

        val reader = Ps4gFileReader(outputFile, setOf("chr1", "chr2"))
        assertEquals(setOf("chr1", "chr2"), reader.contigSet(), "scaf1 should have been filtered out")

        // Filtered contigs have no read map at all, rather than an empty one.
        assertNull(reader.readMapForContig("scaf1"), "Filtered contig should not be in the data map")

        // Data for the retained contigs is unchanged by filtering.
        val chr1Reads = reader.readMapForContig("chr1")
        assertEquals(2, chr1Reads?.size)
        assertContentEquals(intArrayOf(0, 2), chr1Reads!![100]!!.first().gameteIndices)
        assertEquals(5, chr1Reads[100]!!.first().count)
        assertContentEquals(intArrayOf(1, 3), chr1Reads[200]!!.first().gameteIndices)
        assertEquals(3, chr1Reads[200]!!.first().count)

        val chr2Reads = reader.readMapForContig("chr2")
        assertEquals(2, chr2Reads?.size)
        assertContentEquals(intArrayOf(0), chr2Reads!![100]!!.first().gameteIndices)
        assertEquals(7, chr2Reads[100]!!.first().count)

        // The gamete index block is read before the body, so filtering must not change it.
        assertEquals(4, reader.gameteIndexMap().size)
        assertEquals("sampleB:1", reader.gameteIndexMap()[3])
    }

    @Test
    fun testContigFilterWithSingleContig() {
        val outputFile = createMultiContigPs4gFile("singleContigFilter.ps4g")

        val reader = Ps4gFileReader(outputFile, setOf("scaf1"))
        assertEquals(setOf("scaf1"), reader.contigSet())
        assertNull(reader.readMapForContig("chr1"))
        assertNull(reader.readMapForContig("chr2"))
        assertEquals(2, reader.readMapForContig("scaf1")?.size)
    }

    @Test
    fun testContigFilterWithNamesNotInFile() {
        val outputFile = createMultiContigPs4gFile("unknownContigFilter.ps4g")

        // Names in the contig set that are not in the file are simply absent from the result;
        // they are not an error.
        val reader = Ps4gFileReader(outputFile, setOf("chr1", "chr99"))
        assertEquals(setOf("chr1"), reader.contigSet())
        assertNull(reader.readMapForContig("chr99"))

        // A contig set that matches nothing yields an empty reader rather than a failure.
        val emptyReader = Ps4gFileReader(outputFile, setOf("chr99"))
        assertTrue(emptyReader.contigSet().isEmpty(), "Expected no contigs to survive filtering")
        assertEquals(4, emptyReader.gameteIndexMap().size, "Gamete index map is read before filtering")
    }

    @Test
    fun testPs4gGameteSetEquality() {
        // Ps4gGameteSet holds an IntArray, so it overrides equals/hashCode to compare contents.
        // Without those overrides two sets with equal indices would compare unequal.
        val first = Ps4gFileReader.Ps4gGameteSet(intArrayOf(0, 2), 5)
        val second = Ps4gFileReader.Ps4gGameteSet(intArrayOf(0, 2), 5)

        assertEquals(first, second)
        assertEquals(first, first)
        assertEquals(first.hashCode(), second.hashCode())

        // Equal sets must collapse in a hash-based collection.
        assertEquals(1, setOf(first, second).size)

        // Differing counts, differing indices, and differing index order are all unequal.
        assertNotEquals(first, Ps4gFileReader.Ps4gGameteSet(intArrayOf(0, 2), 6))
        assertNotEquals(first, Ps4gFileReader.Ps4gGameteSet(intArrayOf(0, 3), 5))
        assertNotEquals(first, Ps4gFileReader.Ps4gGameteSet(intArrayOf(2, 0), 5))
        assertNotEquals(first, Ps4gFileReader.Ps4gGameteSet(intArrayOf(0), 5))

        assertFalse(first.equals(null))
        assertFalse(first.equals("not a gamete set"))
    }
}