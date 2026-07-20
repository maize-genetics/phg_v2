package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.Position
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertEquals
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
}