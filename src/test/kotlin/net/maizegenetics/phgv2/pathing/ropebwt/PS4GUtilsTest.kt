package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.Position
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File

class PS4GUtilsTest {

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
    fun testBuildOutputFileName() {
        val inputFile = "/path/to/input/file.txt"
        val outputDir = "/output/dir"
        val result = PS4GUtils.buildOutputFileName(inputFile, outputDir)
        assertEquals("/output/dir/file_ps4g.txt", result)

        val inputFileBed = "/path/to/input/file.bed"
        val resultBed = PS4GUtils.buildOutputFileName(inputFileBed, outputDir)
        assertEquals("/output/dir/file_ps4g.txt", resultBed)

        val withSampleGamete = PS4GUtils.buildOutputFileName(inputFile, outputDir, "sampleA:0")
        assertEquals("/output/dir/file_sampleA:0_ps4g.txt", withSampleGamete)
    }

    @Test
    fun testConvertCountMapToPS4GData() {
        val countMap = mapOf(
            Pair(Position("1", 100), listOf(2, 1, 3)) to 5,
            Pair(Position("1", 50), listOf(1, 2)) to 3
        )
        val ps4gData = PS4GUtils.convertCountMapToPS4GData(countMap, sortPositions = true)

        assertEquals(2, ps4gData.size)
        // Should be sorted by position
        assertEquals(Position("1", 50), ps4gData[0].refPos)
        assertEquals(Position("1", 100), ps4gData[1].refPos)
        // Gamete lists should be sorted within each entry
        assertEquals(listOf(1, 2), ps4gData[0].gameteList)
        assertEquals(listOf(1, 2, 3), ps4gData[1].gameteList)
    }

    @Test
    fun testConvertCountMapToPS4GDataUnsorted() {
        val countMap = linkedMapOf(
            Pair(Position("1", 100), listOf(2, 1, 3)) to 5,
            Pair(Position("1", 50), listOf(1, 2)) to 3
        )
        val ps4gData = PS4GUtils.convertCountMapToPS4GData(countMap, sortPositions = false)

        assertEquals(2, ps4gData.size)
        // Should maintain insertion order
        assertEquals(Position("1", 100), ps4gData[0].refPos)
        assertEquals(Position("1", 50), ps4gData[1].refPos)
        // Gamete lists should NOT be sorted when sortPositions is false
        assertEquals(listOf(2, 1, 3), ps4gData[0].gameteList)
        assertEquals(listOf(1, 2), ps4gData[1].gameteList)
    }

    @Test
    fun testWriteOutPS4GFileIncludesAllGametes() {
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

        val outputFile = "$tempTestDir/test_all_gametes.txt"

        PS4GUtils.writeOutPS4GFile(
            ps4gData,
            sampleGameteCount,
            gameteToIdxMap,
            outputFile,
            listOf("test header"),
            "test command"
        )

        // Read and verify the output
        val lines = File(outputFile).readLines()

        // Find the gamete header lines
        val gameteLines = lines.filter { it.startsWith("#") && it.contains("\t") && !it.startsWith("#PS4G") &&
                                          !it.startsWith("#version") && !it.startsWith("#test") &&
                                          !it.startsWith("#Command") && !it.startsWith("#TotalUniqueCounts") &&
                                          !it.startsWith("#gamete\tgameteIndex") }

        // Should have 5 gamete lines (one for each gamete in gameteToIdxMap)
        assertEquals(5, gameteLines.size, "Should include all gametes from gameteToIdxMap")

        // Check that sampleC:0 with index 4 is present with count 0
        val sampleCLine = gameteLines.find { it.contains("sampleC:0") }
        assertNotNull(sampleCLine, "Should include gamete with zero count")
        assertTrue(sampleCLine!!.contains("\t4\t0"), "sampleC:0 should have index 4 and count 0")
    }

    @Test
    fun testWriteOutPS4GFileGametesSortedByIndex() {
        // Create test data
        val ps4gData = listOf(PS4GData(listOf(0, 1), Position("1", 100), 5))

        // Create gameteToIdxMap with indices NOT in alphabetical order by gamete name
        val gameteToIdxMap = mapOf(
            SampleGamete("sampleC", 0) to 0,
            SampleGamete("sampleA", 0) to 1,
            SampleGamete("sampleB", 1) to 2,
            SampleGamete("sampleB", 0) to 3,
            SampleGamete("sampleA", 1) to 4
        )

        val sampleGameteCount = mapOf(
            SampleGamete("sampleC", 0) to 5,
            SampleGamete("sampleA", 0) to 3,
            SampleGamete("sampleB", 1) to 2,
            SampleGamete("sampleB", 0) to 1,
            SampleGamete("sampleA", 1) to 4
        )

        val outputFile = "$tempTestDir/test_sorted_gametes.txt"

        PS4GUtils.writeOutPS4GFile(
            ps4gData,
            sampleGameteCount,
            gameteToIdxMap,
            outputFile,
            listOf(),
            "test command"
        )

        // Read and verify the output
        val lines = File(outputFile).readLines()

        // Find the gamete header lines
        val gameteLines = lines.filter { it.startsWith("#") && it.contains("\t") &&
                                          !it.startsWith("#PS4G") &&
                                          !it.startsWith("#version") &&
                                          !it.startsWith("#Command") &&
                                          !it.startsWith("#TotalUniqueCounts") &&
                                          !it.startsWith("#gamete\tgameteIndex") }

        assertEquals(5, gameteLines.size)

        // Extract indices from gamete lines
        val indices = gameteLines.map { line ->
            val parts = line.split("\t")
            parts[1].toInt()
        }

        // Verify indices are sorted 0, 1, 2, 3, 4
        assertEquals(listOf(0, 1, 2, 3, 4), indices, "Gametes should be sorted by index")

        // Verify the order of gametes matches index order
        assertTrue(gameteLines[0].contains("sampleC:0\t0"), "Index 0 should be first")
        assertTrue(gameteLines[1].contains("sampleA:0\t1"), "Index 1 should be second")
        assertTrue(gameteLines[2].contains("sampleB:1\t2"), "Index 2 should be third")
        assertTrue(gameteLines[3].contains("sampleB:0\t3"), "Index 3 should be fourth")
        assertTrue(gameteLines[4].contains("sampleA:1\t4"), "Index 4 should be fifth")
    }

    @Test
    fun testWriteOutPS4GFileZeroCountGametes() {
        // Create test data with only gametes 0 and 1 having data
        val ps4gData = listOf(
            PS4GData(listOf(0, 1), Position("1", 100), 10)
        )

        // Create gameteToIdxMap with 4 gametes
        val gameteToIdxMap = mapOf(
            SampleGamete("sampleA", 0) to 0,
            SampleGamete("sampleA", 1) to 1,
            SampleGamete("sampleB", 0) to 2,
            SampleGamete("sampleB", 1) to 3
        )

        // Only gametes 0 and 1 have counts
        val sampleGameteCount = mapOf(
            SampleGamete("sampleA", 0) to 10,
            SampleGamete("sampleA", 1) to 10
        )

        val outputFile = "$tempTestDir/test_zero_counts.txt"

        PS4GUtils.writeOutPS4GFile(
            ps4gData,
            sampleGameteCount,
            gameteToIdxMap,
            outputFile,
            listOf(),
            "test command"
        )

        // Read and verify the output
        val lines = File(outputFile).readLines()

        // Find the gamete header lines
        val gameteLines = lines.filter { it.startsWith("#") && it.contains("\t") &&
                                          !it.startsWith("#PS4G") &&
                                          !it.startsWith("#version") &&
                                          !it.startsWith("#Command") &&
                                          !it.startsWith("#TotalUniqueCounts") &&
                                          !it.startsWith("#gamete\tgameteIndex") }

        assertEquals(4, gameteLines.size)

        // Check that gametes 2 and 3 have count 0
        val sampleB0Line = gameteLines.find { it.contains("sampleB:0") }
        assertNotNull(sampleB0Line)
        assertTrue(sampleB0Line!!.endsWith("\t0"), "sampleB:0 should have count 0")

        val sampleB1Line = gameteLines.find { it.contains("sampleB:1") }
        assertNotNull(sampleB1Line)
        assertTrue(sampleB1Line!!.endsWith("\t0"), "sampleB:1 should have count 0")

        // Check that gametes 0 and 1 have their correct counts
        val sampleA0Line = gameteLines.find { it.contains("sampleA:0") }
        assertNotNull(sampleA0Line)
        assertTrue(sampleA0Line!!.endsWith("\t10"), "sampleA:0 should have count 10")

        val sampleA1Line = gameteLines.find { it.contains("sampleA:1") }
        assertNotNull(sampleA1Line)
        assertTrue(sampleA1Line!!.endsWith("\t10"), "sampleA:1 should have count 10")
    }

    @Test
    fun testWriteOutPS4GFileFormat() {
        // Create simple test data
        val ps4gData = listOf(
            PS4GData(listOf(0, 1), Position("1", 100), 5),
            PS4GData(listOf(2), Position("2", 200), 3)
        )

        val gameteToIdxMap = mapOf(
            SampleGamete("sampleA", 0) to 0,
            SampleGamete("sampleA", 1) to 1,
            SampleGamete("sampleB", 0) to 2
        )

        val sampleGameteCount = mapOf(
            SampleGamete("sampleA", 0) to 5,
            SampleGamete("sampleA", 1) to 5,
            SampleGamete("sampleB", 0) to 3
        )

        val outputFile = "$tempTestDir/test_format.txt"
        val header = listOf("test header line 1", "test header line 2")
        val command = "test command line"

        PS4GUtils.writeOutPS4GFile(
            ps4gData,
            sampleGameteCount,
            gameteToIdxMap,
            outputFile,
            header,
            command
        )

        // Read and verify the output
        val lines = File(outputFile).readLines()

        // Check file format
        assertEquals("#PS4G", lines[0])
        assertEquals("#version=2.0", lines[1])
        assertEquals("#test header line 1", lines[2])
        assertEquals("#test header line 2", lines[3])
        assertEquals("#Command: test command line", lines[4])
        assertTrue(lines[5].startsWith("#TotalUniqueCounts: "))
        assertEquals("#gamete\tgameteIndex\tcount", lines[6])

        // Check data header line
        val dataHeaderIndex = lines.indexOfFirst { it == "gameteSet\trefContig\trefPosBinned\tcount" }
        assertTrue(dataHeaderIndex > 6, "Data header should appear after gamete section")

        // Check data lines
        val dataLines = lines.subList(dataHeaderIndex + 1, lines.size)
        assertEquals(2, dataLines.size)
        assertEquals("0,1\t1\t100\t5", dataLines[0])
        assertEquals("2\t2\t200\t3", dataLines[1])
    }
}
