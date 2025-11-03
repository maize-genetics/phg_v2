package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class AlignReadsTest {

    companion object {
        val tempTestDir = "${TestExtension.tempDir}alignReadsTest/"
        val indexFilePrefix = tempTestDir + "testIndex"
        val inputFasta = TestExtension.smallSeqInputDir + "/pangenome/pangenome.fa"
        val testQueryReads = "data/test/kmerReadMapping/simulatedReads/LineA_1.fq"
        val outputBed = tempTestDir + "aligned_reads.bed"

        @JvmStatic
        @BeforeAll
        fun setup() {
            setupDebugLogging()
            File(tempTestDir).mkdirs()

            // Create the index files needed for testing
            val numThreads = 3
            val ropeBwtIndex = RopeBwtIndex()
            RopeBWTUtils.runBuildStep(inputFasta, indexFilePrefix, numThreads, "")
            RopeBWTUtils.convertBWTIndex(indexFilePrefix, "")
            RopeBWTUtils.buildSuffixArray(indexFilePrefix, numThreads, "")
            // Build the chromosome length file needed for -p parameter
            ropeBwtIndex.buildChrLengthFile(inputFasta, indexFilePrefix)
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(tempTestDir).deleteRecursively()
        }
    }

    @Test
    fun testCliktParamsMissingIndex() {
        val alignReads = AlignReads()
        val result = alignReads.test("--query $testQueryReads --output $outputBed")

        assertEquals(1, result.statusCode)
        assertTrue(result.stderr.contains("missing option --index"))
    }

    @Test
    fun testCliktParamsMissingQuery() {
        val alignReads = AlignReads()
        val result = alignReads.test("--index $indexFilePrefix.fmd --output $outputBed")

        assertEquals(1, result.statusCode)
        assertTrue(result.stderr.contains("missing option --query"))
    }

    @Test
    fun testCliktParamsMissingOutput() {
        val alignReads = AlignReads()
        val result = alignReads.test("--index $indexFilePrefix.fmd --query $testQueryReads")

        assertEquals(1, result.statusCode)
        assertTrue(result.stderr.contains("missing option --output"))
    }

    @Test
    fun testValidateInputsNonExistentIndex() {
        val alignReads = AlignReads()
        val nonExistentIndex = tempTestDir + "nonexistent.fmd"
        val testOutput = tempTestDir + "test_output.bed"

        try {
            val result = alignReads.test("--index $nonExistentIndex --query $testQueryReads --output $testOutput")
            // Should fail before reaching here
            assertEquals(1, result.statusCode)
        } catch (e: IllegalArgumentException) {
            assertTrue(e.message!!.contains("Index file does not exist"))
        }
    }

    @Test
    fun testValidateInputsNonExistentQuery() {
        val alignReads = AlignReads()
        val nonExistentQuery = tempTestDir + "nonexistent.fq"
        val testOutput = tempTestDir + "test_output.bed"

        try {
            val result = alignReads.test("--index $indexFilePrefix.fmd --query $nonExistentQuery --output $testOutput")
            // Should fail before reaching here
            assertEquals(1, result.statusCode)
        } catch (e: IllegalArgumentException) {
            assertTrue(e.message!!.contains("Query file does not exist"))
        }
    }

    @Test
    fun testAlignReadsBasic() {
        val alignReads = AlignReads()
        val testOutput = tempTestDir + "test_basic_output.bed"

        val result = alignReads.test("--index $indexFilePrefix.fmd --query $testQueryReads --output $testOutput")

        assertEquals(0, result.statusCode)

        // Verify output file was created
        val outputFile = File(testOutput)
        assertTrue(outputFile.exists(), "Output BED file should exist")
        assertTrue(outputFile.length() > 0, "Output BED file should not be empty")
    }

    @Test
    fun testAlignReadsWithCustomMinSmemLen() {
        val alignReads = AlignReads()
        val testOutput = tempTestDir + "test_custom_smem.bed"

        val result = alignReads.test(
            "--index $indexFilePrefix.fmd --query $testQueryReads --output $testOutput --min-smem-len 20"
        )

        assertEquals(0, result.statusCode)

        // Verify output file was created
        val outputFile = File(testOutput)
        assertTrue(outputFile.exists(), "Output BED file should exist")
        assertTrue(outputFile.length() > 0, "Output BED file should not be empty")
    }

    @Test
    fun testAlignReadsWithMultipleThreads() {
        val alignReads = AlignReads()
        val testOutput = tempTestDir + "test_multithread.bed"

        val result = alignReads.test(
            "--index $indexFilePrefix.fmd --query $testQueryReads --output $testOutput --threads 2"
        )

        assertEquals(0, result.statusCode)

        // Verify output file was created
        val outputFile = File(testOutput)
        assertTrue(outputFile.exists(), "Output BED file should exist")
        assertTrue(outputFile.length() > 0, "Output BED file should not be empty")
    }

    @Test
    fun testAlignReadsAllParameters() {
        val alignReads = AlignReads()
        val testOutput = tempTestDir + "test_all_params.bed"

        val result = alignReads.test(
            "--index $indexFilePrefix.fmd --query $testQueryReads --output $testOutput " +
            "--min-smem-len 25 --threads 3"
        )

        assertEquals(0, result.statusCode)

        // Verify output file was created
        val outputFile = File(testOutput)
        assertTrue(outputFile.exists(), "Output BED file should exist")
        assertTrue(outputFile.length() > 0, "Output BED file should not be empty")
    }

    @Test
    fun testAlignReadsOutputDirectoryCreation() {
        val alignReads = AlignReads()
        val nestedOutput = tempTestDir + "nested/directory/output.bed"

        // Delete the nested directory if it exists to ensure clean test
        File(tempTestDir + "nested").deleteRecursively()

        val result = alignReads.test(
            "--index $indexFilePrefix.fmd --query $testQueryReads --output $nestedOutput"
        )

        assertEquals(0, result.statusCode)

        // Verify output directory was created
        val outputFile = File(nestedOutput)
        assertTrue(outputFile.parentFile.exists(), "Output directory should be created")
        assertTrue(outputFile.exists(), "Output BED file should exist")
    }

    @Test
    fun testAlignReadsDifferentQueryFile() {
        val alignReads = AlignReads()
        val testOutput = tempTestDir + "test_lineB.bed"
        val lineBQuery = "data/test/kmerReadMapping/simulatedReads/LineB_1.fq"

        val result = alignReads.test(
            "--index $indexFilePrefix.fmd --query $lineBQuery --output $testOutput"
        )

        assertEquals(0, result.statusCode)

        // Verify output file was created
        val outputFile = File(testOutput)
        assertTrue(outputFile.exists(), "Output BED file should exist")
        assertTrue(outputFile.length() > 0, "Output BED file should not be empty")
    }

    @Test
    fun testAlignReadsVerifyBedFormat() {
        val alignReads = AlignReads()
        val testOutput = tempTestDir + "test_bed_format.bed"

        val result = alignReads.test(
            "--index $indexFilePrefix.fmd --query $testQueryReads --output $testOutput"
        )

        assertEquals(0, result.statusCode)

        // Read output and verify it has expected BED-like format
        val outputFile = File(testOutput)
        assertTrue(outputFile.exists(), "Output BED file should exist")

        val lines = outputFile.readLines()
        assertTrue(lines.isNotEmpty(), "Output should contain alignment results")

        // BED format from ropebwt3 mem should have tab-separated fields
        // Format: readName \t readStart \t readEnd \t numHits \t strand \t contig:strand:pos ...
        val firstLine = lines.first()
        val fields = firstLine.split("\t")
        assertTrue(fields.size >= 4, "Each line should have at least 4 tab-separated fields")
    }

    @Test
    fun testAlignReadsWithSubsetPositions() {
        val alignReads = AlignReads()
        val testOutput = tempTestDir + "test_subset_positions.bed"

        val result = alignReads.test(
            "--index $indexFilePrefix.fmd --query $testQueryReads --output $testOutput --subset-positions 25"
        )

        assertEquals(0, result.statusCode)

        // Verify output file was created
        val outputFile = File(testOutput)
        assertTrue(outputFile.exists(), "Output BED file should exist")
        assertTrue(outputFile.length() > 0, "Output BED file should not be empty")
    }

    @Test
    fun testAlignReadsWithGapReporting() {
        val alignReads = AlignReads()
        val testOutput = tempTestDir + "test_gaps.bed"

        val result = alignReads.test(
            "--index $indexFilePrefix.fmd --query $testQueryReads --output $testOutput --gap 50"
        )

        assertEquals(0, result.statusCode)

        // Verify output file was created
        val outputFile = File(testOutput)
        assertTrue(outputFile.exists(), "Output BED file should exist")
        // Note: Gap output might be empty if all regions are covered by SMEMs
    }

    @Test
    fun testAlignReadsWithAllNewParameters() {
        val alignReads = AlignReads()
        val testOutput = tempTestDir + "test_all_new_params.bed"

        val result = alignReads.test(
            "--index $indexFilePrefix.fmd --query $testQueryReads --output $testOutput " +
            "--subset-positions 30 --gap 40 --threads 2 --min-smem-len 25"
        )

        assertEquals(0, result.statusCode)

        // Verify output file was created
        val outputFile = File(testOutput)
        assertTrue(outputFile.exists(), "Output BED file should exist")
    }
}
