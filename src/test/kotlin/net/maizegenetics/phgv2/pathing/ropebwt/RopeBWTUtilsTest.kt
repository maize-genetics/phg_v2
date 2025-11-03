package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.RopeBwtIndexTest.Companion.indexFilePrefix
import net.maizegenetics.phgv2.pathing.ropebwt.RopeBwtIndexTest.Companion.inputFasta
import net.maizegenetics.phgv2.pathing.ropebwt.RopeBwtIndexTest.Companion.tempTestDir
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.assertEquals

class RopeBWTUtilsTest {

    companion object {
        @JvmStatic
        @BeforeAll
        fun setup() {
            File(tempTestDir).mkdirs()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(tempTestDir).deleteRecursively()
        }
    }
    //runBuildStep(inputFasta:String, indexFilePrefix:String, numThreads: Int, condaEnvPrefix:String)
    @Test
    fun testRunBuildStep() {
        val numThreads = 3
        RopeBWTUtils.runBuildStep(inputFasta, indexFilePrefix, numThreads, "")

        //verify that the output files exist
        val fmdFile = File("$indexFilePrefix.fmr")
        assert(fmdFile.exists())
    }

    //convertBWTIndex(indexFilePrefix: String, condaEnvPrefix: String)
    //deleteFMRIndex(indexFilePrefix: String)
    @Test
    fun testConvertAndDeleteBWTIndex() {
        val ropebwtIndex = RopeBwtIndex()
        val numThreads = 3

        RopeBWTUtils.runBuildStep(inputFasta, indexFilePrefix, numThreads, "")

        RopeBWTUtils.convertBWTIndex(indexFilePrefix, "")
        val fmdFile = File("$indexFilePrefix.fmd")
        assert(fmdFile.exists())

        RopeBWTUtils.deleteFMRIndex(indexFilePrefix)
        val fmrFile = File("$indexFilePrefix.fmr")
        assert(!fmrFile.exists())
    }

    //buildSuffixArray(indexFilePrefix: String, numThreads: Int, condaEnvPrefix: String)
    @Test
    fun testBuildSuffixArray() {
        val ropebwtIndex = RopeBwtIndex()
        val numThreads = 3

        RopeBWTUtils.runBuildStep(inputFasta, indexFilePrefix, numThreads, "")
        RopeBWTUtils.convertBWTIndex(indexFilePrefix, "")

        RopeBWTUtils.buildSuffixArray(indexFilePrefix, numThreads, "")
        val saFile = File("$indexFilePrefix.fmd.ssa")
        assert(saFile.exists())
    }

    @Test
    fun testParseStringIntoMem() {
        // Test parsing a typical ropebwt3 mem output line
        val testLine = "read1\t10\t50\t2\t40\tcontig1:+:100\tcontig2:-:200"
        val mem = RopeBWTUtils.parseStringIntoMem(testLine)

        assertEquals("read1", mem.readName)
        assertEquals(10, mem.readStart)
        assertEquals(50, mem.readEnd)
        assertEquals(2, mem.numHits)
        assertEquals(2, mem.listMemHits.size)

        // Verify first hit
        assertEquals("contig1", mem.listMemHits[0].contig)
        assertEquals("+", mem.listMemHits[0].strand)
        assertEquals(100, mem.listMemHits[0].pos)

        // Verify second hit
        assertEquals("contig2", mem.listMemHits[1].contig)
        assertEquals("-", mem.listMemHits[1].strand)
        assertEquals(200, mem.listMemHits[1].pos)
    }

    @Test
    fun testParseStringIntoMemWithEmptyHits() {
        // Test parsing a line with empty hits at the end
        val testLine = "read2\t5\t35\t1\t30\tcontig3:+:150\t"
        val mem = RopeBWTUtils.parseStringIntoMem(testLine)

        assertEquals("read2", mem.readName)
        assertEquals(5, mem.readStart)
        assertEquals(35, mem.readEnd)
        assertEquals(1, mem.numHits)
        assertEquals(1, mem.listMemHits.size)

        assertEquals("contig3", mem.listMemHits[0].contig)
        assertEquals("+", mem.listMemHits[0].strand)
        assertEquals(150, mem.listMemHits[0].pos)
    }

    @Test
    fun testParseStringIntoMemMultipleHits() {
        // Test parsing a line with many hits
        val testLine = "read3\t0\t100\t4\t100\thap1:+:50\thap2:+:60\thap3:-:70\thap4:+:80"
        val mem = RopeBWTUtils.parseStringIntoMem(testLine)

        assertEquals("read3", mem.readName)
        assertEquals(0, mem.readStart)
        assertEquals(100, mem.readEnd)
        assertEquals(4, mem.numHits)
        assertEquals(4, mem.listMemHits.size)

        assertEquals("hap1", mem.listMemHits[0].contig)
        assertEquals("hap2", mem.listMemHits[1].contig)
        assertEquals("hap3", mem.listMemHits[2].contig)
        assertEquals("hap4", mem.listMemHits[3].contig)
    }

    @Test
    fun testRunBuildUpdateStep() {
        val numThreads = 3

        // First create initial index
        RopeBWTUtils.runBuildStep(inputFasta, indexFilePrefix, numThreads, "")

        // Verify initial fmr file exists
        val fmrFile = File("$indexFilePrefix.fmr")
        assert(fmrFile.exists())

        // Run the update step (this would normally add sequences to existing index)
        // For testing, we're just verifying it executes without error
        RopeBWTUtils.runBuildUpdateStep(inputFasta, indexFilePrefix, numThreads, "")

        // Verify fmr file still exists after update
        assert(fmrFile.exists())
    }

    @Test
    fun testMEMDataClass() {
        // Test MEM data class creation and properties
        val hit1 = MEMHit("contig1", "+", 100)
        val hit2 = MEMHit("contig2", "-", 200)
        val mem = MEM("testRead", 0, 50, 2, listOf(hit1, hit2))

        assertEquals("testRead", mem.readName)
        assertEquals(0, mem.readStart)
        assertEquals(50, mem.readEnd)
        assertEquals(2, mem.numHits)
        assertEquals(2, mem.listMemHits.size)
        assertEquals("contig1", mem.listMemHits[0].contig)
        assertEquals("+", mem.listMemHits[0].strand)
        assertEquals(100, mem.listMemHits[0].pos)
    }

    @Test
    fun testMEMHitDataClass() {
        // Test MEMHit data class creation
        val hit = MEMHit("chr1", "-", 12345)

        assertEquals("chr1", hit.contig)
        assertEquals("-", hit.strand)
        assertEquals(12345, hit.pos)
    }

    @Test
    fun testRunBuildStepErrorHandling() {
        val numThreads = 3
        // Use an invalid conda environment that will cause conda run to fail
        val invalidCondaEnv = "/this/path/does/not/exist/anywhere"
        val testIndexPrefix = tempTestDir + "testErrorHandling"

        // This should trigger the exit code != 0 error handling path
        var exceptionThrown = false
        try {
            RopeBWTUtils.runBuildStep(inputFasta, testIndexPrefix, numThreads, invalidCondaEnv)
        } catch (e: IllegalStateException) {
            exceptionThrown = true
            // Verify the error message contains information about the failure
            assert(e.message!!.contains("exit code")) {
                "Error message should contain 'exit code', got: ${e.message}"
            }
        }

        assert(exceptionThrown) { "Expected an IllegalStateException to be thrown" }
    }

    @Test
    fun testConvertBWTIndexErrorHandling() {
        // Use an invalid conda environment to trigger failure
        val invalidCondaEnv = "/nonexistent/invalid/conda/path"
        val testIndexPrefix = tempTestDir + "testConvertError"

        // First create a valid .fmr file
        RopeBWTUtils.runBuildStep(inputFasta, testIndexPrefix, 3, "")

        // Now try to convert with invalid conda env, which should trigger error handling
        var exceptionThrown = false
        try {
            RopeBWTUtils.convertBWTIndex(testIndexPrefix, invalidCondaEnv)
        } catch (e: IllegalStateException) {
            exceptionThrown = true
            assert(e.message!!.contains("exit code")) {
                "Error message should contain 'exit code', got: ${e.message}"
            }
        }

        assert(exceptionThrown) { "Expected an IllegalStateException to be thrown" }

        // Cleanup
        File("$testIndexPrefix.fmr").delete()
    }

    @Test
    fun testBuildSuffixArrayErrorHandling() {
        val numThreads = 3
        // Use an invalid conda environment to trigger failure
        val invalidCondaEnv = "/invalid/conda/environment/path"
        val testIndexPrefix = tempTestDir + "testSsaError"

        // First create valid .fmr and .fmd files
        RopeBWTUtils.runBuildStep(inputFasta, testIndexPrefix, numThreads, "")
        RopeBWTUtils.convertBWTIndex(testIndexPrefix, "")

        // Now try to build suffix array with invalid conda env
        var exceptionThrown = false
        try {
            RopeBWTUtils.buildSuffixArray(testIndexPrefix, numThreads, invalidCondaEnv)
        } catch (e: IllegalStateException) {
            exceptionThrown = true
            assert(e.message!!.contains("exit code")) {
                "Error message should contain 'exit code', got: ${e.message}"
            }
        }

        assert(exceptionThrown) { "Expected an IllegalStateException to be thrown" }

        // Cleanup
        File("$testIndexPrefix.fmr").delete()
        File("$testIndexPrefix.fmd").delete()
    }

    @Test
    fun testRunBuildUpdateStepErrorHandling() {
        val numThreads = 3
        // Use an invalid conda environment to trigger failure
        val invalidCondaEnv = "/completely/invalid/path"
        val testIndexPrefix = tempTestDir + "testUpdateError"

        // First create initial index
        RopeBWTUtils.runBuildStep(inputFasta, testIndexPrefix, numThreads, "")

        // Now try update with invalid conda env
        var exceptionThrown = false
        try {
            RopeBWTUtils.runBuildUpdateStep(inputFasta, testIndexPrefix, numThreads, invalidCondaEnv)
        } catch (e: IllegalStateException) {
            exceptionThrown = true
            assert(e.message!!.contains("exit code")) {
                "Error message should contain 'exit code', got: ${e.message}"
            }
        }

        assert(exceptionThrown) { "Expected an IllegalStateException to be thrown" }

        // Cleanup
        File("$testIndexPrefix.fmr").delete()
    }
}