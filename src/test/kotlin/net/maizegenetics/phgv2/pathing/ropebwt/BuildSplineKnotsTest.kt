package net.maizegenetics.phgv2.pathing.ropebwt

import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.BeforeAll
import java.io.File


class BuildSplineKnotsTest {

    companion object {
        val tempTestDir = "${TestExtension.tempDir}BuildSplineKnitsTestDir/"


        //Setup/download  files
        //Resetting on both setup and teardown just to be safe.
        @JvmStatic
        @BeforeAll
        fun setup() {
            resetDirs()
            setupDebugLogging()
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
    fun testCliktParams() {
        val buildSplineKnots = BuildSplineKnots()

        val noVcfDir = buildSplineKnots.test("--output-dir ./knotFiles/")
        assertEquals(1, noVcfDir.statusCode)
        assertEquals("Usage: build-spline-knots [<options>]\n\n" +
                "Error: missing option --vcf-dir\n", noVcfDir.stderr)

        val noOutputFile = buildSplineKnots.test("--vcf-dir testDir")
        assertEquals(1, noOutputFile.statusCode)
        assertEquals("Usage: build-spline-knots [<options>]\n\n" +
                "Error: missing option --output-dir\n", noOutputFile.stderr)
    }

    @Test
    fun testWriteFiles() {
        val buildSplineKnots = BuildSplineKnots()
        val hvcfDir = "data/test/ropebwt/testHVCFs"

        buildSplineKnots.test("--vcf-dir $hvcfDir --output-dir $tempTestDir --disable-spline-downsampling")

        assertEquals(true, File("$tempTestDir/LineA.h_spline_knots.json.gz").exists())
        assertEquals(true, File("$tempTestDir/index_maps.json.gz").exists())
    }

    fun oldSplitter(chrList: String): Set<String> =
        chrList.split(",").map(String::trim).toSet()
    fun newSplitter(chrList: String): Set<String> =
        chrList.split(",").map(String::trim).filter(String::isNotEmpty).toSet()

    @Test
    fun testContigListProcessing() {
        // Test prior set behavior
        val emptyString = ""
        val problematicSet = oldSplitter(emptyString)

        // This set should not be empty, which causes the issue
        assertEquals(1, problematicSet.size)
        assertTrue(problematicSet.contains(""))

        // Test new set creator - verify `""` gets converted to an empty set
        val fixSet = newSplitter(emptyString)
        assertEquals(0, fixSet.size)

        // Test with actual contigs (A)
        val nonEmptyStringA = "chr1,chr2,chr3"
        val nonEmptySetA = newSplitter(nonEmptyStringA)
        assertEquals(3, nonEmptySetA.size)
        assertTrue(nonEmptySetA.contains("chr1"))
        assertTrue(nonEmptySetA.contains("chr2"))
        assertTrue(nonEmptySetA.contains("chr3"))

        // Test with actual contigs (B)
        val nonEmpytStringB = "chr1, chr2, chrA,CHRX,   CHR21A,       htg234a"
        val nonEmptySetB = newSplitter(nonEmpytStringB)
        assertEquals(6, nonEmptySetB.size)
        assertTrue(nonEmptySetB.contains("chr1"))
        assertTrue(nonEmptySetB.contains("chr2"))
        assertTrue(nonEmptySetB.contains("chrA"))
        assertTrue(nonEmptySetB.contains("CHRX"))
        assertTrue(nonEmptySetB.contains("CHR21A"))
        assertTrue(nonEmptySetB.contains("htg234a"))
    }
}