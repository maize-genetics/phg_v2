package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.utils.TileDBCoreVariantQueries
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
class LoadHvcfTest {
    companion object {

        private val multiInputDir = "${TestExtension.tempDir}/multi-input/"
        private val outputHvcfDir = "${TestExtension.tempDir}/output/"
        val lineAhvcf = "data/test/tiledbCoreHvcf/LineA.h.vcf"
        val lineBhvcf = "data/test/tiledbCoreHvcf/LineB.h.vcf"
        val imputedHvcf = "data/test/resequenceHaplotypeVCF/Imputation.h.vcf"
        val dbPath = TestExtension.testTileDBURI
        val altHeaderArray = dbPath + "/alt_header_array"
        val variantsArray = dbPath + "/hvcf_variants_array"

        @BeforeAll
        @JvmStatic
        fun setup() {
            File(multiInputDir).mkdirs()
            File(outputHvcfDir).mkdirs()
            File(dbPath).mkdirs()

            File(lineAhvcf).copyTo(File(multiInputDir + File(lineAhvcf).name))
            File(lineBhvcf).copyTo(File(multiInputDir + File(lineBhvcf).name))
            File(imputedHvcf).copyTo(File(multiInputDir + File(imputedHvcf).name))
        }

        @AfterAll
        @JvmStatic
        fun teardown() {
            File(outputHvcfDir).deleteRecursively()
            File(multiInputDir).deleteRecursively()
            // delete the tempDir
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testCliktParams() {
        val loadHVCF = LoadHvcf()

        // Test missing vcf-dir parameter
        val resultMissingHVCFDir = loadHVCF.test("--db-path ${TestExtension.testTileDBURI} ")
        assertEquals(resultMissingHVCFDir.statusCode, 1)
        assertEquals("Usage: load-hvcf [<options>]\n" +
                "\n" +
                "Error: missing option --hvcf-dir\n",resultMissingHVCFDir.output)

    }

    @Test
    fun testLoadHvcf() {
        // First create the arrays:
        val initHvcfArray = InitHvcfArray()
        val resultInit = initHvcfArray.test("--db-path ${TestExtension.testTileDBURI}")
        assertEquals(resultInit.statusCode, 0)

        // Now load the hvcf files
        val loadHVCF = LoadHvcf()
        val resultLoad = loadHVCF.test("--hvcf-dir $multiInputDir --db-path ${TestExtension.testTileDBURI}")
        assertEquals(resultLoad.statusCode, 0)

        // Add some asserts to check the output
        // Get number of distinct sample names from the variantsArray
        val names = TileDBCoreVariantQueries.queryDistinctSampleNames(variantsArray)
        println("Distinct sample names: $names")
        assertEquals(3,names.size)
        assertTrue(names.contains("LineA"))
        assertTrue(names.contains("LineB"))
        assertTrue(names.contains("TestLine2"))

        // get number of distinct samples from the altHeaderArray
        // TestLine2 will not show up as a sample name in the altHeaderArray
        val names2 = TileDBCoreVariantQueries.queryDistinctSampleNames(altHeaderArray)
        println("Distinct sample names: $names2")
        assertEquals(2,names2.size)
        assertTrue(names2.contains("LineA"))
        assertTrue(names2.contains("LineB"))

        // Get number of distinct refRanges from the variantsArray
        val refRanges = TileDBCoreVariantQueries.queryDistinctRefRanges(variantsArray)
        assertEquals(38,refRanges.size)

        // GEt number of distinct refRanges from the altHeaderArray
        val refRangeAltHeaders = TileDBCoreVariantQueries.queryDistinctRefRanges(altHeaderArray)
        assertEquals(38,refRangeAltHeaders.size)

        // verify the number of refRanges in the altHeaderArray is the same as in the variantsArray
        assertEquals(refRanges.size,refRangeAltHeaders.size)

    }

}