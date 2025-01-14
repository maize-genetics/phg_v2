package net.maizegenetics.phgv2.utils

import biokotlin.util.bufferedReader
import htsjdk.variant.vcf.VCFFileReader
import io.tiledb.java.api.Array
import io.tiledb.java.api.QueryType
import io.tiledb.java.api.*
import net.maizegenetics.phgv2.cli.TestExtension
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File

@ExtendWith(TestExtension::class)
class TileDBCoreHvcfUtilsTest {

    private val myLogger = LogManager.getLogger(TileDBCoreHvcfUtilsTest::class.java)

    companion object {

        private val multiInputDir = "${TestExtension.tempDir}/multi-input/"
        private val outputHvcfDir = "${TestExtension.tempDir}/output/"
        private val imputedHvcfDir = "${TestExtension.tempDir}/imputed/"
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
            File(imputedHvcfDir).mkdirs()

            File(lineAhvcf).copyTo(File(multiInputDir + File(lineAhvcf).name))
            File(lineBhvcf).copyTo(File(multiInputDir + File(lineBhvcf).name))
            File(imputedHvcf).copyTo(File(imputedHvcfDir + File(imputedHvcf).name))
        }

        @AfterAll
        @JvmStatic
        fun teardown() {
            File(outputHvcfDir).deleteRecursively()
            File(multiInputDir).deleteRecursively()
            File(imputedHvcfDir).deleteRecursively()
            // delete the tempDir
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testTiledbArrayExists() {
        // This is mostly to  verify the commands that test for
        // tiledb array existence
        // make sure tiledb main folder exists
        File(dbPath).mkdirs()
        println("Creating tiledb array")
        TileDBCoreHvcfUtils.createTileDBCoreArrays(dbPath)
        val array = Array(Context(), altHeaderArray, QueryType.TILEDB_READ)
        assertTrue(array.schema != null)
        array.close()

        // verify  the variants array exists
        val variantsArray = Array(Context(), variantsArray, QueryType.TILEDB_READ)
        assertTrue(variantsArray.schema != null)
        variantsArray.close()

        // Verify an exception is thrown if name exists but is not a tiledb array
        assertThrows<TileDBError>{
            Array(Context(), lineAhvcf, QueryType.TILEDB_READ)
        }

        // Verify exception is thrown if the array does not exist
        assertThrows<TileDBError>{
            Array(Context(), "badArray", QueryType.TILEDB_READ)
        }

        // delete the tiledbArray so next tests can recreate what they need
        File(dbPath).deleteRecursively()

    }

    @Test
    fun testParseTIledbVariantData() {
        val vcfReader = VCFFileReader(File(lineAhvcf), false)
        val variantData = TileDBCoreHvcfUtils.parseTiledbVariantData(vcfReader)
        println("Finished parsing lineAhvcf variant data")
        assertEquals(38, variantData.size)
        // Verify the first and lsat entries in the variant data
        // they should be:
        //    RefRange=1:1-1000, ID1=12f0cec9102e84a161866e37072443b7, SampleName=LineA, ID2=12f0cec9102e84a161866e37072443b7
        //    RefRange=2:49501-50500, ID1=0eb9029f3896313aebc69c8489923141, SampleName=LineA, ID2=0eb9029f3896313aebc69c8489923141

        val firstEntry = variantData.first()
        assertEquals(4, firstEntry.size)
        assertEquals("1:1-1000", firstEntry["RefRange"])
        assertEquals("12f0cec9102e84a161866e37072443b7", firstEntry["ID1"])
        assertEquals("LineA", firstEntry["SampleName"])
        assertEquals("12f0cec9102e84a161866e37072443b7", firstEntry["ID2"])

        val lastEntry = variantData.last()
        assertEquals(4, lastEntry.size)
        assertEquals("2:49501-50500", lastEntry["RefRange"])
        assertEquals("0eb9029f3896313aebc69c8489923141", lastEntry["ID1"])
        assertEquals("LineA", lastEntry["SampleName"])
        assertEquals("0eb9029f3896313aebc69c8489923141", lastEntry["ID2"])

    }

    @Test
    fun testParseAltHeadersTiledb() {
        // testing output from parseTiledbAltHeaders
        val vcfReader = VCFFileReader(File(lineAhvcf), false)
        val altHeaders = TileDBCoreHvcfUtils.parseTiledbAltHeaders(vcfReader)
        println("Finished parsing smallSeqLineAHvcfFile")

        // Read the altheader data into a list to compare to what
        // is stored from the parseTiledbAltHeaders function
        // Use buffered reader to read the lines from lineAhvcf, but keep only the lines
        // that start with "##ALT"
        val altHeaderListFromHvcf = bufferedReader(lineAhvcf).readLines().filter { it.startsWith("##ALT") }
        assertTrue(altHeaderListFromHvcf.size == 38)
        assertEquals(altHeaders.size, altHeaderListFromHvcf.size)

        // Extract data
        // Note the ALT headers are stored in the h.vcf file ordered lexicographically by checksum (Hapid)
        val ids = altHeaders.map { it["ID"].orEmpty() }
        val sampleNames = altHeaders.map { it["SampleName"].orEmpty() }
        val regions = altHeaders.map { it["Regions"].orEmpty() }
        val refChecksums = altHeaders.map { it["RefChecksum"].orEmpty() }
        val refRanges = altHeaders.map { it["RefRange"].orEmpty() }

        // verify the number of entries from parseTiledbAltHeaders matches the number of ALT header lines in the VCF file
        assertEquals(altHeaders.size, altHeaderListFromHvcf.size)

        // THe first ALT header line has 1 region in it: 2:49501-50300, verify it is stored correctly
        assertEquals("2:49501-50300", regions[0])
        // The second ALT header line has 2 regions in it: Regions="1:1-950,1:952-1000"
        // Verify that the data is extracted correctly
        assertEquals("1:1-950,1:952-1000", regions[1])
        // Verify the sample name is LineA
        assertEquals("LineA", sampleNames[1])
        // Verify the hapID of the second entry is 12f0cec9102e84a161866e37072443b7
        assertEquals("12f0cec9102e84a161866e37072443b7", ids[1])
        // THe last line in the file is this: verify the values in the parsed header data
        //##ALT=<ID=f50fe6d6b3d9a9d305889db977969916,Description="haplotype data for line: LineA",
        // Source="/Users/lcj34/temp/phgv2Tests/tempDir/testTileDBURI//assemblies.agc",SampleName=LineA,
        // Regions=1:11001-12000,Checksum=f50fe6d6b3d9a9d305889db977969916,
        // RefChecksum=2b4590f722ef9229c15d29e0b4e51a0e,RefRange=1:11001-12000>

        assertEquals("f50fe6d6b3d9a9d305889db977969916", ids.last())
        assertEquals("LineA", sampleNames.last())
        assertEquals("1:11001-12000", regions.last())
        assertEquals("2b4590f722ef9229c15d29e0b4e51a0e", refChecksums.last())
        assertEquals("1:11001-12000", refRanges.last())

        println("Finished extracting data")
    }

}