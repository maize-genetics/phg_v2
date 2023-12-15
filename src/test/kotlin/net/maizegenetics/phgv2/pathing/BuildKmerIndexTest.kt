package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.AgcCompress
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import org.junit.jupiter.api.extension.ExtendWith
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class BuildKmerIndexTest {

    @Test
    fun testKmerUpdating() {
        val testString = "ACACGTGTAACCGGTTGTGACTGACGGTAACGTCGAATGACGTAACCGTCGA"

        var hashValues = Pair(0L, 0L)
        repeat(32) {
            hashValues = BuildKmerIndex.updateKmerHashAndReverseCompliment(hashValues, testString[it])
        }

        // The first 32-mer of testString
        // A  C  A  C  G  T  G  T  A  A  C  C  G  G  T  T  G  T  G  A  C  T  G  A  C  G  G  T  A  A  C  G
        // 00 01 00 01 10 11 10 11 00 00 01 01 10 10 11 11 10 11 10 00 01 11 10 00 01 10 10 11 00 00 01 10

        // The reverse compliment
        // G  C  A  A  T  G  G  C  A  G  T  C  A  G  T  G  T  T  G  G  C  C  A  A  T  G  T  G  C  A  C  A
        // 01 10 11 11 00 01 01 10 11 01 00 10 11 01 00 01 00 00 01 01 10 10 11 11 00 01 00 01 10 11 10 11

        val expected: Long = 0b0001000110111011000001011010111110111000011110000110101100000110.toLong()
        assertEquals(expected, hashValues.first, "Discrepancy in Kmer hash")

        val expectedRC: Long = 0b0110111100010110110100101101000100000101101011110001000110111011.toLong()
        assertEquals(expectedRC, hashValues.second, "Discrepancy in reverse compliment Kmer hash")

        //test trying to convert an invalid character
        val badTestString = "ACACGTGTAACCGGTTGTMACTGACGGTAACGTCGAATGACGTAACCGTCGA"
        var badHashValues = Pair(0L, 0L)

        val exception = assertThrows<IllegalArgumentException> {
            repeat(32) {
                badHashValues = BuildKmerIndex.updateKmerHashAndReverseCompliment(badHashValues, badTestString[it])
            }
        }

        assertEquals("Attempted to update kmer hash with an invalid nucleotide character(M). Must be one of A,G,C,T", exception.message)

    }

    @Test
    fun testCliktParams() {
        //parameters:
        //tiledbPath - defaults to ""
        //agcPath - required
        //maxHaplotypeProportion - defaults to 0.75
        //hashMask - default to 3L
        //hashFilterValue - defauls to 1L
        //hvcfDir - defaults to ""

        val testBuild = BuildKmerIndex()
        val noargResult = testBuild.test("")
        assertEquals(1, noargResult.statusCode)
        assertEquals("Usage: build-kmer-index [<options>]\n\n" +
                "Error: missing option --agc-path\n", noargResult.stderr)

    }


    @Test
    fun testProcessGraphKmers() {

        //populate the AGC database
        setupAgc()

        //set up temporary file names
        val tempTestDir = "${TestExtension.tempDir}kmerTest/"
        val tempHvcfDir = "${tempTestDir}hvcfDir/"
        val tempAGCDir = "${TestExtension.testOutputFastaDir}/dbPath"

        //copy hvcf files to temp directory
        listOf(TestExtension.smallseqLineAHvcfFile,TestExtension.smallseqLineBHvcfFile).forEach { hvcfFile ->
            val dst = File("$tempHvcfDir${File(hvcfFile).name}")
            if (!dst.exists()) {
                File(hvcfFile).copyTo(dst)
            }
        }

        //create a HaplotypeGraph from the hvcf files
        val buildIndexResult = BuildKmerIndex().test("--agc-path $tempAGCDir --hvcf-dir $tempHvcfDir")

        //Was the index created?
        assertEquals(0, buildIndexResult.statusCode)
        assert(File("${tempHvcfDir}/kmerIndex.txt").exists())
    }

    @Test
    fun testTiledb() {
        //Setting the tiledb path but not the hvcf should generate a not implemented error
        val tempTestDir = "${TestExtension.tempDir}kmerTest/"
        val tempAGCDir = "${TestExtension.testOutputFastaDir}/dbPath"

        //try to build a graph from a (non-existent) tiledb database
        val exception = assertThrows<NotImplementedError> {
            BuildKmerIndex().test("--agc-path $tempAGCDir --tiledb-path ${TestExtension.testTileDBURI}")
        }
        assertEquals("An operation is not implemented: TileDB VCF Reader Not implemented yet.  Please run with --hvcf-dir", exception.message)

    }


    private fun setupAgc() {
        //create an AGC record with the Ref in it
        val altFileListFile = TestExtension.testOutputFastaDir+"/agc_altList.txt"
        BufferedWriter(FileWriter(altFileListFile)).use { writer ->
            writer.write("data/test/smallseq/LineA.fa\n")
            writer.write("data/test/smallseq/LineB.fa\n")
        }

        val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"
        File(dbPath).mkdirs()

        //Call AGCCompress to create the AGC file
        val agcCompress = AgcCompress()
        agcCompress.processAGCFiles(dbPath,altFileListFile,"data/test/smallseq/Ref.fa")
    }
}