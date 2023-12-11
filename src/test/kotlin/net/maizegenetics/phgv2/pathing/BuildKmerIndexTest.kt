package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
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

    //Tests to run
    //test buildHaplotypeGraph method
    //   throws error both of --tiledb-path and --hvcf-dir are empty strings
    //   when --hvcf-dir is provided returns not yet implemented
    //
    //test processGraphKmers method
    //test saveKmerHashesAndHapids method
}