package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith

@ExtendWith(TestExtension::class)
class MapKmersTest {

    companion object {
        @JvmStatic
        @AfterAll
        fun tearDown() {
        }
    }

    @Test
    fun testCliktParams() {
        val mapKmers = MapKmers()

        val resultMissingKmerIndex =
            mapKmers.test("--hvcf-dir ${TestExtension.testVCFDir} --read-files ${TestExtension.testReads} --output-dir ${TestExtension.testOutputDir}")
        assertEquals(resultMissingKmerIndex.statusCode, 1)
        assertEquals(
            "Usage: map-kmers [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --kmer-index: --kmer-index must not be blank\n", resultMissingKmerIndex.output
        )

        val resultMissingReadsAndKeyFile =
            mapKmers.test("--hvcf-dir ${TestExtension.testVCFDir} --kmer-index ${TestExtension.testKmerIndex} --output-dir ${TestExtension.testOutputDir}")
        assertEquals(resultMissingReadsAndKeyFile.statusCode, 1)
        assertEquals(
            "Usage: map-kmers [<options>]\n" +
                    "\n" +
                    "Error: must provide one of --key-file, --read-files\n", resultMissingReadsAndKeyFile.output
        )

        val resultHavingBothReadsAndKeyFile =
            mapKmers.test("--hvcf-dir ${TestExtension.testVCFDir} --kmer-index ${TestExtension.testKmerIndex} --output-dir ${TestExtension.testOutputDir} --key-file ${TestExtension.testReads} --read-files ${TestExtension.testReads}")

        assertEquals(resultHavingBothReadsAndKeyFile.statusCode, 1)
        //This returns the same error message regardless of ordering between key-file and read-files
        assertEquals(
            "Usage: map-kmers [<options>]\n" +
                    "\n" +
                    "Error: option --key-file cannot be used with --read-files\n", resultHavingBothReadsAndKeyFile.output
        )


        val resultMissingOutputDir =
            mapKmers.test("--hvcf-dir ${TestExtension.testVCFDir} --kmer-index ${TestExtension.testKmerIndex} --read-files ${TestExtension.testReads}")
        assertEquals(resultMissingOutputDir.statusCode, 1)
        assertEquals(
            "Usage: map-kmers [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --output-dir: --output-dir/-o must not be blank\n", resultMissingOutputDir.output
        )

        val testMissingHVCFDir = mapKmers.test("--kmer-index ${TestExtension.testKmerIndex} --read-files ${TestExtension.testReads} --output-dir ${TestExtension.testOutputDir}")
        assertEquals(testMissingHVCFDir.statusCode, 1)
        assertEquals(
            "Usage: map-kmers [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --hvcf-dir: --hvcf-dir must not be blank\n", testMissingHVCFDir.output
        )
    }


    @Test
    fun testImportKmerMap() {
//        val kmerIndexFile = "data/test/kmerReadMapping/SimpleIndex.txt"
//        val kmerMapData = importKmerIndex(kmerIndexFile)
//
//        val kmerMap = kmerMapData.kmerHashToLongMap
//
//
//        println(kmerMap.keys)


    }

}