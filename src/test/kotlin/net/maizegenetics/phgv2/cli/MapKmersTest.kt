package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.utils.importKmerIndex
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File

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
            mapKmers.test("--read-files ${TestExtension.testReads} --output-dir ${TestExtension.testOutputDir}")
        assertEquals(resultMissingKmerIndex.statusCode, 1)
        assertEquals(
            "Usage: map-kmers [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --kmer-index: --kmer-index must not be blank\n", resultMissingKmerIndex.output
        )

        val resultMissingReads =
            mapKmers.test("--kmer-index ${TestExtension.testKmerIndex} --output-dir ${TestExtension.testOutputDir}")
        assertEquals(resultMissingReads.statusCode, 1)
        assertEquals(
            "Usage: map-kmers [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --read-files: --read-files must not be blank\n", resultMissingReads.output
        )

        val resultMissingOutputDir =
            mapKmers.test("--kmer-index ${TestExtension.testKmerIndex} --read-files ${TestExtension.testReads}")
        assertEquals(resultMissingOutputDir.statusCode, 1)
        assertEquals(
            "Usage: map-kmers [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --output-dir: --output-dir/-o must not be blank\n", resultMissingOutputDir.output
        )
    }


    @Test
    fun testImportKmerMap() {
        val kmerIndexFile = "data/test/kmerReadMapping/SimpleIndex.txt"
        val kmerMapData = importKmerIndex(kmerIndexFile)

        val kmerMap = kmerMapData.kmerHashToLongMap


        println(kmerMap.keys)


    }

}