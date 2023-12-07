package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.utils.getChecksum
import net.maizegenetics.phgv2.utils.testMergingMAF
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
class AlignAssembliesTest {

    @Test
    fun testNumThreadsAndRuns() {
        // ok - this is now working.  It is giving me 3/3
        // The numbers with 10 threads and 20 total assemblies are:
        // 1 run: 10 threads
        // 2 runs: 5 threads
        // 3 runs: 3 threads
        // 4 runs: 2 threads
        // 5 runs: 2 threads
        // 6 or higher runs: 1 thread

        // BUT .. the calculations do NOT take into account the number of assemblies
        // the user wants to align.  In AlignAssemblies, this is called whenever the user
        // has not specified an inParallel amount.
        val threadsToAssembliesMap = mutableMapOf<Int, Int>()
        val totalConcurrentThreads = 10
       // val totalAssemblies = 20
       // val assembliesMax = 10
        // This loop says if each assembly gets "numThreads", how many concurrent runs can we do?
        for (numThreads in 1..totalConcurrentThreads) {
            //val numRuns = assembliesMax / numThreads
            val numRuns = totalConcurrentThreads/numThreads
            val currentThreads = threadsToAssembliesMap[numRuns]
            // if currentThreads is not null and is > than numThreads, ignore.
            // otherwise, replace this entry
            if (currentThreads == null || currentThreads < numThreads) {
                threadsToAssembliesMap[numRuns] = numThreads
            }
        }
        // we should now have a map with the highest number of threads for each number of runs
        //At this point, we pick
        // 1.  if only 1 entry, use that
        // 2.  if are only 2 entries, use the one with the highest number of threads
        // 3.  if there are > 3 entries, drop the one with the lowest number of runs and the one with the highest number of runs.
        // Repeat then there are 2 entries or fewer entries left.

        // 1.  if only 1 entry, use that
        if (threadsToAssembliesMap.size == 1) {
            val entry = threadsToAssembliesMap.entries.first()
            println("Using ${entry.value} threads for ${entry.key} runs")
        } else if (threadsToAssembliesMap.size == 2) {
            // 2.  if are only 2 entries, use the one with the highest number of threads
            val entry = threadsToAssembliesMap.entries.maxByOrNull { it.value }
            println("Using ${entry!!.value} threads for ${entry.key} runs")
        } else {
            // 3.  if there are > 3 entries, drop the one with the lowest number of runs and the one with the highest number of runs.
            // Repeat then there are 2 entries or fewer entries left.
            while (threadsToAssembliesMap.size > 2) {
                val minEntry = threadsToAssembliesMap.entries.minByOrNull { it.key }
                val maxEntry = threadsToAssembliesMap.entries.maxByOrNull { it.key }
                threadsToAssembliesMap.remove(minEntry!!.key)
                threadsToAssembliesMap.remove(maxEntry!!.key)
            }
            // 2.  if are only 2 entries, use the one with the highest number of threads
            val entry = threadsToAssembliesMap.entries.maxByOrNull { it.value }
            println("Using ${entry!!.value} threads for ${entry.key} runs")
        }


    }
    @Test
    fun testCliktParams() {
        val alignAssemblies = AlignAssemblies()

        // Testing the "good" case happens in the actual test case below
        // A good test case contains a reference file, a gff file, an assemblies list file, and an output directory
        // It may optionally contain values for the --total-threads and --in-parallel parameters.

        // Test missing gff file parameter,
        val resultMissingGff =
            alignAssemblies.test(" --reference-file ${TestExtension.testRefFasta} -a ${TestExtension.smallseqAssembliesListFile} -o ${TestExtension.tempDir}")
        assertEquals(resultMissingGff.statusCode, 1)
        assertEquals(
            "Usage: align-assemblies [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --gff: --gff must not be blank\n", resultMissingGff.output
        )

        // Test missing reference file parameter,
        val resultMissingRef =
            alignAssemblies.test(" --gff ${TestExtension.smallseqAnchorsGffFile} -a ${TestExtension.smallseqAssembliesListFile} -o ${TestExtension.tempDir}")
        assertEquals(resultMissingRef.statusCode, 1)
        assertEquals(
            "Usage: align-assemblies [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --reference-file: --reference-file must not be blank\n", resultMissingRef.output
        )

        // Test missing assemblies list file parameter,
        val resultMissingAssembliesList =
            alignAssemblies.test(" --gff ${TestExtension.smallseqAnchorsGffFile} --reference-file ${TestExtension.smallseqRefFile} -o ${TestExtension.tempDir}")
        assertEquals(resultMissingAssembliesList.statusCode, 1)
        assertEquals(
            "Usage: align-assemblies [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --assemblies: --assemblies must not be blank\n", resultMissingAssembliesList.output
        )

        // Test missing output directory parameter
        val resultMissingOutputDir =
            alignAssemblies.test(" --gff ${TestExtension.smallseqAnchorsGffFile} --reference-file ${TestExtension.smallseqRefFile} -a ${TestExtension.smallseqAssembliesListFile}")
        assertEquals(resultMissingOutputDir.statusCode, 1)
        assertEquals(
            "Usage: align-assemblies [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --output-dir: --output-dir must not be blank\n", resultMissingOutputDir.output
        )
    }

    @Test
    fun testRunningAlignAssemblies() {

        // phg align-assemblies --gff /workdir/tmc46/AlignAssemblies/smallSeq_data/anchors.gff
        // --ref /workdir/tmc46/AlignAssemblies/smallSeq_data/Ref.fa
        // -a /workdir/tmc46/AlignAssemblies/assembliesList.txt
        // -o /workdir/tmc46/AlignAssemblies/temp

        val alignAssemblies = AlignAssemblies()

        val result = alignAssemblies.test(
            "--gff ${TestExtension.smallseqAnchorsGffFile} --reference-file ${TestExtension.smallseqRefFile} " +
                    "-a ${TestExtension.smallseqAssembliesListFile} -o ${TestExtension.tempDir}"
        )

        println("testRunningAlignAssemblies: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        val lineAMAF = TestExtension.tempDir + "LineA.maf"
        assertTrue(File(lineAMAF).exists(), "File $lineAMAF does not exist")

        val lineBMAF = TestExtension.tempDir + "LineB.maf"
        assertTrue(File(lineBMAF).exists(), "File $lineBMAF does not exist")

        val mafOutputA1 = TestExtension.tempDir + "LineA_unsplitA1.maf"
        val mafOutputA2 = TestExtension.tempDir + "LineA_unsplitA2.maf"
        testMergingMAF(TestExtension.smallseqLineAMafFile, mafOutputA1)
        testMergingMAF(lineAMAF, mafOutputA2)

        var checksum1 = getChecksum(mafOutputA1)
        var checksum2 = getChecksum(mafOutputA2)

        println("LineA expected checksum1: $checksum1")
        println("LineA actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "LineA.maf checksums do not match")

        val mafOutputB1 = TestExtension.tempDir + "LineA_unsplitB1.maf"
        val mafOutputB2 = TestExtension.tempDir + "LineA_unsplitB2.maf"
        testMergingMAF(TestExtension.smallseqLineBMafFile, mafOutputB1)
        testMergingMAF(lineBMAF, mafOutputB2)

        checksum1 = getChecksum(mafOutputB1)
        checksum2 = getChecksum(mafOutputB2)

        println("LineB expected checksum1: $checksum1")
        println("LineB actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "LineB.maf checksums do not match")

    }

}