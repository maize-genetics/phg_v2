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
    fun testSystemMemory() {
        val availMemory = AlignAssemblies().getSystemMemory()
        // Verify that memory is not 0.  It is difficult to know precisely
        // what should be returned as we do not know the system running the test.
        // This has been tested manually on multiple systems.
        assertTrue(availMemory > 0, "System memory is 0")
    }

    @Test
    fun testNumThreadsAndRuns() {
        // Test more assemblies than threads, it will pick the
        // middle option that has the highest number of parallel alignments
        // The numbers with 10 threads and 20 total assemblies are:
        // 1 run: 10 threads
        // 2 runs: 5 threads
        // 3 runs: 3 threads
        // 4 runs: 2 threads
        // 5 runs: 2 threads
        // 6 or higher runs: 1 thread

        var totalConcurrentThreads = 10
        var totalAssemblies = 20

        var alignmentsToThreads = AlignAssemblies().maximizeRunsAndThreads(totalConcurrentThreads, totalAssemblies)
        println("\nAlignAssembliesTest: alignmentsToThreads: $alignmentsToThreads")
        assertEquals(alignmentsToThreads, Pair(3, 3))

        totalAssemblies = 5
        // This picks 3/3
        // Options will be:
        // 1 run: 10 threads
        // 2 runs: 5 threads
        // 3 runs: 3 threads
        // 4 runs: 2 threads
        // 5 runs: 2 threads
        alignmentsToThreads = AlignAssemblies().maximizeRunsAndThreads(totalConcurrentThreads, totalAssemblies)
        println("\nAlignAssembliesTest: alignmentsToThreads: $alignmentsToThreads")
        assertEquals(alignmentsToThreads, Pair(3, 3))

        // test higher number of threads:
        // This picks 6/7
        // Options will be: (only highest number of runs with the same thread count is kept)
        // 1 run:  45 threads
        // 2 runs: 22 threads
        // 3 runs: 15 threads
        // 4 runs: 11 threads
        // 5 runs: 9 threads
        // 6 runs: 7 threads
        // 7 runs: 6 threads
        // 9 runs: 5 threads
        // 11 runs: 4 threads
        // 15 runs: 3 threads
        // 22 runs: 2 threads
        totalConcurrentThreads = 45
        totalAssemblies = 25
        alignmentsToThreads = AlignAssemblies().maximizeRunsAndThreads(totalConcurrentThreads, totalAssemblies)
        println("\nAlignAssembliesTest: alignmentsToThreads: $alignmentsToThreads")
        assertEquals( Pair(6,7), alignmentsToThreads)

        // test with just 1 thread and 1 run. This is to hit
        // code coverage for the case where there is only 1 option
        // for the number of runs and threads.
        totalConcurrentThreads = 1
        totalAssemblies = 1
        alignmentsToThreads = AlignAssemblies().maximizeRunsAndThreads(totalConcurrentThreads, totalAssemblies)
        println("\nAlignAssembliesTest: alignmentsToThreads: $alignmentsToThreads")
        assertEquals( Pair(1,1), alignmentsToThreads)
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
    fun testAssembliesWithTrailingSpaces() {
        // This test verifies the code handles an assembly list that
        // has blank spaces following the assembly name on the line.
        // ".trim()" was added when reading the assembly list to handle
        // this case.  This test is to verify that the code works.
        val assemblyList = TestExtension.smallseqAssembliesListFile
        val assemblyListWithSpaces = "${TestExtension.tempDir}assembliesListWithSpaces.txt"
        // create a new list by taking assembliList values and adding a space to the end of each line,
        // writing the new list to assemblyListWithSpaces
        File(assemblyList).forEachLine {
            File(assemblyListWithSpaces).appendText("$it \n")
        }

        // run assemblies, verify we have a failure
        val alignAssemblies = AlignAssemblies()

        val result = alignAssemblies.test(
            "--gff ${TestExtension.smallseqAnchorsGffFile} --reference-file ${TestExtension.smallseqRefFile} " +
                    "-a $assemblyListWithSpaces -o ${TestExtension.tempDir} --total-threads 1 --in-parallel 1"
        )

        val lineAMAF = TestExtension.tempDir + "LineA.maf"
        assertTrue(File(lineAMAF).exists(), "File $lineAMAF does not exist")

        val lineBMAF = TestExtension.tempDir + "LineB.maf"
        assertTrue(File(lineBMAF).exists(), "File $lineBMAF does not exist")

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
                    "-a ${TestExtension.smallseqAssembliesListFile} -o ${TestExtension.tempDir} --total-threads 1 --in-parallel 1"
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

    @Test
    fun testSystemDefinedThreadsAndRuns() {

        val alignAssemblies = AlignAssemblies()

        val result = alignAssemblies.test(
            "--gff ${TestExtension.smallseqAnchorsGffFile} --reference-file ${TestExtension.smallseqRefFile} " +
                    "-a ${TestExtension.smallseqAssembliesListFile} -o ${TestExtension.tempDir} "
        )

        println("testRunningAlignAssemblies: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        val lineAMAF = TestExtension.tempDir + "LineA.maf"
        assertTrue(File(lineAMAF).exists(), "File $lineAMAF does not exist")

        val lineBMAF = TestExtension.tempDir + "LineB.maf"
        assertTrue(File(lineBMAF).exists(), "File $lineBMAF does not exist")

        // Checksums of files are not verified here.  Running in parallel means
        // the entries in the maf files are in unspecified order.  The files, while containing
        // the same data, will have different checksums.  This is not a problem, as the
        // order of the entries in the maf files is not important.


    }

    @Test
    fun testRequestedThreadsGreaterThanAvail() {
        // I'm assuming all of the machines which will run this test will
        // have less than 300 threads available for processing.  The code
        // should default to the maximum number of threads available (minus 2 for IO)
        // This test ups our code coverage.
        val alignAssemblies = AlignAssemblies()

        val result = alignAssemblies.test(
            "--gff ${TestExtension.smallseqAnchorsGffFile} --reference-file ${TestExtension.smallseqRefFile} " +
                    "-a ${TestExtension.smallseqAssembliesListFile} -o ${TestExtension.tempDir} --total-threads 300 --in-parallel 1"
        )

        println("testRunningAlignAssemblies: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        val lineAMAF = TestExtension.tempDir + "LineA.maf"
        assertTrue(File(lineAMAF).exists(), "File $lineAMAF does not exist")

        val lineBMAF = TestExtension.tempDir + "LineB.maf"
        assertTrue(File(lineBMAF).exists(), "File $lineBMAF does not exist")

    }

    @Test
    fun testInParallelGreaterThanNumberOfAssemblies() {
        // There are only 2 assemblies that we are aligning, will give it
        // in-parallel of 4.  It should adjust to 2.
        // This test ups our code coverage.
        val alignAssemblies = AlignAssemblies()

        val result = alignAssemblies.test(
            "--gff ${TestExtension.smallseqAnchorsGffFile} --reference-file ${TestExtension.smallseqRefFile} " +
                    "-a ${TestExtension.smallseqAssembliesListFile} -o ${TestExtension.tempDir} --total-threads 4 --in-parallel 4"
        )

        println("testRunningAlignAssemblies: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        val lineAMAF = TestExtension.tempDir + "LineA.maf"
        assertTrue(File(lineAMAF).exists(), "File $lineAMAF does not exist")

        val lineBMAF = TestExtension.tempDir + "LineB.maf"
        assertTrue(File(lineBMAF).exists(), "File $lineBMAF does not exist")

    }

}