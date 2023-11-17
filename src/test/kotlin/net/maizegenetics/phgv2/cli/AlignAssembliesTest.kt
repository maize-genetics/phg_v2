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