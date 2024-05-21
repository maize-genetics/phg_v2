package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.extension.ExtendWith
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.assertEquals

/**
 * This class created to test the PrepareSlurmAlignFile class
 */

@ExtendWith(TestExtension::class)
class PrepareSlurmAlignFileTest {

    companion object {

        //val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"
        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.testTileDBURI).mkdirs()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testCliktCommands() {
        val prepareAssemblies = PrepareSlurmAlignFile()
        val refFasta = TestExtension.smallseqRefFile
        val refSam = TestExtension.smallseqRefFile
        val refCDSfasta = TestExtension.smallseqRefFile
        val refCDSSam = TestExtension.smallseqRefFile
        val assembliesList = TestExtension.smallseqAssembliesListFile

        // Test missing gff parameter
        // LCJ - this needs ALL the parameters except gff or you get error message on all that are missing.
        val resultMissingGFF = prepareAssemblies.test("--reference-file ${refFasta} --reference-cds-sam ${refSam} --reference-cds-fasta ${refCDSfasta} --assemblies ${assembliesList}")
        assertEquals(1, resultMissingGFF.statusCode )
        assertEquals("Usage: prepare-slurm-align-file [<options>]\n" +
                "\n" +
                "Error: invalid value for --gff: --gff must not be blank\n", resultMissingGFF.output)

//        // Test missing reference-file parameter
//        val resultMissingRefFile = prepareAssemblies.test("--gff ${TestExtension.smallseqAnchorsGffFile} --reference-cds-sam ${refSam} --reference-cds-fasta ${refCDSfasta} --assemblies ${assembliesList}")
//        assertEquals(1, resultMissingRefFile.statusCode )
//        assertEquals("Usage: prepare-slurm-align-file [<options>]\n" +
//                "\n" +
//                "Error: invalid value for --reference-file: --reference-file must not be blank\n", resultMissingRefFile.output)
//
//        // Test missing reference-cds-sam parameter
//        val resultMissingRefSam = prepareAssemblies.test("--gff ${TestExtension.smallseqAnchorsGffFile} --reference-file ${refFasta} --reference-cds-fasta ${refCDSfasta} --assemblies ${assembliesList}")
//        assertEquals(1, resultMissingRefSam.statusCode )
//        assertEquals("Usage: prepare-slurm-align-file [<options>]\n" +
//                "\n" +
//                "Error: invalid value for --reference-cds-sam: --reference-cds-sam must not be blank\n", resultMissingRefSam.output)

        // Test missing reference-cds-fasta parameter
        //

    }
    @Test
    fun testPrepareSlurmAlignFile() {

        // FIrst, we need to run alignAssemblies with the just-ref-prep option to create the reference cds sam file
        val alignAssemblies = AlignAssemblies()
        val outputDir = TestExtension.tempDir
        val referenceFile = TestExtension.smallseqRefFile
        val prepResult = alignAssemblies.test(
            "--gff ${TestExtension.smallseqAnchorsGffFile} --reference-file ${TestExtension.smallseqRefFile} " +
                    "--assembly-file-list ${TestExtension.smallseqAssembliesListFile} -o ${outputDir} --total-threads 1 --in-parallel 1 --just-ref-prep"
        )
        // verify good result
        assertEquals(0, prepResult.statusCode, "status code not 0 for AlignAssemblies with just-ref-prep: ${prepResult.statusCode}")
        // verify no MAF files in the output directory
        val mafFiles = File(TestExtension.tempDir).listFiles { _, name -> name.endsWith(".maf") }
        assertEquals(mafFiles.size, 0, "MAF files found in output directory")

        // Verify the reference cds sam and refCds fasta files were created
        val justNameRef = File(referenceFile).nameWithoutExtension
        val samOutFile = "${justNameRef}.sam"
        val refCDSfasta = "${outputDir}/ref.cds.fasta"
        val refCDSsam = "${outputDir}/${samOutFile}"

        assertEquals(File(refCDSfasta).exists(), true, "Reference CDS fasta file not created")
        assertEquals(File(refCDSsam).exists(), true, "Reference CDS sam file not created")

        // Using the files created above, create the slurm align file
        val prepareAssemblies = PrepareSlurmAlignFile()

        println("Finished just-ref-prep, now running PrepareSlurmAlignFile")
        val result = prepareAssemblies.test(
            "--gff ${TestExtension.smallseqAnchorsGffFile} --reference-file ${TestExtension.smallseqRefFile} " +
                    "--reference-cds-sam ${refCDSsam} --reference-cds-fasta ${refCDSfasta} " +
                    "--slurm-command-file ${TestExtension.tempDir}/slurm_align.sh " +
                    "--assemblies ${TestExtension.smallseqAssembliesListFile} -o ${TestExtension.tempDir} --total-threads 1 --in-parallel 1"
        )
        // verify good result
        assertEquals(0, result.statusCode,  "status code not 0 for PrepareSlurmAlignFile: ${result.statusCode}")

        // verify the slurm align file was created
        val slurmFile = File("${TestExtension.tempDir}/slurm_align.sh")
        assertEquals(slurmFile.exists(), true, "Slurm align file not created")

        // Veriy the slurm file contains 2 lines
        val lines = slurmFile.readLines()
        assertEquals(lines.size, 2, "Slurm file does not contain 2 lines")
        // Verify both lines start with "phg "
        lines.forEach {
            assertEquals(it.startsWith("phg "), true, "Line does not start with 'phg '")
        }
    }

    @Test
    fun testPairOfBlank() {
        val refCDSfasta = ""
        val refCDSsam = ""
        var anchorwaveRefFiles = Pair(refCDSfasta,refCDSsam)
        println("Anchorwave reference files: $anchorwaveRefFiles")
    }

}