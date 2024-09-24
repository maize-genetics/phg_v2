package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import com.google.common.io.Files
import net.maizegenetics.phgv2.brapi.createSmallSeqTiledb
import net.maizegenetics.phgv2.brapi.resetDirs
import net.maizegenetics.phgv2.utils.bgzipAndIndexGVCFfile
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class ListSamplesTest {
    companion object {

        @JvmStatic
        @BeforeAll
        fun setup() {
            // delete, reset the directories
            resetDirs()

            // create the tiledb datasets, load them with from the vcf files
            // This will also create the AGC compressed file
            // tiledbURI:String = TestExtension.testTileDBURI
            createSmallSeqTiledb()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testRunningListSamples() {
        val outputFile = "${TestExtension.tempDir}list_samples.txt"
        var result = ListSamples().test(
            "--db-path ${TestExtension.testTileDBURI} --output-file $outputFile"
        )
        assertEquals(result.statusCode, 0)

        // Read the output file.
        // Verify the file has a header line with columns SampleName, InAGC, InGVCF, InHVCF
        // Verify the file Has a "Y" for all samples in the InAGC and InHVCF columns, and "N" for all samples in the InGVCF column.
        val outFileLines = File(outputFile).readLines()
        assertEquals(outFileLines.size, 4)
        assertEquals(outFileLines[0], "SampleName\tInAGC\tInGVCF\tInHVCF")
        assertEquals(outFileLines[1], "LineA\tY\tN\tY")
        assertEquals(outFileLines[2], "LineB\tY\tN\tY")
        assertEquals(outFileLines[3], "Ref\tY\tN\tY")

        // On setup, we called createSmallSeqTiledb() which calls loadVcfFiles(tiledbURI)
        // This is from the brapiSetup code, and it only loads LineA.h.vcf and LineB.h.vcf
        // files.  No gvcf files were loaded, so load one now.

        // Copy the LineA gvcf file.
        val origGvcfFile = "data/test/smallseq/LineA.g.vcf"
        val testGvcfFile = "${TestExtension.testVCFDir}LineA.g.vcf"
        Files.copy(File(origGvcfFile), File(testGvcfFile))

        // call bgzipAndIndexGVCFfile to zip and index the file
        bgzipAndIndexGVCFfile(testGvcfFile)

        // Load the LineA gvcf file to tileDB
        val loadVCF = LoadVcf()
        val vcfDir = TestExtension.testVCFDir
        result = loadVCF.test("--vcf-dir ${vcfDir} --db-path ${TestExtension.testTileDBURI} ")
        assertEquals(result.statusCode, 0)

        val outputFile2 = "${TestExtension.tempDir}list_samples2.txt"
        // Re-run ListSamples command.  LineA should now have Y in the InGVCF column.
        result = ListSamples().test(
            "--db-path ${TestExtension.testTileDBURI} --output-file $outputFile2"
        )
        assertEquals(result.statusCode, 0)

        // Verify the new file
        val outFileLines2 = File(outputFile2).readLines()
        assertEquals(outFileLines2.size, 4)
        assertEquals(outFileLines2[0], "SampleName\tInAGC\tInGVCF\tInHVCF")
        assertEquals(outFileLines2[1], "LineA\tY\tY\tY")
        assertEquals(outFileLines2[2], "LineB\tY\tN\tY")
        assertEquals(outFileLines2[3], "Ref\tY\tN\tY")

        println("testRunningListSamples: result output: ${result.output}")
    }
}