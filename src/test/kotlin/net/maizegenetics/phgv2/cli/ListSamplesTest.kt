package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import com.google.common.io.Files
import net.maizegenetics.phgv2.brapi.createSmallSeqTiledb
import net.maizegenetics.phgv2.brapi.resetDirs
import net.maizegenetics.phgv2.utils.bgzipAndIndexGVCFfile
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertThrows
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
    fun testRunningListSamplesAll() {
        // all and gvcf are tested together because the "beforeAll" sets up the database
        // with only hvcf files.  This db persists through all the tests.  Initially there
        // are no gvcf files in the db.  We test, then add a gvcf filea and test again.
        // If the entire suite of tests in this file are run at once, the db for subsequent
        // tests has the gvcf file in it, making testing with/without gvcf samples more time intensive
        // as it would be necessary to delete from the db.  But that would fail if the test
        // was then run individually.  So, the tests are run in two groups, all/gvcf and hvcf.
        val outputFile = "${TestExtension.tempDir}list_samples.txt"
        var result = ListSamples().test(
            "--db-path ${TestExtension.testTileDBURI} --output-file $outputFile --data-set all"
        )
        assertEquals(0,result.statusCode)

        // Read the output file.
        // Verify the file has a header line with columns SampleName, InAGC, InGVCF, InHVCF
        // Verify the file Has a "Y" for all samples in the InAGC and InHVCF columns, and "N" for all samples in the InGVCF column.
        var outFileLines = File(outputFile).readLines()
        assertEquals(outFileLines.size, 4)
        assertEquals(outFileLines[0], "SampleName\tInAGC\tInGVCF\tInHVCF")
        assertEquals(outFileLines[1], "LineA\tY\tN\tY")
        assertEquals(outFileLines[2], "LineB\tY\tN\tY")
        assertEquals(outFileLines[3], "Ref\tY\tN\tY")

        // Test the gvcf option
        val outputFileGvcf = "${TestExtension.tempDir}sampleGvcf.txt"
        result = ListSamples().test(
            "--db-path ${TestExtension.testTileDBURI} --output-file $outputFileGvcf --data-set gvcf"
        )
        outFileLines = File(outputFileGvcf).readLines()
        assertEquals(0,result.statusCode)
        assertEquals(0, outFileLines.size) // no gvcf samples yet

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
        // Re-run ListSamples command for "all".  LineA should now have Y in the InGVCF column.
        result = ListSamples().test(
            "--db-path ${TestExtension.testTileDBURI} --output-file $outputFile2 --data-set all"
        )
        assertEquals(result.statusCode, 0)

        // Verify the new file
        val outFileLines2 = File(outputFile2).readLines()
        assertEquals(outFileLines2.size, 4)
        assertEquals(outFileLines2[0], "SampleName\tInAGC\tInGVCF\tInHVCF")
        assertEquals(outFileLines2[1], "LineA\tY\tY\tY")
        assertEquals(outFileLines2[2], "LineB\tY\tN\tY")
        assertEquals(outFileLines2[3], "Ref\tY\tN\tY")

        // re-run the gvcf option
        val outputFileGvcf2 = "${TestExtension.tempDir}sampleGvcf2.txt"
        result = ListSamples().test(
            "--db-path ${TestExtension.testTileDBURI} --output-file $outputFileGvcf2 --data-set gvcf"
        )
        outFileLines = File(outputFileGvcf2).readLines()
        assertEquals(0,result.statusCode)
        assertEquals(1, outFileLines.size) // LineA should now be in the gvcf file
        assertEquals(outFileLines[0], "LineA")

    }

    @Test
    fun testListSamplesForHvcf() {
        val outputFile = "${TestExtension.tempDir}listHvcfSamples.txt"
        var result = ListSamples().test(
            "--db-path ${TestExtension.testTileDBURI} --output-file $outputFile --data-set hvcf"
        )
        assertEquals(0,result.statusCode)

        val outFileLines = File(outputFile).readLines()
        assertEquals(outFileLines.size, 3)
        assertEquals(outFileLines[0], "LineA")
        assertEquals(outFileLines[1], "LineB")
        assertEquals(outFileLines[2], "Ref")
    }

    @Test
    fun testListSamplesDefault() {
        val outputFile = "${TestExtension.tempDir}listHvcfSamples.txt"
        var result = ListSamples().test(
            "--db-path ${TestExtension.testTileDBURI} --output-file $outputFile "
        )
        assertEquals(0,result.statusCode)

        // Should run and print the samples in the hvcf file (which is the default)
        val outFileLines = File(outputFile).readLines()
        assertEquals(outFileLines.size, 3)
        assertEquals(outFileLines[0], "LineA")
        assertEquals(outFileLines[1], "LineB")
        assertEquals(outFileLines[2], "Ref")
    }

    @Test
    fun testBadDatasetOption() {
        val outputFile = "${TestExtension.tempDir}listHvcfSamples.txt"
        var result = ListSamples().test(
            "--db-path ${TestExtension.testTileDBURI} --output-file $outputFile --data-set bad"
        )
        assertEquals(1,result.statusCode)

        // Cannot get this error message correct !!!!
//        assertEquals("Usage: list-samples [<options>]\n" +
//                "\n" ,result.output)
//        assertEquals("Usage: list-samples [<options>]\n" +
//                "\n" +
//                "Error: invalid value for --data-set: Invalid dataset option. Available options are: agc, gvcf, hvcf, all\n",result.output)


    }
}