package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.brapi.createSmallSeqTiledb
import net.maizegenetics.phgv2.brapi.resetDirs
import net.maizegenetics.phgv2.utils.getChecksum
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
class ExportVcfTest {

    companion object {

        private val myLogger = LogManager.getLogger(ExportVcfTest::class.java)

        private val exportHvcfDir = "${TestExtension.tempDir}export-vcfs/"
        private val outputHvcfDir = "${exportHvcfDir}output/"
        private val inputHvcfDir = "${exportHvcfDir}input/"
        private val dbPath = "${exportHvcfDir}/tiledb_export_hvcf"
        private val testHvcfFile = inputHvcfDir + "/" + File(TestExtension.smallseqRefHvcfFile).name

        @BeforeAll
        @JvmStatic
        fun setup() {
            // delete, reset the directories
            resetDirs()

            // create the tiledb datasets, load them from the vcf files
            // This will also create the AGC compressed file
            createSmallSeqTiledb(dbPath)

            File(exportHvcfDir).mkdirs()
            File(outputHvcfDir).mkdirs()
            File(inputHvcfDir).mkdirs()

        }

        @AfterAll
        @JvmStatic
        fun teardown() {
            File(TestExtension.tempDir).deleteRecursively()
        }

    }

    @Test
    fun testRunningExportHvcf() {

        // phg export-hvcf --db-path tiledb --sample-names Ref -o exported-vcfs

        val result = ExportVcf().test(
            "--db-path $dbPath --sample-names Ref -o $outputHvcfDir"
        )

        println("testRunningExportHvcf: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        var checksum1 = getChecksum(TestExtension.smallseqRefHvcfFile)
        var checksum2 = getChecksum("$outputHvcfDir/Ref.h.vcf")

        println("Ref.h.vcf expected checksum1: $checksum1")
        println("Ref.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "Ref.h.vcf checksums do not match")

        //test using a regions-file without .vcf or .bed extension
        assertThrows<IllegalArgumentException>() {
            val regionResult = ExportVcf().test(
            "--db-path $dbPath --sample-names Ref -o $outputHvcfDir --regions-file ${TestExtension.testKeyFile}"
        )}

    }

    @Test
    fun testMultipleSamplesFromList() {
        val result = ExportVcf().test(
            "--db-path $dbPath --sample-names Ref,LineA -o $outputHvcfDir"
        )

        println("testRunningExportHvcf: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        // Verify Ref.vcf checksum
        var checksum1 = getChecksum(TestExtension.smallseqRefHvcfFile)
        var checksum2 = getChecksum("$outputHvcfDir/Ref.h.vcf")

        println("Ref.h.vcf expected checksum1: $checksum1")
        println("Ref.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "Ref.h.vcf checksums do not match")

        // Get checksum for LineA
        checksum1 = getChecksum(TestExtension.smallseqLineAHvcfFile)
        checksum2 = getChecksum("$outputHvcfDir/LineA.h.vcf")
        // verify checksums match
        assertEquals(checksum1, checksum2, "LineA.h.vcf checksums do not match")
    }

    @Test
    fun testMultipleSamplesFromFile() {
        // write a test file that has 2 sample names in it: Ref and LineA, each on a separate line
        val sampleFile = File("$exportHvcfDir/sample-names.txt")
        sampleFile.writeText("Ref\nLineA")


        val result = ExportVcf().test(
            "--db-path $dbPath --sample-file $sampleFile -o $outputHvcfDir"
        )

        println("testRunningExportHvcf: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        // Verify ref vcf
        var checksum1 = getChecksum(TestExtension.smallseqRefHvcfFile)
        var checksum2 = getChecksum("$outputHvcfDir/Ref.h.vcf")

        println("Ref.h.vcf expected checksum1: $checksum1")
        println("Ref.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "Ref.h.vcf checksums do not match")

        // Get checksum for LineA
        checksum1 = getChecksum(TestExtension.smallseqLineAHvcfFile)
        checksum2 = getChecksum("$outputHvcfDir/LineA.h.vcf")
        // verify checksums match
        assertEquals(checksum1, checksum2, "LineA.h.vcf checksums do not match")
    }

    @Test
    fun testCliktParams() {

        // Test missing sample-names parameter
        val resultMissingSampleNames = ExportVcf().test(
            "--db-path $dbPath -o $outputHvcfDir"
        )
        assertEquals(resultMissingSampleNames.statusCode, 1)
        assertEquals(
            "Usage: export-vcf [<options>]\n" +
                    "\n" +
                    "Error: must provide one of --sample-names, --sample-file\n",
            resultMissingSampleNames.output
        )

        // Test missing output dir parameter
        val resultMissingOutputDir = ExportVcf().test(
            "--db-path $dbPath --sample-names Ref"
        )
        assertEquals(resultMissingOutputDir.statusCode, 1)
        assertEquals(
            "Usage: export-vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --outputDir\n",
            resultMissingOutputDir.output
        )

    }

}