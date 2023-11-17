package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import com.google.common.io.Files
import net.maizegenetics.phgv2.utils.bgzipAndIndexGVCFfile
import net.maizegenetics.phgv2.utils.getChecksum
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class ExportVcfTest {

    companion object {

        private val myLogger = LogManager.getLogger(ExportVcfTest::class.java)

        private val exportHvcfDir = "${TestExtension.tempDir}/export-vcfs/"
        private val outputHvcfDir = "${exportHvcfDir}/output/"
        private val inputHvcfDir = "${exportHvcfDir}/input/"
        private val dbPath = "${exportHvcfDir}/tiledb_export_hvcf"
        private val testHvcfFile = inputHvcfDir + "/" + File(TestExtension.smallseqRefHvcfFile).name

        @BeforeAll
        @JvmStatic
        fun setup() {

            File(exportHvcfDir).mkdirs()
            File(outputHvcfDir).mkdirs()
            File(inputHvcfDir).mkdirs()

            Initdb().createDataSets(dbPath)

            Files.copy(File(TestExtension.smallseqRefHvcfFile), File(testHvcfFile))

            bgzipAndIndexGVCFfile(testHvcfFile)

            // Load the vcf file into the tiledb database, so that the export has
            // something to export. This is the cli version, but using
            // programmatic way here.
            // phg load-vcf --vcf-dir /Users/tmc46/phg_v2/ --db-path tiledb/ --temp-dir tiledb/temp/

            val result = LoadVcf().test(
                "--vcf-dir $inputHvcfDir --db-path $dbPath "
            )

            println("setup: load vcf output: ${result.output}")

            assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

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
        var checksum2 = getChecksum("$outputHvcfDir/Ref.vcf")

        println("Ref.h.vcf expected checksum1: $checksum1")
        println("Ref.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "Ref.h.vcf checksums do not match")

    }

    @Test
    fun testCliktParams() {

        // Test missing db-path parameter
        val resultMissingDbpath = ExportVcf().test(
            "--sample-names Ref -o $outputHvcfDir"
        )
        assertEquals(resultMissingDbpath.statusCode, 1)
        assertEquals(
            "Usage: export-hvcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --db-path\n",
            resultMissingDbpath.output
        )

        // Test missing sample-names parameter
        val resultMissingSampleNames = ExportVcf().test(
            "--db-path $dbPath -o $outputHvcfDir"
        )
        assertEquals(resultMissingSampleNames.statusCode, 1)
        assertEquals(
            "Usage: export-hvcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --sample-names\n",
            resultMissingSampleNames.output
        )

        // Test missing output dir parameter
        val resultMissingOutputDir = ExportVcf().test(
            "--db-path $dbPath --sample-names Ref"
        )
        assertEquals(resultMissingOutputDir.statusCode, 1)
        assertEquals(
            "Usage: export-hvcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --outputDir\n",
            resultMissingOutputDir.output
        )

    }

}