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
class ExportHvcfTest {

    companion object {

        private val myLogger = LogManager.getLogger(ExportHvcfTest::class.java)

        val exportHvcfDir = "${TestExtension.tempDir}/export-vcfs/"
        val outputHvcfDir = "${exportHvcfDir}/output/"
        val inputHvcfDir = "${exportHvcfDir}/input/"
        val dbPath = "${exportHvcfDir}/tiledb_export_hvcf"
        val testHvcfFile = inputHvcfDir + "/" + File(TestExtension.smallseqRefHvcfFile).name

        @BeforeAll
        @JvmStatic
        fun setup() {

            File(exportHvcfDir).mkdirs()
            File(outputHvcfDir).mkdirs()
            File(inputHvcfDir).mkdirs()

            Initdb().createDataSets(dbPath)

            Files.copy(File(TestExtension.smallseqRefHvcfFile), File(testHvcfFile))

            bgzipAndIndexGVCFfile(testHvcfFile)

            // phg load-vcf --vcf-dir /Users/tmc46/phg_v2/ --db-path tiledb/ --temp-dir tiledb/temp/

            val result = LoadVcf().test(
                "--vcf-dir $inputHvcfDir --db-path $dbPath --temp-dir $exportHvcfDir"
            )

            println("setup: load vcf output: ${result.output}")

            assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        }

    }

    @Test
    fun testRunningExportHvcf() {

        // phg export-hvcf --dbpath tiledb --sample-names Ref -o exported-vcfs

        val result = ExportHvcf().test(
            "--dbpath $dbPath --sample-names Ref -o $outputHvcfDir"
        )

        println("testRunningExportHvcf: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        var checksum1 = getChecksum(TestExtension.smallseqRefHvcfFile)
        var checksum2 = getChecksum("$outputHvcfDir/Ref.vcf")

        println("Ref.h.vcf expected checksum1: $checksum1")
        println("Ref.vcf actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "Ref.h.vcf checksums do not match")

    }

}