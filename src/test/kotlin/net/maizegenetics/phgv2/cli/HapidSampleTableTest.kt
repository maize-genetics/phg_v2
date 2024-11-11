package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.assertEquals

class HapidSampleTableTest {

    companion object {

        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.testVCFDir).mkdirs()
            val hvcfFile = File("data/test/smallseq/LineImpute.h.vcf")
            // Copy hvcf file to testVCFDir
            hvcfFile.copyTo(File(TestExtension.testVCFDir, hvcfFile.name))

        }

        @JvmStatic
        @AfterAll
        fun tearDown() {
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testCliktClass() {

        // This file has a samplename of "LineImpute".
        // All of the haplotypes for chrom1 are from LineA, all of the haplotypes for chrom2 are from LineB
        val hapidSampleTable = HapidSampleTable()
        val outputDir = TestExtension.tempDir
        val resultMissingHvcfDir =
            hapidSampleTable.test(" --output-dir ${outputDir}  ")
        assertEquals(1, resultMissingHvcfDir.statusCode )
        assertEquals(
            "Usage: hapid-sample-table [<options>]\n" +
                    "\n" +
                    "Error: missing option --hvcf-dir\n", resultMissingHvcfDir.output
        )
    }

    @Test
    fun testHapidSampleTable() {
        val hapidSampleTable = HapidSampleTable()
        val outputDir = TestExtension.tempDir
        val result = hapidSampleTable.test(" --hvcf-dir ${TestExtension.testVCFDir} --output-dir ${outputDir} ")
        assertEquals(0, result.statusCode)

        // THis should result in either LineA or LineB as the sample mapping to the hapid, but in all
        // cases it is coming back as LineImpute

        // Check entries in the table
        println("FInished!")
    }
}