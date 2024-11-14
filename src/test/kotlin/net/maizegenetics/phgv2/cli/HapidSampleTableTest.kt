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
            // These are the good hvcf files created from aligned assemblies
            val fileList = listOf("LineA","LineB","LineD","LineE")
            val hvcfFolder = "data/test/smallseq/"
            // copy the files in fileList that live in the hvcfFolder to
            // the testVCFDir
            fileList.forEach {
                val hvcfFile = File(hvcfFolder, "${it}.h.vcf")
                hvcfFile.copyTo(File(TestExtension.testVCFDir, hvcfFile.name))
            }

        }

        @JvmStatic
        @AfterAll
        fun tearDown() {
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testCliktClass() {

        // test missing hvcf dir
        val hapidSampleTable = HapidSampleTable()
        val outputFile = "${TestExtension.tempDir}/hapid_sample_table.tsv"
        val resultMissingHvcfDir =
            hapidSampleTable.test(" --output-file ${outputFile}  ")
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
        val outputFile = "${TestExtension.tempDir}/hapid_sample_table.tsv"
        val result = hapidSampleTable.test(" --hvcf-dir ${TestExtension.testVCFDir} --output-file ${outputFile} ")
        assertEquals(0, result.statusCode)

        // Check the haplotype table for the expected entries
        // Read the tab-delimited output file into a map where the first column is the key
        // and the second column is the value
        val hapidSampleMap = File(outputFile).readLines().map { line ->
            val (key, value) = line.split("\t")
            key to value
        }.toMap()

        val hapidSet = mutableSetOf<String>()
        // Read all the hvcf files in the testVCFDir
        // For each one, skip any lines that begin with #
        // For the other lines in this tab-delimited file, grab the value of the ALT field (the 5th
        // column) and add it to the set.
        File(TestExtension.testVCFDir).listFiles { file -> file.name.endsWith(".h.vcf") || file.name.endsWith(".h.vcf.gz") }
            .forEach { file ->
                file.forEachLine { line ->
                    if (!line.startsWith("#")) {
                        val hapid = line.split("\t")[4].removeSurrounding("<", ">")
                        hapidSet.add(hapid)
                    }
                }
            }

        // Check the number of entries in the hapidSet is the same as the number of haplotypes in the hvcf files
        assertEquals(hapidSet.size, hapidSampleMap.size)

        // verify select haplotype entries
        assertEquals("LineA,LineD", hapidSampleMap["13417ecbb38b9a159e3ca8c9dade7088"])
        assertEquals("LineA,LineD,LineE", hapidSampleMap["12f0cec9102e84a161866e37072443b7"])
        assertEquals("LineB", hapidSampleMap["4fc7b8af32ddd74e07cb49d147ef1938"])
        assertEquals("LineA",hapidSampleMap["d915f009b3e02030ab68baf1c43c55ad"])

        println("Finished!")
    }
}