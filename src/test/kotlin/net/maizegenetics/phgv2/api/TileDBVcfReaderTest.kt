package net.maizegenetics.phgv2.api

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.Initdb
import net.maizegenetics.phgv2.cli.LoadVcf
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.bgzipAndIndexGVCFfile
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.BufferedReader
import java.io.File
import java.io.FileInputStream
import java.io.InputStreamReader
import java.util.zip.GZIPInputStream
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class TileDBVcfReaderTest {

    @Test
    fun testGvcfReading() {
        //create a gvcf tiledb-vcf database

        //copy the gvcf files to a temp vcf directory
        val sampleGvcf = File(TestExtension.testVCFDir).resolve("sample.g.vcf")
        val samplegz = File(TestExtension.testVCFDir).resolve("sample.g.vcf.gz")
        val lineAGvcf = File(TestExtension.testVCFDir).resolve(File(TestExtension.smallSeqLineAGvcfFile).name)
        File(TestExtension.smallSeqLineAGvcfFile).copyTo(lineAGvcf, true)
        File(TestExtension.smallseqGvcfFile).copyTo(sampleGvcf, true)

        println("bgzipping and indexing files")
        bgzipAndIndexGVCFfile(sampleGvcf.absolutePath)
        bgzipAndIndexGVCFfile(lineAGvcf.absolutePath)

        val dbPath = TestExtension.testTileDBURI

        //delete the database, if it exists
        val dbFile = File(dbPath)
        if (dbFile.exists()) {
            dbFile.deleteRecursively()
        }

        //create the database
        Initdb().createDataSets(dbPath)
        val loader = LoadVcf()

        val command = "--vcf-dir ${TestExtension.testVCFDir} --db-path $dbPath"
        val loaderResult = loader.test(command)
        if (loaderResult.statusCode != 0) {
            println("loader result status code = ${loaderResult.statusCode}")
            println("stderr: ${loaderResult.stderr}")
        }

        //create a VCFReader for vcf_headers
        val path = "${dbPath}gvcf_dataset"
        println("$path is valid: ${File(path).exists()}")
        val reader = TileDBVcfReader(path)

        //check the sampleNames
        val sampleNames = reader.sampleNames().sorted()
        assertEquals("LineA",sampleNames[0])
        assertEquals("SampleLine",sampleNames[1])

        //check the headers
        //remove the ##FILTER line that is inserted by TileDB
        //set the initial buffer size to a small value for testing
        val headerFromTileDB = reader.headerForSample("SampleLine", 500).lines().filter { !it.startsWith("##FILTER") }
        val headerFromFile = BufferedReader(InputStreamReader(GZIPInputStream(FileInputStream(samplegz)))).use { vcfReader ->
            vcfReader.lines().limit(10).toList()
        }

        //Test the header lines
        headerFromFile.indices.forEach { ndx ->
            assertEquals(headerFromFile[ndx], headerFromTileDB[ndx], "headers not equal at line index $ndx")
            println(headerFromTileDB[ndx])
        }

        //set a range and read some data
        reader.ranges(listOf("1:1-1000"))
        val dataResult = reader.data()

        //expected sample names
        val expectedNames = listOf("LineA", "SampleLine")
        val names = dataResult.map { it.sampleName }.distinct().sorted()
        (0..1).forEach { assertEquals(expectedNames[it], names[it], "sample names do not match") }

        //expected number of records per sample
        //LineA 60, SampleLine 55
        val lineAExpectedCount = 60
        val sampleExpectedCount = 55
        val lineACount = dataResult.filter { it.sampleName == "LineA" }.count()
        val sampleCount = dataResult.filter { it.sampleName == "SampleLine" }.count()

        println("name, contig, start, end, genotype, AD, DP")
        dataResult.filter { it.sampleName == "LineA" }.forEach { println("${it.sampleName}, ${it.contig}, ${it.startPos}, ${it.endPos}, ${it.genotype}, ${it.AD}, ${it.DP}") }
        dataResult.filter { it.sampleName == "SampleLine" }.forEach { println("${it.sampleName}, ${it.contig}, ${it.startPos}, ${it.endPos}, ${it.genotype}, ${it.AD}, ${it.DP}") }
        assertEquals(lineAExpectedCount, lineACount, "counts of data for lineA")
        assertEquals(sampleExpectedCount, sampleCount, "counts of data for SampleLine")

        //create a new reader to test setting the sample name
        val oneSampleReader = TileDBVcfReader(path, listOf("SampleLine"))
        oneSampleReader.ranges(listOf("1:1-1000"))
        val oneSampleResult = oneSampleReader.data()
        val lineARecount = oneSampleResult.filter { it.sampleName == "LineA" }.count()
        val sampleRecount = oneSampleResult.filter { it.sampleName == "SampleLine" }.count()
        assertEquals(0, lineARecount, "recounts of data for lineA")
        assertEquals(sampleExpectedCount, sampleRecount, "recounts of data for SampleLine")

    }
}