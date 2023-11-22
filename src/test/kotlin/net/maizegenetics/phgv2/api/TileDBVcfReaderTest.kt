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


        //create the database
        if (!File("${dbPath}gvcf_dataset").exists()) {
            Initdb().createDataSets(dbPath)
            val loader = LoadVcf()

            val command = "--vcf-dir ${TestExtension.testVCFDir} --db-path $dbPath"
            val loaderResult = loader.test(command)
            if (loaderResult.statusCode != 0) {
                println("loader result status code = ${loaderResult.statusCode}")
                println("stderr: ${loaderResult.stderr}")
            }
        }

        //create a VCFReader for vcf_headers for gvcf data
        val path = "${dbPath}gvcf_dataset"
        println("$path is valid: ${File(path).exists()}")
        val reader = TileDBVcfReader(path)

        //check the sampleNames
        val sampleNames = reader.sampleNames().sorted()
        assertEquals("LineA",sampleNames[0])
        assertEquals("SampleLine",sampleNames[1])

        //check the headers
        //remove the ##FILTER line that is inserted by TileDB
        val headerFromTileDB = reader.headerForSample("SampleLine").lines().filter { !it.startsWith("##FILTER") }
        val headerFromFile = BufferedReader(InputStreamReader(GZIPInputStream(FileInputStream(samplegz)))).use { vcfReader ->
            vcfReader.lines().limit(10).toList()
        }

        //Test the header lines
        headerFromFile.indices.forEach { ndx ->
            assertEquals(headerFromFile[ndx], headerFromTileDB[ndx], "headers not equal at line index $ndx")
            println(headerFromTileDB[ndx])
        }

        //set a range and read some data


    }
}