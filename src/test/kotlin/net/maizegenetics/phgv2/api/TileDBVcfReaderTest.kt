package net.maizegenetics.phgv2.api

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.LoadVcf
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import java.nio.file.Files
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class TileDBVcfReaderTest {

    @Test
    fun testGvcfReading() {
        //create a gvcf tiledb-vcf database

        //copy the gvcf files to a temp vcf directory
        File(TestExtension.smallSeqLineAGvcfFile).copyTo(File(TestExtension.testVCFDir))
        File(TestExtension.smallseqGvcfFile).copyTo(File(TestExtension.testVCFDir))
        val sampleGvcf = File(TestExtension.testVCFDir).resolve(File(TestExtension.smallseqGvcfFile).name)

        //create the database
        val loader = LoadVcf()

        loader.test("--vcf-dir ${TestExtension.testVCFDir} --db-path ${TestExtension.testTileDBURI}")

        //create a VCFReader
        val reader = TileDBVcfReader(TestExtension.testTileDBURI)

        //check the sampleNames
        val sampleNames = reader.sampleNames().sorted()
        assertEquals("LineA",sampleNames[0])
        assertEquals("SampleLine",sampleNames[1])

        //check the headers
        val headerFromTileDB = reader.headerForSample("SampleLine").lines()
        val headerFromFile = Files.newBufferedReader(sampleGvcf.toPath()).use { vcfReader ->
            vcfReader.lines().limit(10).toList()
        }

        headerFromFile.indices.forEach { ndx ->
            assertEquals(headerFromFile[ndx], headerFromTileDB[ndx], "headers not equal at line index $ndx")
        }

        //set a range and read some data


    }
}