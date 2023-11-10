package net.maizegenetics.phgv2.cli

import biokotlin.featureTree.Genome
import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import java.io.File
import kotlin.test.assertEquals
import kotlin.test.Test
import kotlin.test.assertFails

class CreateRangesTest {


    companion object {
        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"

        @JvmStatic
        @BeforeAll
        fun setup() {
            File(tempDir).mkdirs()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(tempDir).deleteRecursively()
        }
    }
    @Test
    fun evaluateMethods() {
        assertEquals(2, 2)
        val testGffPath = "src/test/resources/net/maizegenetics/phgv2/cli/zm_b73v5_test.gff3.gz"
        val genome = Genome.fromFile(testGffPath)
        val genes = genome.iterator().asSequence().filter { it.type == "gene" }.toList()


        val cr = CreateRanges()

        val obsIdList01 = cr.idMinMaxBounds(genes, "gene", 0)
        val obsIdList02 = cr.idMinMaxBounds(genes, "cds", 0)
        val obsIdList03 = cr.idMinMaxBounds(genes, "gene", 100)

        val obsBedList01 = cr.generateBedRows(obsIdList01, genes)
        val obsBedList02 = cr.generateBedRows(obsIdList01, genes, ",")
        val obsBedList03 = cr.generateBedRows(obsIdList02, genes, featureId = "biotype")

        assertEquals(2, obsIdList01.size)
        assertEquals(Pair(34616, 40204), obsIdList01[0]) // should be 34616
        assertEquals(Pair(41213, 46762), obsIdList01[1])

        assertEquals(2, obsIdList02.size)
        assertEquals(Pair(34721, 38366), obsIdList02[0])
        assertEquals(Pair(41526, 45913), obsIdList02[1])

        assertEquals(2, obsIdList03.size)
        assertEquals(Pair(34516, 40304), obsIdList03[0])
        assertEquals(Pair(41113, 46862), obsIdList03[1])

        assertEquals(2, obsBedList01.size)
        assertEquals("chr1\t34616\t40204\tZm00001eb000010\t0\t+", obsBedList01[0] )
        assertEquals("chr1,34616,40204,Zm00001eb000010,0,+", obsBedList02[0] )
        assertEquals("chr1\t34721\t38366\tprotein_coding\t0\t+", obsBedList03[0])

        assertFails {
            cr.idMinMaxBounds(genes, "geeeene", 0)
        }
    }

    @Test
    fun testCreateRangesCli() {
        val testGffPath = "src/test/resources/net/maizegenetics/phgv2/cli/zm_b73v5_test.gff3.gz"
        val command = CreateRanges()

        val result = command.test("--gff $testGffPath --make-only-genic")
        assertEquals(result.statusCode, 0)
        assertEquals(command.gff, testGffPath)
        assertEquals(command.boundary, "gene")
        assertEquals(command.pad, 0)
    }

    @Test
    fun testMissingGFFCLI() {
        val command = CreateRanges()

        val result = command.test("")
        assertEquals(result.statusCode, 1)
        assertEquals("Usage: create-ranges [<options>]\n" +
                "\n" +
                "Error: invalid value for --gff: --gff must not be blank\n",result.output)

    }

    @Test
    fun testFileOutput() {

        val testGffPath = "src/test/resources/net/maizegenetics/phgv2/cli/zm_b73v5_test.gff3.gz"
        val command = CreateRanges()

        val outputFileName = "${tempDir}test.bed"

        val result = command.test("--gff $testGffPath --output $outputFileName --make-only-genic")
        assertEquals(result.statusCode, 0)
        assertEquals(command.gff, testGffPath)
        assertEquals(command.boundary, "gene")
        assertEquals(command.pad, 0)
        assertEquals(command.output, outputFileName)

        val lines = File(outputFileName).bufferedReader().readLines()
        assertEquals("chr1\t34616\t40204\tZm00001eb000010\t0\t+", lines[0])
        assertEquals("chr1\t41213\t46762\tZm00001eb000020\t0\t-", lines[1])
    }
}
