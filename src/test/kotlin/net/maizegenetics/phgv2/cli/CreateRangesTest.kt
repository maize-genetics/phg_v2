package net.maizegenetics.phgv2.cli

import biokotlin.featureTree.Genome
import biokotlin.seqIO.NucSeqIO
import com.github.ajalt.clikt.testing.test
import com.google.common.collect.Range
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.createFlankingList
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.assertThrows
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
        val refBad = "data/test/smallseq/Ref.fa" // will this work with Zack's gff3 above?  It works with the smallSeq anchors.gff
        val refGood = "src/test/resources/net/maizegenetics/phgv2/cli/RefChr.fa"
        val genome = Genome.fromFile(testGffPath)
        val genes = genome.iterator().asSequence().filter { it.type == "gene" }.toList()

        val cr = CreateRanges()

        val obsIdList01 = cr.idMinMaxBounds(NucSeqIO(refGood).readAll(),genes, "gene",0)
        val obsIdList02 = cr.idMinMaxBounds(NucSeqIO(refGood).readAll(), genes, "cds",0)

        val obsIdList01Keys = obsIdList01.asMapOfRanges().keys
        var key1 = Range.closed(Position("chr1",34616),Position("chr1",40204))
        var key2 = Range.closed(Position("chr1",41213),Position("chr1",46762))

        assertEquals(2, obsIdList01.asMapOfRanges().keys.size)
        assertTrue(obsIdList01Keys.contains(key1))
        assertTrue(obsIdList01Keys.contains(key2))

        val obsIdList02Keys = obsIdList02.asMapOfRanges().keys
        key1 = Range.closed(Position("chr1",34721),Position("chr1",38366))
        key2 = Range.closed(Position("chr1",41526),Position("chr1",45913))
        assertEquals(2, obsIdList02.asMapOfRanges().keys.size)
        assertTrue(obsIdList02Keys.contains(key1))
        assertTrue(obsIdList02Keys.contains(key2))

        // This one had flanking
        // First check with a bad reference fasta - here, the chrom names don't match those in gff
        assertThrows<IllegalArgumentException> {
            //Check that an error is thrown when the bed file has overlapping intervals
            cr.idMinMaxBounds(NucSeqIO(refBad).readAll(),genes, "gene",100)
        }

        // Run again with the good reference fasta - chrom names match those in gff
        val obsIdList03 = cr.idMinMaxBounds(NucSeqIO(refGood).readAll(),genes, "gene",100)
        val obsIdList03Keys = obsIdList03.asMapOfRanges().keys
        key1 = Range.closed(Position("chr1",34516),Position("chr1",40304))
        key2 = Range.closed(Position("chr1",41113),Position("chr1",46862))
        assertEquals(2, obsIdList03.asMapOfRanges().keys.size)
        assertTrue(obsIdList03Keys.contains(key1))
        assertTrue(obsIdList03Keys.contains(key2))

        val obsBedList01 = cr.generateBedRecords(obsIdList01)
        val obsBedList02 = cr.generateBedRecords(obsIdList01)

        assertEquals(2, obsBedList01.size)

        assertEquals(BedRecord("chr1", 34616,40204,"Zm00001eb000010",0,"."), obsBedList01[0] )
        assertEquals(BedRecord("chr1",34616,40204,"Zm00001eb000010",0,"."), obsBedList02[0] )

        // Verify a bad boundary type throws an error
        assertFails {
            cr.idMinMaxBounds(NucSeqIO(refGood).readAll(), genes, "geeeene",0)
        }
    }

    @Test
    fun testCreateRangesCli() {
        //val testGffPath = "src/test/resources/net/maizegenetics/phgv2/cli/zm_b73v5_test.gff3.gz"
        val testGffPath = "data/test/smallseq/anchors.gff"
        val ref = "data/test/smallseq/Ref.fa"
        val command = CreateRanges()

        val result = command.test("--gff $testGffPath --reference-file $ref --make-only-genic")
        assertEquals(result.statusCode, 0)
        assertEquals(command.gff, testGffPath)
        assertEquals(command.boundary, "gene")
        assertEquals(command.pad, 0)
    }

    @Test
    fun testMissingGFFCLI() {
        val command = CreateRanges()

        val result = command.test("--reference-file data/test/smallseq/Ref.fa --make-only-genic")
        assertEquals(result.statusCode, 1)
        assertEquals("Usage: create-ranges [<options>]\n" +
                "\n" +
                "Error: invalid value for --gff: --gff must not be blank\n",result.output)

    }

    @Test
    fun testMissingRefCLI() {
        val command = CreateRanges()

        val result = command.test("--gff data/test/smallseq/anchors.gff --make-only-genic")
        assertEquals(result.statusCode, 1)
        assertEquals("Usage: create-ranges [<options>]\n" +
                "\n" +
                "Error: invalid value for --reference-file: --reference-file must not be blank\n",result.output)

    }

    @Test
    fun testFileOutput() {

        val testGffPath = "src/test/resources/net/maizegenetics/phgv2/cli/zm_b73v5_test.gff3.gz"
        val ref = "data/test/smallseq/Ref.fa"
        val command = CreateRanges()

        val outputFileName = "${tempDir}test.bed"

        val result = command.test("--gff $testGffPath --reference-file $ref --output $outputFileName --make-only-genic")
        assertEquals(result.statusCode, 0)
        assertEquals(command.gff, testGffPath)
        assertEquals(command.boundary, "gene")
        assertEquals(command.pad, 0)
        assertEquals(command.output, outputFileName)

        val lines = File(outputFileName).bufferedReader().readLines()
        assertEquals("chr1\t34616\t40204\tZm00001eb000010\t0\t.", lines[0])
        assertEquals("chr1\t41213\t46762\tZm00001eb000020\t0\t.", lines[1])
    }

    @Test
    fun testFillIntergenicRegions() {
        val createRanges = CreateRanges()
        //Make a list of Genic Bed records
        val records = listOf(BedRecord("chr1", 16,200,"Zm00001eb000010",0,"+"),
            BedRecord("chr1",300,400,"Zm00001eb000020",0,"+"),
            BedRecord("chr1",500,600,"Zm00001eb000030",0,"+"),
            BedRecord("chr2",200,600,"Zm00001eb000040",0,"+"))

        //Fill In the Intergenic regions
        val filledInRecords = createRanges.fillIntergenicRegions(records, NucSeqIO("data/test/createRanges/testFasta.fa").readAll())
        //Check the output
        assertEquals(10, filledInRecords.size)
        assertEquals(BedRecord("chr1", 0,16,"intergenic_chr1:0-16",0,"+"), filledInRecords[0])
        assertEquals(BedRecord("chr1", 16,200,"Zm00001eb000010",0,"+"), filledInRecords[1])
        assertEquals(BedRecord("chr1", 200,300,"intergenic_chr1:200-300",0,"+"), filledInRecords[2])
        assertEquals(BedRecord("chr1", 300,400,"Zm00001eb000020",0,"+"), filledInRecords[3])
        assertEquals(BedRecord("chr1", 400,500,"intergenic_chr1:400-500",0,"+"), filledInRecords[4])
        assertEquals(BedRecord("chr1", 500,600,"Zm00001eb000030",0,"+"), filledInRecords[5])
        assertEquals(BedRecord("chr1", 600,700,"intergenic_chr1:600-700",0,"+"), filledInRecords[6])
        assertEquals(BedRecord("chr2", 0,200,"intergenic_chr2:0-200",0,"+"), filledInRecords[7])
        assertEquals(BedRecord("chr2", 200,600,"Zm00001eb000040",0,"+"), filledInRecords[8])
        assertEquals(BedRecord("chr2", 600,700,"intergenic_chr2:600-700",0,"+"), filledInRecords[9])

    }
    @Test
    fun testConvertBedRecordsIntoOutputStrings() {
        val createRanges = CreateRanges()

        val records = listOf(BedRecord("chr1", 16,200,"Zm00001eb000010",0,"+"),
            BedRecord("chr1",300,400,"Zm00001eb000020",0,"+"),
            BedRecord("chr1",500,600,"Zm00001eb000030",0,"+"),
            BedRecord("chr2",200,600,"Zm00001eb000040",0,"+"))

        val outputStrings = createRanges.convertBedRecordsIntoOutputStrings(records)
        assertEquals(4, outputStrings.size)
        assertEquals("chr1\t16\t200\tZm00001eb000010\t0\t+", outputStrings[0])
        assertEquals("chr1\t300\t400\tZm00001eb000020\t0\t+", outputStrings[1])
        assertEquals("chr1\t500\t600\tZm00001eb000030\t0\t+", outputStrings[2])
        assertEquals("chr2\t200\t600\tZm00001eb000040\t0\t+", outputStrings[3])

    }

}
