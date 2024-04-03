package net.maizegenetics.phgv2.cli

import biokotlin.featureTree.Genome
import biokotlin.seqIO.NucSeqIO
import com.github.ajalt.clikt.testing.test
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import com.google.common.collect.TreeRangeMap
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.addRange
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

        val obsIdList01 = cr.idMinMaxBounds(NucSeqIO(refGood).readAll(),genes, "gene")
        val obsIdList02 = cr.idMinMaxBounds(NucSeqIO(refGood).readAll(), genes, "cds")

        val obsIdList01Keys = obsIdList01.asMapOfRanges().keys
        var key1 = Range.closed(Position("chr1",34617),Position("chr1",40204))
        var key2 = Range.closed(Position("chr1",41214),Position("chr1",46762))

        assertEquals(2, obsIdList01.asMapOfRanges().keys.size)
        assertTrue(obsIdList01Keys.contains(key1))
        assertTrue(obsIdList01Keys.contains(key2))

        val obsIdList02Keys = obsIdList02.asMapOfRanges().keys
        key1 = Range.closed(Position("chr1",34722),Position("chr1",38366))
        key2 = Range.closed(Position("chr1",41527),Position("chr1",45913))
        assertEquals(2, obsIdList02.asMapOfRanges().keys.size)
        assertTrue(obsIdList02Keys.contains(key1))
        assertTrue(obsIdList02Keys.contains(key2))



        // Run again with the good reference fasta - chrom names match those in gff
        val obsIdList03 = cr.idMinMaxBounds(NucSeqIO(refGood).readAll(),genes, "gene")
        val obsIdList03Keys = obsIdList03.asMapOfRanges().keys
        key1 = Range.closed(Position("chr1",34617),Position("chr1",40204))
        key2 = Range.closed(Position("chr1",41214),Position("chr1",46762))
        assertEquals(2, obsIdList03.asMapOfRanges().keys.size)
        assertTrue(obsIdList03Keys.contains(key1))
        assertTrue(obsIdList03Keys.contains(key2))

        // Add flanking
        // First test with a bad reference fasta - here, the chrom names don't match those in gff
        assertThrows<IllegalArgumentException> {
            //Check that an error is thrown when the bed file has overlapping intervals
            createFlankingList(obsIdList03, 100,NucSeqIO(refBad).readAll())
        }

        // Run again with the good reference fasta - chrom names match those in gff
        val obsIdList03Flanking = createFlankingList(obsIdList03, 100,NucSeqIO(refGood).readAll())
        val obsIdList03FlankingKeys = obsIdList03Flanking.asMapOfRanges().keys
        key1 = Range.closed(Position("chr1",34517),Position("chr1",40304))
        key2 = Range.closed(Position("chr1",41114),Position("chr1",46862))
        assertEquals(2, obsIdList03.asMapOfRanges().keys.size)
        assertTrue(obsIdList03FlankingKeys.contains(key1))
        assertTrue(obsIdList03FlankingKeys.contains(key2))

        // Note: this is where we switch from 1-based indexing to 0-based indexing
        val obsBedList01 = cr.generateBedRecords(obsIdList01)
        val obsBedList02 = cr.generateBedRecords(obsIdList01)

        assertEquals(2, obsBedList01.size)

        assertEquals(BedRecord("chr1", 34616,40204,"Zm00001eb000010",0,"."), obsBedList01[0] )
        assertEquals(BedRecord("chr1",34616,40204,"Zm00001eb000010",0,"."), obsBedList02[0] )

        // Verify a bad boundary type throws an error
        assertFails {
            cr.idMinMaxBounds(NucSeqIO(refGood).readAll(), genes, "geeeene")
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
        assertEquals(command.rangeMinSize, 500)
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
    fun testFileOutputWithOverlappingRegions() {
        val testGffPath = "src/test/resources/net/maizegenetics/phgv2/cli/zm_b73v5_test.gff3.gz"
        val ref = "src/test/resources/net/maizegenetics/phgv2/cli/RefChr.fa"
        val command = CreateRanges()

        val outputFileName = "${tempDir}test_extrapadding.bed"

        val result = command.test("--gff $testGffPath --reference-file $ref --pad 600 --output $outputFileName")
        assertEquals(result.statusCode, 0)
        assertEquals(command.gff, testGffPath)
        assertEquals(command.boundary, "gene")
        assertEquals(command.pad, 600)
        assertEquals(command.output, outputFileName)

        val lines = File(outputFileName).bufferedReader().readLines()
        assertEquals("chr1\t0\t34016\tintergenic_chr1:0-34016\t0\t+", lines[0])
        assertEquals("chr1\t34016\t40709\tZm00001eb000010\t0\t.", lines[1])
        assertEquals("chr1\t40709\t47362\tZm00001eb000020\t0\t.", lines[2])
        assertEquals("chr1\t47362\t55000\tintergenic_chr1:47362-55000\t0\t+", lines[3])
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

    @Test
    fun testCreateFlankingList() {

        // This test verifies the functionality of the addRange() and createFlankingList() functions
        // For addRange, it verifies that the ranges are added correctly, embedded genes are tossed and that overlapping
        // ranges are merged
        // For flankingList it verifies flanking handles going past the end of the chromosome, and that when there
        // are not enough bases to flank, the flanking is reduced to the available bases, split between 2 genes
        val refFile = "src/test/resources/net/maizegenetics/phgv2/cli/testRefGenome.fa"
        val geneRangeMap: RangeMap<Position, String> = TreeRangeMap.create()
        addRange(geneRangeMap, Range.closed(Position("chr1",442),Position("chr1",498)),"gene1")
        // These 2 genes overlap - they will be merged and the name will be gene2a-gene2b
        addRange(geneRangeMap, Range.closed(Position("chr1",708), Position("chr1",740)),"gene2a");
        addRange(geneRangeMap, Range.closed(Position("chr1",738), Position("chr1",757)),"gene2b");

        addRange(geneRangeMap, Range.closed(Position("chr1",922), Position("chr1",951)),"gene3");
        // This gene is embedded in the one gene above
        addRange(geneRangeMap, Range.closed(Position("chr1",930), Position("chr1",950)),"geneEmbedded");

        // ALlowing the next 2 to overlap when we flank with 100 bps
        addRange(geneRangeMap, Range.closed(Position("chr1",1216),
            Position("chr1",1283)),"gene4");
        addRange(geneRangeMap, Range.closed(Position("chr1",1486),
            Position("chr1",1497)),"gene5");
        addRange(geneRangeMap, Range.closed(Position("chr1",2221),
            Position("chr1",2226)),"gene6");
        addRange(geneRangeMap, Range.closed(Position("chr1",2483),
            Position("chr1",2490)),"gene7");
        addRange(geneRangeMap, Range.closed(Position("chr1",2869),
            Position("chr1",2891)),"gene8");
        addRange(geneRangeMap, Range.closed(Position("chr1",3003),
            Position("chr1",3057)),"gene9");

        // This range cannot be flanked by 100 as it hits the end of the ref chrom - will
        // verify the end is 4959 (the end of the chromosome) is 4959
        addRange(geneRangeMap, Range.closed(Position("chr1",4904),
            Position("chr1",4945)),"gene10");

        println("Last gene added, size of range: " + geneRangeMap.asMapOfRanges().size);
        // There were 12 genes added, but one was embedded, and 2 overlapped so we should only have 10
        assertEquals(geneRangeMap.asMapOfRanges().size,10);


        // verify map with lowerBound = 708 has value gene2a-gene2b
        var gene2 = Position("chr1",708);
        assertEquals(geneRangeMap.getEntry(gene2)!!.value,"gene2a-gene2b");

        // verify map with lowerBound = 922 has value gene3-geneEmbedded
        var gene3 = Position("chr1",922);
        assertEquals(geneRangeMap.getEntry(gene3)!!.value,"gene3-geneEmbedded");

        // Create ref sequence. Needed for createFlankingLIst to find the length of each chromosome
        // when adding flanking to last gene.
        println("Create genome sequence ...");
        val refSeq = NucSeqIO(refFile).readAll();
        val chromLen = refSeq["chr1"]!!.size()

        println("Creating gene-with-flanking map");

        val flankingGeneMap:RangeMap<Position,String> = createFlankingList( geneRangeMap,  100, refSeq);

        assertEquals(flankingGeneMap.asMapOfRanges().size,10);

        // gene1 - check range based on original lower bound
        // The getEntry call returns the range containing this position
        var  gene = Position("chr1",442);

        var geneLowerFlank = flankingGeneMap.getEntry(gene)!!.key.lowerEndpoint().position;
        var geneUpperFlank = flankingGeneMap.getEntry(gene)!!.key.upperEndpoint().position;
        assertEquals(342, geneLowerFlank);
        assertEquals(598, geneUpperFlank);

        //gene2
        // there is not 200 bps between the end of gene2 and start of gene3, so the code
        // splits the difference when flanking, which gets the end here as 839 and the
        // beginning of gene 3 (flanked) as 840)
        gene = Position("chr1",708);
        geneLowerFlank = flankingGeneMap.getEntry(gene)!!.key.lowerEndpoint().position;
        geneUpperFlank = flankingGeneMap.getEntry(gene)!!.key.upperEndpoint().position;
        assertEquals(608,geneLowerFlank);
        assertEquals(839,geneUpperFlank);

        //gene3
        gene = Position("chr1",922);
        geneLowerFlank = flankingGeneMap.getEntry(gene)!!.key.lowerEndpoint().position;
        geneUpperFlank = flankingGeneMap.getEntry(gene)!!.key.upperEndpoint().position;
        assertEquals(840,geneLowerFlank);
        assertEquals(1051, geneUpperFlank);

        // gene4
        gene = Position("chr1",1216);
        geneLowerFlank = flankingGeneMap.getEntry(gene)!!.key.lowerEndpoint().position;
        geneUpperFlank = flankingGeneMap.getEntry(gene)!!.key.upperEndpoint().position;
        assertEquals(1116,geneLowerFlank);
        assertEquals(1383, geneUpperFlank);

        // gene5
        gene = Position("chr1",1486);
        geneLowerFlank = flankingGeneMap.getEntry(gene)!!.key.lowerEndpoint().position;
        geneUpperFlank = flankingGeneMap.getEntry(gene)!!.key.upperEndpoint().position;
        assertEquals(1386,geneLowerFlank);
        assertEquals(1597,geneUpperFlank);

        // gene6 - checking with original upperbound now
        gene = Position("chr1",2226);
        geneLowerFlank = flankingGeneMap.getEntry(gene)!!.key.lowerEndpoint().position;
        geneUpperFlank = flankingGeneMap.getEntry(gene)!!.key.upperEndpoint().position;
        assertEquals(2121,geneLowerFlank);
        assertEquals(2326,geneUpperFlank);

        // gene7 - checking with original upperbound now
        gene = Position("chr1",2490);
        geneLowerFlank = flankingGeneMap.getEntry(gene)!!.key.lowerEndpoint().position;
        geneUpperFlank = flankingGeneMap.getEntry(gene)!!.key.upperEndpoint().position;
        assertEquals(2383,geneLowerFlank);
        assertEquals(2590,geneUpperFlank);

        // gene8 - checking with original upperbound now
        gene = Position("chr1",2891);
        geneLowerFlank = flankingGeneMap.getEntry(gene)!!.key.lowerEndpoint().position;
        geneUpperFlank = flankingGeneMap.getEntry(gene)!!.key.upperEndpoint().position;
        assertEquals(2769,geneLowerFlank);
        assertEquals(2947,geneUpperFlank);

        // gene9 - checking with original upperbound now
        // another overlap.  splitting the difference between the previous gene and this one,
        // there was only 112 bps between the two
        gene = Position("chr1",3057);
        geneLowerFlank = flankingGeneMap.getEntry(gene)!!.key.lowerEndpoint().position;
        geneUpperFlank = flankingGeneMap.getEntry(gene)!!.key.upperEndpoint().position;
        assertEquals(2948,geneLowerFlank);
        assertEquals(3157,geneUpperFlank);

        // gene10 - checking with a value from the middle of the original range
        // This range cannot be flanked by 100 as it hits the end of the ref chrom - will
        // verify the end is 4959 (the end of the chromosome) is 4959
        gene = Position("chr1",4950);
        geneLowerFlank = flankingGeneMap.getEntry(gene)!!.key.lowerEndpoint().position;
        geneUpperFlank = flankingGeneMap.getEntry(gene)!!.key.upperEndpoint().position;
        assertEquals(4804,geneLowerFlank);
        assertEquals(4959,geneUpperFlank);

        println("\nFinished testCreateFlankingList!\n");


    }


    @Test
    fun testMergeShortRangesUsuallyLeft() {
        // This test verifies the functionality of the mergeShortRanges() function.
        // Create a bed file with 16 ranges, some entries are less than 10 bps long
        // The function should merge the short ranges to the adjacent range which is
        //  shortest in length.  Merges "left" is adjacent ranges have equal length.
        // Verify if the first region is less than 10 bps, it is merged with the next region

        // create a list of BedRecord based on the beeFile coordinates below
        // ARe the BedRecord values 0-based or 1-based?  And inclusive or exclusive on the ends?
        val origBedList = mutableListOf<BedRecord>()
        origBedList.add(BedRecord("chr1", 0,10,"gene1",0,"+"))
        origBedList.add(BedRecord("chr1", 10,20,"gene2",0,"+"))
        origBedList.add(BedRecord("chr1", 20,30,"gene3",0,"+"))
        origBedList.add(BedRecord("chr1", 30,35,"gene4",0,"+")) // less than 10
        origBedList.add(BedRecord("chr1", 40,50,"gene5",0,"+"))
        origBedList.add(BedRecord("chr1", 50,60,"gene6",0,"+"))
        origBedList.add(BedRecord("chr1", 60,66,"gene7",0,"+")) // less than 10
        origBedList.add(BedRecord("chr1", 70,80,"gene8",0,"+"))
        origBedList.add(BedRecord("chr1", 80,90,"gene9",0,"+"))
        origBedList.add(BedRecord("chr1", 90,100,"gene10",0,"+"))
        origBedList.add(BedRecord("chr2", 0,3,"gene6",0,"+")) // less than 10
        origBedList.add(BedRecord("chr2", 30,60,"gene7",0,"+"))
        origBedList.add(BedRecord("chr2", 60,66,"gene8",0,"+")) // less than 10, but will merge right
        origBedList.add(BedRecord("chr2", 80,90,"gene9",0,"+"))
        origBedList.add(BedRecord("chr2", 90,100,"gene10",0,"+"))
        origBedList.add(BedRecord("chr2", 100,105,"gene11",0,"+")) // less than 10

        val createRanges = CreateRanges()
        // needs a bedRecords: List<BedRecord>), which is also what it will return
        val newBedLines = createRanges.mergeShortRanges(origBedList,10)

        println("newBedLines size: " + newBedLines.size)
        for (line in newBedLines) {
            println( line.toString())
        }

        // 5 records are smaller than 10 bps so they should be merged into the previous record.
        // This includes the first record of chr2, which is only 3 bps long.  It should be merged
        // with the succeeding record.
        assertEquals(11, newBedLines.size)
        // Verify the 4th record (chr1:30-35) is merged with the 3rd, which is now chr1:20-35
        assertEquals(BedRecord("chr1", 20,35,"gene3,gene4",0,"+"), newBedLines[2])
        // Verify the 6th record(chr1:60-66) is merged with the 5th, which is now chr1:50-66
        assertEquals(BedRecord("chr1", 50,66,"gene6,gene7",0,"+"), newBedLines[4])
        // For chrom 2, the first and 3rd entries are both short, resulting in the first 3 entries being merged.
        // This means the 9th record of the new bed file is chr2:0-66
        assertEquals(BedRecord("chr2", 0,60,"gene6,gene7",0,"+"), newBedLines[8])
        // Because the first 2 regions in chr2 were merged, the original region 3, which is also shorter than 10bps,
        // will merge with original region 4.  Merging "right" vs "left" because the "left" region now has size of 66
        // which is larger than the "right" region which has size of 10.
        assertEquals(BedRecord("chr2", 60,90,"gene8,gene9",0,"+"), newBedLines[9])
        // verify the last record in the new bedfile is a merge of the last 2 records of the original bed file
        assertEquals(BedRecord("chr2", 90,105,"gene10,gene11",0,"+"), newBedLines[10])

    }

    @Test
    fun testMergeConsecutiveShortRanges() {
        // Test many short ranges in successions.  First on chr2 is short
        // and must be merged with the next region.

        val origBedList = mutableListOf<BedRecord>()
        origBedList.add(BedRecord("chr1", 0,150,"gene1",0,"+"))
        origBedList.add(BedRecord("chr1", 150,250,"gene2",0,"+"))
        origBedList.add(BedRecord("chr1", 250,255,"gene3",0,"+"))
        origBedList.add(BedRecord("chr1", 255,260,"gene4",0,"+"))
        origBedList.add(BedRecord("chr1", 260,265,"gene5",0,"+"))
        origBedList.add(BedRecord("chr1", 365, 465,"gene6",0,"+"))
        origBedList.add(BedRecord("chr1", 465,660,"gene7",0,"+"))
        origBedList.add(BedRecord("chr2", 0,3,"gene6",0,"+"))
        origBedList.add(BedRecord("chr2", 30,60,"gene7",0,"+"))
        origBedList.add(BedRecord("chr2", 60,66,"gene8",0,"+"))
        origBedList.add(BedRecord("chr2", 80,90,"gene9",0,"+"))
        origBedList.add(BedRecord("chr2", 90,100,"gene10",0,"+"))
        origBedList.add(BedRecord("chr2", 100,105,"gene11",0,"+"))

        val createRanges = CreateRanges()
        // needs a bedRecords: List<BedRecord>), which is also what it will return
        val newBedLines = createRanges.mergeShortRanges(origBedList,50)

        println("newBedLines size: " + newBedLines.size)
        for (line in newBedLines) {
            println( line.toString())
        }
        // verify the new bed file has 6 records
        assertEquals(6, newBedLines.size)
        // verify the first record of the newBedLines is the same as the first record of the original bed file
        assertEquals(BedRecord("chr1", 0,150,"gene1",0,"+"), newBedLines[0])
        // verify the first bedrecord for chr2 is merged with the second record
        assertEquals(BedRecord("chr2", 0,60,"gene6,gene7",0,"+"), newBedLines[4])
        // verify the last 4 records of the original bed file are merged into 1 record
        assertEquals(BedRecord("chr2", 60,105,"gene8,gene9,gene10,gene11",0,"+"), newBedLines[5])


    }

    @Test
    fun testMergeRightThenLeft() {
        val origBedList = mutableListOf<BedRecord>()
        // For chr3, the 3rd=5th regions are small and shoudl be merged together. The resulting region is
        // still small and should be merged with the second region as it is smaller than the 6th region
        // with the previous record as that is smaller than the
        origBedList.add(BedRecord("chr3", 0,150,"gene1",0,"+"))
        origBedList.add(BedRecord("chr3", 150,250,"gene2",0,"+"))
        origBedList.add(BedRecord("chr3", 250,255,"gene3",0,"+"))
        origBedList.add(BedRecord("chr3", 255,260,"gene4",0,"+"))
        origBedList.add(BedRecord("chr3", 260,265,"gene5",0,"+"))
        origBedList.add(BedRecord("chr3", 365, 665,"gene6",0,"+"))
        origBedList.add(BedRecord("chr3", 665,750,"gene7",0,"+"))

        val createRanges = CreateRanges()
        // needs a bedRecords: List<BedRecord>), which is also what it will return
        val newBedLines = createRanges.mergeShortRanges(origBedList,50)

        println("newBedLines size: " + newBedLines.size)
        for (line in newBedLines) {
            println( line.toString())
        }
        // verify the new bed file has 4 records
        assertEquals(4, newBedLines.size)
        // verify the 2nd, 3rd, 4th and 5th records are merged into 1 record
        assertEquals(BedRecord("chr3", 150,265,"gene2,gene3,gene4,gene5",0,"+"), newBedLines[1])
        // verify the first record of the newBedLines is the same as the first record of the original bed file
        assertEquals(BedRecord("chr3", 0,150,"gene1",0,"+"), newBedLines[0])
        // verify the last 2 records in the newBedLines are equal to the last 2 records in the original bed file
        assertEquals(BedRecord("chr3", 365, 665,"gene6",0,"+"), newBedLines[2])
        // verify the last record in the new bedfile is a merge of the last 2 records of the original bed file
        assertEquals(BedRecord("chr3", 665,750,"gene7",0,"+"), newBedLines[3])
    }

    @Test
    fun testLastRegionIsShort() {
        val origBedList = mutableListOf<BedRecord>()
        // For chr3, the 3rd=5th regions are small and shoudl be merged together. The resulting region is
        // still small and should be merged with the second region as it is smaller than the 6th region
        // with the previous record as that is smaller than the
        origBedList.add(BedRecord("chr3", 0,150,"gene1",0,"+"))
        origBedList.add(BedRecord("chr3", 150,250,"gene2",0,"+"))
        origBedList.add(BedRecord("chr3", 250,365,"gene3",0,"+"))
        origBedList.add(BedRecord("chr3", 365, 370,"gene6",0,"+"))

        val createRanges = CreateRanges()
        // needs a bedRecords: List<BedRecord>), which is also what it will return
        val newBedLines = createRanges.mergeShortRanges(origBedList,50)

        println("newBedLines size: " + newBedLines.size)
        for (line in newBedLines) {
            println( line.toString())
        }
        // verify the last record in the new bedfile is a merge of the last 2 records of the original bed file
        assertEquals(BedRecord("chr3", 250,370,"gene3,gene6",0,"+"), newBedLines[2])
        // verify newBedLines has 3 records
        assertEquals(3, newBedLines.size)
    }

}
