package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedReader
import htsjdk.samtools.*
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.assertThrows
import java.io.BufferedReader
import java.io.File
import java.io.FileReader
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class ExtractEdgeReadsTest {

    private fun createSAMRecord(factory: SAMRecordFactory, header: SAMFileHeader, readName : String,
                                readString : String, mappingQuality : Int, contig : String, BS : String,
                                startPos: Int = 1, firstInPair: Boolean = true, cigarString: String = ""): SAMRecord {
        val samRecord = factory.createSAMRecord(header)
        samRecord.readName = readName
        samRecord.readString = readString
        samRecord.mappingQuality = mappingQuality
        samRecord.referenceName = contig
        samRecord.baseQualityString = BS
        samRecord.alignmentStart = startPos
        if(firstInPair) {
            samRecord.firstOfPairFlag = true
        } else {
            samRecord.secondOfPairFlag = true
        }
        if(cigarString == "") samRecord.cigar = Cigar.fromCigarString(cigarString)
        return samRecord
    }

    @Test
    fun testOutputReadsAndHapIdSets() {
        // Create a SAMFileHeader
        val samHeader = SAMFileHeader()
        // Create a SAMRecordFactory to generate SAM records
        val samRecordFactory = DefaultSAMRecordFactory()
        // Generate four SAMRecords
        val sam1 = createSAMRecord(samRecordFactory, samHeader, "read1", "AAAA", 30, "chr1", "IIII")
        val sam2 = createSAMRecord(samRecordFactory, samHeader, "read2", "TTTT", 40, "chr2", "IIII")
        val sam3 = createSAMRecord(samRecordFactory, samHeader, "read3", "CCCC", 0, "chr3", "IIII")
        val sam4 = createSAMRecord(samRecordFactory, samHeader, "read4", "GGGG", 50, "chr4", "IIII")

        // construct input data
        val list : List<Pair<SAMRecord, SAMRecord>> = listOf(Pair(sam1, sam2), Pair(sam3, sam4))
        val alignments : Map<String, List<Pair<SAMRecord, SAMRecord>>> = mapOf("readID1" to list)

        val extractEdgeReads = ExtractEdgeReads()
        // create output files
        extractEdgeReads.outputReadsAndHapIdSets(tableFileName = "table.txt", fastqFileName = "fastq.fq", alignments)

        // read in data from output files and check if correct
        bufferedReader("table.txt").use { reader ->
            val header : String = reader.readLine()
            assertEquals(header, "readID\tHapIDHits")
            val line1 : String = reader.readLine()
            assertEquals(line1, "readID1\tchr1, chr2, chr3, chr4")
        }
        bufferedReader("1fastq.fq").use { reader ->
            val name1 : String = reader.readLine()
            assertEquals(name1, "@read1")
            val seq1 : String = reader.readLine()
            assertEquals(seq1, "AAAA")
            val plus1 : String = reader.readLine()
            assertEquals(plus1, "+")
            val score1 : String = reader.readLine()
            assertEquals(score1, "IIII")

            val name2 : String = reader.readLine()
            assertEquals(name2, "@read3")
            val seq2 : String = reader.readLine()
            assertEquals(seq2, "CCCC")
            val plus2 : String = reader.readLine()
            assertEquals(plus2, "+")
            val score2 : String = reader.readLine()
            assertEquals(score2, "IIII")
        }
        bufferedReader("2fastq.fq").use { reader ->

            val name1 : String = reader.readLine()
            assertEquals(name1, "@read2")
            val seq1 : String = reader.readLine()
            assertEquals(seq1, "TTTT")
            val plus1 : String = reader.readLine()
            assertEquals(plus1, "+")
            val score1 : String = reader.readLine()
            assertEquals(score1, "IIII")

            val name2 : String = reader.readLine()
            assertEquals(name2, "@read4")
            val seq2 : String = reader.readLine()
            assertEquals(seq2, "GGGG")
            val plus2 : String = reader.readLine()
            assertEquals(plus2, "+")
            val score2 : String = reader.readLine()
            assertEquals(score2, "IIII")
        }
        File("table.txt").delete()
        File("1fastq.fq").delete()
        File("2fastq.fq").delete()
    }

    @Test
    fun testIsSingleAndPaired() {
        val extractEdgeReads = ExtractEdgeReads()
        val nullRecords = listOf<Pair<SAMRecord?, SAMRecord?>>(Pair(null, null))
        //In reality this should not happen at all as unaligned is checked first
        assertTrue(extractEdgeReads.isSingle(nullRecords), "Null Records are not passing isSingle correctly")

        val pairedRecords = mutableListOf<Pair<SAMRecord?,SAMRecord?>>()
        val samRecordFactory = DefaultSAMRecordFactory()
        val samHeader = SAMFileHeader()

        val samRecord1 = createSAMRecord(samRecordFactory, samHeader, "read1", "AAAAAAAAAA", 60,"hap1", "IIIIIIIIII",1, true,"10M", )

        val samRecord2 = createSAMRecord(samRecordFactory, samHeader, "read1", "AAAAAAAAAA", 60,"hap1",  "IIIIIIIIII", 50, false, "10M")

        pairedRecords.add(Pair(samRecord1, samRecord2))

        val samRecord3 = createSAMRecord(samRecordFactory, samHeader, "read2", "AAAAAAAAAA", 60, "hap2", "IIIIIIIIII", 1, true, "10M")

        val samRecord4 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA", 60, "hap2", "IIIIIIIIII", 50, false, "10M")

        pairedRecords.add(Pair(samRecord3, samRecord4))

        assertFalse(extractEdgeReads.isSingle(pairedRecords), "Fully paired records are being called as single")
        assertTrue(extractEdgeReads.isPaired(pairedRecords), "Fully paired records are not being called as paired")

        pairedRecords.add(Pair(samRecord1, null))
        assertFalse(extractEdgeReads.isSingle(pairedRecords), "Partially paired records are being called as single")
        assertTrue(extractEdgeReads.isPaired(pairedRecords), "Partially paired records are not being called as paired")

        //Add in a single ended one to the list to verify that it is still failing the test
        val singleSamRecords = mutableListOf<Pair<SAMRecord?,SAMRecord?>>()
        singleSamRecords.add(Pair(samRecord1, null))
        singleSamRecords.add(Pair(null, samRecord4))
        assertTrue(extractEdgeReads.isSingle(singleSamRecords), "Single ended records are not being called as single")
        assertFalse(extractEdgeReads.isPaired(singleSamRecords), "Single ended records are being called as paired")
    }

    @Test
    fun testIsSingleAlignSplit() {
        val extractEdgeReads = ExtractEdgeReads()
        val samRecordFactory = DefaultSAMRecordFactory()
        val samHeader = SAMFileHeader()

        //Test null first  Should be false as the set will be empty
        assertFalse(extractEdgeReads.isSingleAlignSplit(listOf(), mapOf()), "Empty list is being called as split")

        //Assuming we only have single reads
        val samRecord1 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, true, "10M")
        val samRecord3 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")

        val records = listOf(Pair(samRecord1, null), Pair(samRecord2, null), Pair(null,samRecord3))

        val refRange1 = ReferenceRange("1", 1, 10)
        val hapIdToRefRangeSameRefRange = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange1)
        )

        assertFalse(extractEdgeReads.isSingleAlignSplit(records, hapIdToRefRangeSameRefRange), "Single align split is not being called correctly")

        val refRange2 = ReferenceRange("2", 1, 10)
        val hapIdToRefRangeDiffRefRange = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange2)
        )
        assertTrue(extractEdgeReads.isSingleAlignSplit(records, hapIdToRefRangeDiffRefRange), "Single align split is not being called correctly")

        //Hap Shared across refRanges
        val hapIdToRefRangeSharedHap = mapOf("hap1" to listOf(refRange1, refRange2),
            "hap2" to listOf(refRange1)
        )

        //This is align split as some of the haplotypes are shared across ref ranges
        assertTrue(extractEdgeReads.isSingleAlignSplit(records, hapIdToRefRangeSharedHap), "Single align split is not being called correctly")
    }

    @Test
    fun testIsSingleAlignConsec() {
        val extractEdgeReads = ExtractEdgeReads()
        val samRecordFactory = DefaultSAMRecordFactory()
        val samHeader = SAMFileHeader()

        //Test null first  Should be false as the set will be empty
        assertFalse(extractEdgeReads.isSingleAlignConsec(listOf(), mapOf(),mapOf()), "Empty list is being called as split")

        //Assuming we only have single reads
        val samRecord1 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, true, "10M")
        val samRecord3 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")

        val records = listOf(Pair(samRecord1, null), Pair(samRecord2, null), Pair(null,samRecord3))

        val refRange1 = ReferenceRange("1", 1, 10)
        val refRange2 = ReferenceRange("1", 11, 20)
        val refRange3 = ReferenceRange("2", 1, 10)
        //refRangeToIndexMap: Map<String,Int>
        val refRangeToIndexMap = mapOf(refRange1.toString() to 0, refRange2.toString() to 1, refRange3.toString() to 2)

        val hapIdToRefRangeSameRefRange = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange1)
        )

        assertFalse(extractEdgeReads.isSingleAlignConsec(records, hapIdToRefRangeSameRefRange, refRangeToIndexMap), "Single align split is not being called correctly")


        val hapIdToRefRangeDiffRefRange = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange2)
        )
        assertTrue(extractEdgeReads.isSingleAlignConsec(records, hapIdToRefRangeDiffRefRange, refRangeToIndexMap), "Single align split is not being called correctly")

        //Hap Shared across refRanges
        val hapIdToRefRangeSharedHap = mapOf("hap1" to listOf(refRange1, refRange2),
            "hap2" to listOf(refRange3)
        )

        assertFalse(extractEdgeReads.isSingleAlignConsec(records, hapIdToRefRangeSharedHap, refRangeToIndexMap), "Single align split is not being called correctly")


        //Reverse the order of the records
        val hapIdToRefRangeReverse = mapOf("hap1" to listOf(refRange2),
            "hap2" to listOf(refRange1)
        )
        assertTrue(extractEdgeReads.isSingleAlignConsec(records, hapIdToRefRangeReverse, refRangeToIndexMap), "Single align split is not being called correctly")
    }

    @Test
    fun testIsSingleOffASM() {
        val extractEdgeReads = ExtractEdgeReads()
        val samRecordFactory = DefaultSAMRecordFactory()
        val samHeader = SAMFileHeader()

        //Test null first  Should be false as the set will be empty
        assertThrows<IllegalStateException> { extractEdgeReads.isSingleOffASM("sample1", listOf(), mapOf())}

        //Assuming we only have single reads
        val samRecord1 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, true, "10M")
        val samRecord3 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")


        val records = listOf(Pair(samRecord1, null), Pair(samRecord2, null), Pair(null,samRecord3))

        val hapIdToSampleGamete = mapOf("hap1" to listOf(SampleGamete("sample1", 0)),
            "hap2" to listOf(SampleGamete("sample1", 1))
        )

        assertFalse(extractEdgeReads.isSingleOffASM("sample1", records, hapIdToSampleGamete), "Single align split is not being called correctly")

        val hapIdToSampleGameteOffASM = mapOf("hap1" to listOf(SampleGamete("sample1", 0)),
            "hap2" to listOf(SampleGamete("sample2", 0))
        )
        assertTrue(extractEdgeReads.isSingleOffASM("sample3", records, hapIdToSampleGameteOffASM), "Single align split is not being called correctly")
    }


    //isPairedReadSplit
    @Test
    fun testIsPairedAlignConsec() {
        val extractEdgeReads = ExtractEdgeReads()
        val samRecordFactory = DefaultSAMRecordFactory()
        val samHeader = SAMFileHeader()

        //Test null first  Should be false as the set will be empty
        assertFalse(extractEdgeReads.isPairedAlignConsec(listOf(), mapOf(),mapOf()), "Empty list is being called as split")

        //Assuming we only have single reads
        val samRecord1_1 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord1_2 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")
        val samRecord2_1 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, true, "10M")
        val samRecord2_2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, false, "10M")
        val samRecord3_1 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord3_2 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")

        val records = listOf(Pair(samRecord1_1, samRecord1_2), Pair(samRecord2_1, samRecord2_2), Pair(samRecord3_1,samRecord3_2))

        val refRange1 = ReferenceRange("1", 1, 10)
        val refRange2 = ReferenceRange("1", 11, 20)
        val refRange3 = ReferenceRange("2", 1, 10)
        val refRangeToIndexMap = mapOf(refRange1.toString() to 0, refRange2.toString() to 1, refRange3.toString() to 2)

        val hapIdToRefRangeSameRefRange = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange1)
        )

        assertFalse(extractEdgeReads.isPairedAlignConsec(records, hapIdToRefRangeSameRefRange, refRangeToIndexMap), "Paired align split is not being called correctly")

        val hapIdToRefRangeDiffRefRange = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange2)
        )
        assertTrue(extractEdgeReads.isPairedAlignConsec(records, hapIdToRefRangeDiffRefRange, refRangeToIndexMap), "Paired align split is not being called correctly")

        //Hap Shared across refRanges
        val hapIdToRefRangeSharedHap = mapOf("hap1" to listOf(refRange1, refRange2),
            "hap2" to listOf(refRange3)
        )

        assertFalse(extractEdgeReads.isPairedAlignConsec(records, hapIdToRefRangeSharedHap, refRangeToIndexMap), "Paired align split is not being called correctly")

        val hapIdToRefRangeTooFar = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange3)
        )

        assertFalse(extractEdgeReads.isPairedAlignConsec(records, hapIdToRefRangeTooFar, refRangeToIndexMap), "Paired align split is not being called correctly")
    }

    @Test
    fun testIsPairedAlignSplit() {
        val extractEdgeReads = ExtractEdgeReads()
        val samRecordFactory = DefaultSAMRecordFactory()
        val samHeader = SAMFileHeader()

        //Test null first  Should be false as the set will be empty
        assertFalse(extractEdgeReads.isPairedAlignSplit(listOf(), mapOf()), "Empty list is being called as split")

        //Assuming we only have single reads
        val samRecord1_1 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord1_2 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")
        val samRecord2_1 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, true, "10M")
        val samRecord2_2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, false, "10M")
        val samRecord3_1 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord3_2 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")

        val records = listOf(Pair(samRecord1_1, samRecord1_2), Pair(samRecord2_1, samRecord2_2), Pair(samRecord3_1,samRecord3_2))

        val refRange1 = ReferenceRange("1", 1, 10)
        val refRange2 = ReferenceRange("1", 11, 20)
        val refRange3 = ReferenceRange("2", 1, 10)

        val hapIdToRefRangeSameRefRange = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange1)
        )

        assertFalse(extractEdgeReads.isPairedAlignSplit(records, hapIdToRefRangeSameRefRange), "Paired align split is not being called correctly")

        val hapIdToRefRangeDiffRefRange = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange2)
        )
        assertTrue(extractEdgeReads.isPairedAlignSplit(records, hapIdToRefRangeDiffRefRange), "Paired align split is not being called correctly")

        //Hap Shared across refRanges
        val hapIdToRefRangeSharedHap = mapOf("hap1" to listOf(refRange1, refRange2),
            "hap2" to listOf(refRange3)
        )

        assertTrue(extractEdgeReads.isPairedAlignSplit(records, hapIdToRefRangeSharedHap), "Paired align split is not being called correctly")

        val hapIdToRefRangeTooFar = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange3)
        )

        assertTrue(extractEdgeReads.isPairedAlignSplit(records, hapIdToRefRangeTooFar), "Paired align split is not being called correctly")

    }

    @Test
    fun testIsPairedOffASM() {
        val extractEdgeReads = ExtractEdgeReads()
        val samRecordFactory = DefaultSAMRecordFactory()
        val samHeader = SAMFileHeader()

        //Test null first  Should be false as the set will be empty
        assertThrows<IllegalStateException> { extractEdgeReads.isPairedOffASM("sample1", listOf(), mapOf())}

        //Assuming we only have single reads
        val samRecord1_1 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord1_2 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")
        val samRecord2_1 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, true, "10M")
        val samRecord2_2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, false, "10M")
        val samRecord3_1 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord3_2 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")

        val records = listOf(Pair(samRecord1_1, samRecord1_2), Pair(samRecord2_1, samRecord2_2), Pair(samRecord3_1,samRecord3_2))

        val hapIdToSampleGamete = mapOf("hap1" to listOf(SampleGamete("sample1", 0)),
            "hap2" to listOf(SampleGamete("sample1", 1))
        )

        assertFalse(extractEdgeReads.isPairedOffASM("sample1", records, hapIdToSampleGamete), "Single align split is not being called correctly")

        val hapIdToSampleGameteOffASM = mapOf("hap1" to listOf(SampleGamete("sample1", 0)),
            "hap2" to listOf(SampleGamete("sample2", 0))
        )

        assertTrue(extractEdgeReads.isPairedOffASM("sample3", records, hapIdToSampleGameteOffASM), "Single align split is not being called correctly")

    }

    @Test
    fun testIsPairedReadSplitConsec() {
        val extractEdgeReads = ExtractEdgeReads()
        val samRecordFactory = DefaultSAMRecordFactory()
        val samHeader = SAMFileHeader()

        //Test null first  Should be false as the set will be empty
        assertFalse(extractEdgeReads.isPairedReadSplitConsec(listOf(), mapOf(),mapOf()), "Empty list is being called as split")

        //Assuming we only have single reads
        val samRecord1_1 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord1_2 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")
        val samRecord2_1 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, true, "10M")
        val samRecord2_2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, false, "10M")
        val samRecord3_1 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord3_2 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")

        val records = listOf(Pair(samRecord1_1, samRecord1_2), Pair(samRecord2_1, samRecord2_2), Pair(samRecord3_1,samRecord3_2))

        val refRange1 = ReferenceRange("1", 1, 10)
        val refRange2 = ReferenceRange("1", 11, 20)
        val refRange3 = ReferenceRange("2", 1, 10)

        val refRangeToIndexMap = mapOf(refRange1.toString() to 0, refRange2.toString() to 1, refRange3.toString() to 2)

        val hapIdToRefRangeSameRefRange = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange1)
        )

        assertFalse(extractEdgeReads.isPairedReadSplitConsec(records, hapIdToRefRangeSameRefRange, refRangeToIndexMap), "Paired align split is not being called correctly")

        val hapIdToRefRangeDiffRefRange = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange2)
        )

        assertFalse(extractEdgeReads.isPairedReadSplitConsec(records, hapIdToRefRangeDiffRefRange, refRangeToIndexMap), "Paired align split is not being called correctly")

        val records2 = listOf(Pair(samRecord1_1, samRecord1_2), Pair(samRecord2_1, samRecord1_2), Pair(samRecord3_1,samRecord3_2))
        assertTrue(extractEdgeReads.isPairedReadSplitConsec(records2, hapIdToRefRangeDiffRefRange, refRangeToIndexMap), "Paired align split is not being called correctly")

        //Hap Shared across refRanges
        val hapIdToRefRangeSharedHap = mapOf("hap1" to listOf(refRange1, refRange2),
            "hap2" to listOf(refRange3)
        )

        assertFalse(extractEdgeReads.isPairedReadSplitConsec(records, hapIdToRefRangeSharedHap, refRangeToIndexMap), "Paired align split is not being called correctly")

        val hapIdToRefRangeTooFar = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange3)
        )

        assertFalse(extractEdgeReads.isPairedReadSplitConsec(records, hapIdToRefRangeTooFar, refRangeToIndexMap), "Paired align split is not being called correctly")
    }

    @Test
    fun testIsPairedReadSplit() {
        val extractEdgeReads = ExtractEdgeReads()
        val samRecordFactory = DefaultSAMRecordFactory()
        val samHeader = SAMFileHeader()

        //Test null first  Should be false as the set will be empty
        assertFalse(extractEdgeReads.isPairedReadSplit(listOf(), mapOf()), "Empty list is being called as split")

        //Assuming we only have single reads
        val samRecord1_1 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord1_2 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")
        val samRecord2_1 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, true, "10M")
        val samRecord2_2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, false, "10M")
        val samRecord3_1 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord3_2 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")

        val records = listOf(Pair(samRecord1_1, samRecord1_2), Pair(samRecord2_1, samRecord2_2), Pair(samRecord3_1,samRecord3_2))

        val refRange1 = ReferenceRange("1", 1, 10)
        val refRange2 = ReferenceRange("1", 11, 20)
        val refRange3 = ReferenceRange("2", 1, 10)

        val hapIdToRefRangeSameRefRange = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange1)
        )

        assertFalse(extractEdgeReads.isPairedReadSplit(records, hapIdToRefRangeSameRefRange), "Paired align split is not being called correctly")

        val hapIdToRefRangeDiffRefRange = mapOf("hap1" to listOf(refRange1),
            "hap2" to listOf(refRange2)
        )

        assertFalse(extractEdgeReads.isPairedReadSplit(records, hapIdToRefRangeDiffRefRange), "Paired align split is not being called correctly")

        val records2 = listOf(Pair(samRecord1_1, samRecord1_2), Pair(samRecord2_1, samRecord1_2), Pair(samRecord3_1,samRecord3_2))
        assertTrue(extractEdgeReads.isPairedReadSplit(records2, hapIdToRefRangeDiffRefRange), "Paired align split is not being called correctly")

    }

    @Test
    fun testClassifyAlignments() {
        TODO("Implement this test")
    }

    @Test
    fun testClassifySingleAlignments() {
        TODO("Implement this test")
    }

    @Test
    fun testClassifyPairedAlignments() {
        TODO("Implement this test")
    }
}