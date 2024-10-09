package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.testing.test
import htsjdk.samtools.*
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.assertThrows
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertFalse
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
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
        TODO("Implement readToClass map")
        // create output files
        extractEdgeReads.outputReadsAndHapIdSets(tableFileName = "table.txt", fastqFileName = "fastq", alignments, mapOf())

        // read in data from output files and check if correct
        bufferedReader("table.txt").use { reader ->
            val header : String = reader.readLine()
            assertEquals(header, "readID\tHapIDHits")
            val line1 : String = reader.readLine()
            assertEquals(line1, "readID1\tchr1, chr2, chr3, chr4")
        }
        bufferedReader("fastq_0.fastq").use { reader ->
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
        bufferedReader("fastq_1.fastq").use { reader ->

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
        File("fastq_0.fastq").delete()
        File("fastq_1.fastq").delete()
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
    fun testPairOffAlignments() {
        val extractEdgeReads = ExtractEdgeReads()
        val factory = DefaultSAMRecordFactory()
        val header = SAMFileHeader()
        val sam1 = createSAMRecord(factory, header, "first", "AAAA", 30, "chr1", "IIII", firstInPair = true, startPos = 10)
        val sam2 = createSAMRecord(factory, header, "second", "TTTT", 30, "chr1", "IIII", firstInPair = false, startPos = 20)

        sam1.readPairedFlag = true
        sam1.mateReferenceName = "chr1"
        sam1.mateAlignmentStart = 20
        sam1.mateUnmappedFlag = false

        sam2.readPairedFlag = true
        sam2.mateReferenceName = "chr1"
        sam2.mateAlignmentStart = 10
        sam2.mateUnmappedFlag = false

        val bestAlignments = extractEdgeReads.pairOffAlignments(listOf(sam1, sam2))
        assertEquals(bestAlignments, Pair(sam1, sam2))
    }

    @Test
    fun testKeepBestAlignments() {
        val extractEdgeReads = ExtractEdgeReads()
        val factory = DefaultSAMRecordFactory()
        val header = SAMFileHeader()
        val sam1 = createSAMRecord(factory, header, "first", "AAAA", 30, "chr1", "IIII", firstInPair = true, startPos = 10)
        val sam2 = createSAMRecord(factory, header, "second", "TTTT", 30, "chr1", "IIII", firstInPair = false, startPos = 20)
        val sam3 = createSAMRecord(factory, header, "third", "TTTT", 30, "chr1", "IIII", firstInPair = false, startPos = 20)

        sam1.setAttribute("NM", 10)
        sam2.setAttribute("NM", 20)
        sam3.setAttribute("NM", 0)

        assertEquals(extractEdgeReads.keepBestAlignment(listOf(sam1, sam2)), sam1)
        assertEquals(extractEdgeReads.keepBestAlignment(listOf(sam2, sam3)), sam3)
        assertEquals(extractEdgeReads.keepBestAlignment(listOf(sam1, sam3)), sam3)
    }

    @Test
    fun testClassifySingleAlignments() {
        val extractEdgeReads = ExtractEdgeReads()
        val samRecordFactory = DefaultSAMRecordFactory()
        val samHeader = SAMFileHeader()

        // test SINGLEUNIQUE
        val sam1 = createSAMRecord(samRecordFactory, samHeader, "read1", "AAAA", 30, "hap1", "IIII")
        val ref1 = ReferenceRange("1", 1, 10)
        val gam1 = SampleGamete("gam1", 1)
        val singleUnique = extractEdgeReads.classifySingleAlignments("gam1", 1,
            listOf(Pair(sam1, null)), mapOf("hap1" to listOf(ref1)), mapOf("hap1" to listOf(gam1)), mapOf("hap1" to 0))
        assertEquals(AlignmentClass.SINGLEUNIQUE, singleUnique, "fail SINGLEUNIQUE")

        // SINGLEALIGNSPLIT
        val samRecord1 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, true, "10M")
        val samRecord3 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")
        val records = listOf(Pair(samRecord1, null), Pair(samRecord2, null), Pair(null,samRecord3))
        val refRange1 = ReferenceRange("1", 1, 10)
        val refRange2 = ReferenceRange("2", 1, 10)
        val refRangeMap = mapOf("hap1" to listOf(refRange1), "hap2" to listOf(refRange2))
        val gameteMap = mapOf("hap1" to listOf(SampleGamete("1", 1)),
            "hap2" to listOf(SampleGamete("2", 2)))
        val indexMap = mapOf(refRange1.toString() to 0, refRange2.toString() to 2) // non consecutive ref ranges
        val singleAlignSplit = extractEdgeReads.classifySingleAlignments("1", 2,
            records, refRangeMap, gameteMap, indexMap)
        assertEquals(AlignmentClass.SINGLEALIGNSPLIT, singleAlignSplit, "fail SINGLEALIGNSPLIT")

        // SINGLEALIGNSPLITCONSEC
        val refRange4 = ReferenceRange("1", 11, 20)
        val refRange3 = ReferenceRange("2", 1, 10)
        val refRangeToIndexMap = mapOf(refRange1.toString() to 0, refRange4.toString() to 1, refRange3.toString() to 2) // consecutive
        val hapIdToRefRange = mapOf("hap1" to listOf(refRange1), "hap2" to listOf(refRange4))
        val singleAlignSplitConsec = extractEdgeReads.classifySingleAlignments("1",
            2, records, hapIdToRefRange, mapOf("hap1" to listOf(SampleGamete("1", 1)),
                "hap2" to listOf(SampleGamete("2", 2))), refRangeToIndexMap)
        assertEquals(AlignmentClass.SINGLEALIGNSPLITCONSEC, singleAlignSplitConsec, "fail SINGLEALIGNSPLITCONSEC")

        // SINGLEOFFASM
        val hapIdToSampleGameteOffASM = mapOf("hap1" to listOf(SampleGamete("sample1", 0)),
            "hap2" to listOf(SampleGamete("sample2", 0)), "hap3" to listOf(SampleGamete("sample3", 0)))
        val recordsOffASM = listOf(Pair(samRecord1, null), Pair(samRecord2, null))
        val refRangeMapOffASM = mapOf("hap1" to listOf(refRange1), "hap2" to listOf(refRange1))
        val indexOffASM = mapOf(refRange1.toString() to 0, refRange2.toString() to 0)
        val singleOffASM = extractEdgeReads.classifySingleAlignments("sample3", 2,
            recordsOffASM, refRangeMapOffASM, hapIdToSampleGameteOffASM, indexOffASM)
        assertEquals(AlignmentClass.SINGLEOFFASM, singleOffASM, "fail SINGLEOFFASM")

        // SINGLERARE
        val rareRecords = listOf(Pair(samRecord1, null), Pair(samRecord2, null), Pair(null,samRecord3))
        val rareRefRanges = mapOf("hap1" to listOf(refRange1), "hap2" to listOf(refRange1), "hap3" to listOf(refRange1))
        val rareSampleGametes = mapOf("hap1" to listOf(SampleGamete("rare1", 1)),
            "hap2" to listOf(SampleGamete("rare1", 1)), "hap3" to listOf(SampleGamete("rare1", 1)))
        val rareIndexMap = mapOf(refRange1.toString() to 0, refRange2.toString() to 0, refRange3.toString() to 0)
        val singleRare = extractEdgeReads.classifySingleAlignments("rare1", 100,
            rareRecords, rareRefRanges, rareSampleGametes, rareIndexMap)
        assertEquals(AlignmentClass.SINGLERARE, singleRare, "fail SINGLERARE")

        // SINGLECOMMON
        val singleCommon = extractEdgeReads.classifySingleAlignments("rare1", 4,
            rareRecords, rareRefRanges, rareSampleGametes, rareIndexMap)
        assertEquals(AlignmentClass.SINGLECOMMON, singleCommon, "fail SINGLECOMMON")
    }

    @Test
    fun testClassifyPairedAlignments() {
        val extractEdgeReads = ExtractEdgeReads()
        val samRecordFactory = DefaultSAMRecordFactory()
        val samHeader = SAMFileHeader()

        // PAIRUNIQUE
        val samRecord1_1 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord1_2 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")
        val uniqueRecords = listOf(Pair(samRecord1_1, samRecord1_2))
        val refRange1 = ReferenceRange("1", 1, 10)
        val pairUnique = extractEdgeReads.classifyPairedAlignments("sample1", 2, uniqueRecords,
            mapOf("hap1" to listOf(refRange1)), mapOf("hap1" to listOf(SampleGamete("sample1", 0))), mapOf(refRange1.toString() to 0))
        assertEquals(AlignmentClass.PAIRUNIQUE, pairUnique, "fail PAIRUNIQUE")

        // PAIRRARE
        val samRecord2_1 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, true, "10M")
        val samRecord2_2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, false, "10M")
        val samRecord3_1 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord3_2 = createSAMRecord(samRecordFactory, samHeader, "read3",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")
        val rareRecords = listOf(Pair(samRecord1_1, samRecord1_2), Pair(samRecord2_1, samRecord2_2), Pair(samRecord3_1, samRecord3_2))
        val rareRefRanges = mapOf("hap1" to listOf(refRange1), "hap2" to listOf(refRange1))
        val rareSampleGametes = mapOf("hap1" to listOf(SampleGamete("sample1", 1)),
            "hap2" to listOf(SampleGamete("sample1", 1)))
        val rareIndexMap = mapOf(refRange1.toString() to 0)
        val pairRare = extractEdgeReads.classifyPairedAlignments("sample1", 100, rareRecords, rareRefRanges, rareSampleGametes, rareIndexMap)
        assertEquals(AlignmentClass.PAIRRARE, pairRare, "fail PAIRRARE")

        // PAIRCOMMON
        val pairCommon = extractEdgeReads.classifyPairedAlignments("sample1", 3,
            rareRecords, rareRefRanges, rareSampleGametes, rareIndexMap)
        assertEquals(AlignmentClass.PAIRCOMMON, pairCommon, "fail PAIRCOMMON")

        // PAIRALIGNSPLIT
        val records = listOf(Pair(samRecord1_1, samRecord1_2), Pair(samRecord2_1, samRecord2_2), Pair(samRecord3_1, samRecord3_2))
        val refRange2 = ReferenceRange("2", 1, 10)
        val refRangeMap = mapOf("hap1" to listOf(refRange1), "hap2" to listOf(refRange2))
        val gameteMap = mapOf("hap1" to listOf(SampleGamete("sample1", 1)),
            "hap2" to listOf(SampleGamete("sample1", 1)))
        val nonconsecSplit = mapOf(refRange1.toString() to 0, refRange2.toString() to 2) // non consecutive ref ranges
        val pairAlignSplit = extractEdgeReads.classifyPairedAlignments("sample1", 2,
            records, refRangeMap, gameteMap, nonconsecSplit)
        assertEquals(AlignmentClass.PAIRALIGNSPLIT, pairAlignSplit, "fail PAIRALIGNSPLIT")

        // PAIRALIGNSPLITCONSEC
        val consecSplit = mapOf(refRange1.toString() to 0, refRange2.toString() to 1) // consecutive ref ranges
        val pairAlignSplitConsec = extractEdgeReads.classifyPairedAlignments("sample1", 2,
            records, refRangeMap, gameteMap, consecSplit)
        assertEquals(AlignmentClass.PAIRALIGNSPLITCONSEC, pairAlignSplitConsec, "fail PAIRALIGNSPLITCONSEC")

        // PAIROFFASM
        val hapIdToSampleGameteOffASM = mapOf("hap1" to listOf(SampleGamete("sample1", 0)),
            "hap2" to listOf(SampleGamete("sample2", 0)), "hap3" to listOf(SampleGamete("sample3", 0)))
        val recordsOffASM = listOf(Pair(samRecord1_1, samRecord1_2), Pair(samRecord2_1, samRecord2_2))
        val refRangeMapOffASM = mapOf("hap1" to listOf(refRange1), "hap2" to listOf(refRange1))
        val indexOffASM = mapOf(refRange1.toString() to 0, refRange2.toString() to 0)
        val pairOffASM = extractEdgeReads.classifyPairedAlignments("sample3", 2,
            recordsOffASM, refRangeMapOffASM, hapIdToSampleGameteOffASM, indexOffASM)
        assertEquals(AlignmentClass.PAIROFFASM, pairOffASM, "fail PAIROFFASM")

        // PAIRREADSPLIT
        val samRecord4_1 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, true, "10M")
        val samRecord4_2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")
        val readSplitRecords = listOf(Pair(samRecord4_1, samRecord4_2))
        val pairReadSplit = extractEdgeReads.classifyPairedAlignments("sample1", 4, readSplitRecords,
            refRangeMap, gameteMap, nonconsecSplit)
        assertEquals(AlignmentClass.PAIRREADSPLIT, pairReadSplit, "fail PAIRREADSPLIT")

        // PAIRREADSPLITCONSEC
        val pairReadSplitConsec = extractEdgeReads.classifyPairedAlignments("sample1", 4, readSplitRecords,
            refRangeMap, gameteMap, consecSplit)
        assertEquals(AlignmentClass.PAIRREADSPLITCONSEC, pairReadSplitConsec, "fail PAIRREADSPLITCONSEC")

    }

    @Test
    fun testClassifyAlignments() {
        val extractEdgeReads = ExtractEdgeReads()
        val samRecordFactory = DefaultSAMRecordFactory()
        val samHeader = SAMFileHeader()

        val refRange1 = ReferenceRange("1", 1, 10)
        val refRange2 = ReferenceRange("2", 1, 10)
        val samRecord1 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, true, "10M")
        val samRecord2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, true, "10M")
        val hapIdToSampleGameteOffASM = mapOf("hap1" to listOf(SampleGamete("sample1", 0)),
            "hap2" to listOf(SampleGamete("sample2", 0)), "hap3" to listOf(SampleGamete("sample3", 0)))
        val recordsOffASM = listOf(Pair(samRecord1, null), Pair(samRecord2, null))
        val refRangeMapOffASM = mapOf("hap1" to listOf(refRange1), "hap2" to listOf(refRange1))
        val indexOffASM = mapOf(refRange1.toString() to 0, refRange2.toString() to 0)
        val single = extractEdgeReads.classifyAlignments("sample3", 2,
            recordsOffASM, refRangeMapOffASM, hapIdToSampleGameteOffASM, indexOffASM)
        assertEquals(AlignmentClass.SINGLEOFFASM, single, "failed single classify alignments")

        val samRecord4_1 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, true, "10M")
        val samRecord4_2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")
        val readSplitRecords = listOf(Pair(samRecord4_1, samRecord4_2))
        val refRangeMap = mapOf("hap1" to listOf(refRange1), "hap2" to listOf(refRange2))
        val gameteMap = mapOf("hap1" to listOf(SampleGamete("sample1", 1)),
            "hap2" to listOf(SampleGamete("sample2", 2)))
        val nonconsecSplit = mapOf(refRange1.toString() to 0, refRange2.toString() to 2) // non consecutive ref ranges
        val paired = extractEdgeReads.classifyAlignments("sample1", 4, readSplitRecords,
            refRangeMap, gameteMap, nonconsecSplit)
        assertEquals(AlignmentClass.PAIRREADSPLIT, paired, "failed paired classify alignments")
    }

    @Test
    fun testCliktParams() {
        val extractEdgeReads = ExtractEdgeReads()
        val missingAll = extractEdgeReads.test("")
        assertEquals(missingAll.statusCode, 1)
        assertEquals("Usage: extract-edge-reads [<options>]\n" +
                "\n" +
                "Error: missing option --bam-dir\n" +
                "Error: missing option --hvcf-dir\n" +
                "Error: missing option --sample-name\n" +
                "Error: missing option --output-file-dir\n", missingAll.output)
    }

    @Test
    fun testProcessReads() {
        val extractEdgeReads = ExtractEdgeReads()
        val samRecordFactory = DefaultSAMRecordFactory()
        val samHeader = SAMFileHeader()

        // unaligned classes
        val emptyReads = extractEdgeReads.processReads("empty", 0, emptyList<SAMRecord>(),
            emptyMap<String, List<ReferenceRange>>(), emptyMap<String, List<SampleGamete>>(), emptyMap<String, Int>())
        assertEquals((Pair(AlignmentClass.UNALIGN, listOf(Pair(null,null)))), emptyReads)

        // both reads are either null or are unaligned
        val samRecord1_1 = createSAMRecord(samRecordFactory, samHeader, "read1",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", firstInPair = true, startPos = 1)
        val samRecord1_2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", firstInPair = false, startPos = 30)
        samRecord1_1.readPairedFlag = true
        samRecord1_2.readPairedFlag = true
        samRecord1_1.readUnmappedFlag = true
        samRecord1_1.mateUnmappedFlag = true
        val unalignedRecords = listOf(samRecord1_1, samRecord1_2)
        val refRange1 = ReferenceRange("1", 1, 10)
        val refRange2 = ReferenceRange("2", 1, 10)
        val hapIdToRefRangeMap = mapOf("hap1" to listOf(refRange1), "hap2" to listOf(refRange2))
        val hapIdToSampleGamete = mapOf("hap1" to listOf(SampleGamete("sample1", 1)),
            "hap2" to listOf(SampleGamete("sample2", 2)))
        val refRangeToIndexMap = mapOf(refRange1.toString() to 0, refRange2.toString() to 2) // non consecutive ref ranges
        val unmapped = extractEdgeReads.processReads("unmapped", 2, unalignedRecords, hapIdToRefRangeMap, hapIdToSampleGamete, refRangeToIndexMap)
        assertEquals(
            (Pair(AlignmentClass.UNALIGN, listOf(Pair(unalignedRecords[0],unalignedRecords[1])))),
            unmapped)

        // classify reads
        val samRecord4_1 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap2", "IIIIIIIIII", 1, true, "10M")
        val samRecord4_2 = createSAMRecord(samRecordFactory, samHeader, "read2",  "AAAAAAAAAA",60, "hap1", "IIIIIIIIII", 1, false, "10M")
        samRecord4_1.readPairedFlag = true
        samRecord4_2.readPairedFlag = true
        val recordsForRead = listOf(samRecord4_1, samRecord4_2)
        val refRangeMap = mapOf("hap1" to listOf(refRange1), "hap2" to listOf(refRange2))
        val gameteMap = mapOf("hap1" to listOf(SampleGamete("sample1", 1)),
            "hap2" to listOf(SampleGamete("sample1", 1)))
        val consecSplit = mapOf(refRange1.toString() to 0, refRange2.toString() to 1) // consecutive ref ranges
        val processedReads = extractEdgeReads.processReads("sample1", 2, recordsForRead,
            refRangeMap, gameteMap, consecSplit)
        println(processedReads)
        val recordsGroupedByContig = listOf(Pair(samRecord4_1, samRecord4_2))
        assertEquals(extractEdgeReads.pairOffAlignments(recordsForRead), Pair(samRecord4_1, samRecord4_2))
        assertEquals(Pair(AlignmentClass.PAIRREADSPLITCONSEC, recordsGroupedByContig), processedReads)

    }

}