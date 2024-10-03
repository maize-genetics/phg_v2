package net.maizegenetics.phgv2.pathing

import htsjdk.samtools.*
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.Assertions.assertEquals
import java.io.BufferedReader
import java.io.File
import java.io.FileReader

class ExtractEdgeReadsTest {

    private fun createSAMRecord(factory: SAMRecordFactory, header: SAMFileHeader, readName : String,
                                readString : String, mappingQuality : Int, contig : String, BS : String): SAMRecord {
        val samRecord = factory.createSAMRecord(header)
        samRecord.readName = readName
        samRecord.readString = readString
        samRecord.mappingQuality = mappingQuality
        samRecord.referenceName = contig
        samRecord.baseQualityString = BS
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
        BufferedReader(FileReader("table.txt")).use { reader ->
            val header : String = reader.readLine()
            assertEquals(header, "readID\tHapID hits")
            val line1 : String = reader.readLine()
            assertEquals(line1, "readID1\tchr1, chr2, chr3, chr4")
        }
        BufferedReader(FileReader("1fastq.fq")).use { reader ->
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
        BufferedReader(FileReader("2fastq.fq")).use { reader ->

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
}