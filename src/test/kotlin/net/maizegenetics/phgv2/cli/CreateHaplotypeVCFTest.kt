package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.brapi.resetDirs
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.assertEquals
import net.maizegenetics.phgv2.utils.Position
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.extension.ExtendWith

@ExtendWith(TestExtension::class)
class CreateHaplotypeVCFTest {
    companion object {
        private val haplotypeHvcfDir = "${TestExtension.tempDir}haplotype-vcfs/"

        @BeforeAll
        @JvmStatic
        fun setup() {
            // delete, reset the directories
            resetDirs()
            File(haplotypeHvcfDir).mkdirs()
        }

        @AfterAll
        @JvmStatic
        fun teardown() {
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testFullCreateHaplotypeVCF() {

        // The pathVCF file reflects the h.vcf created from running imputation.
        // From this file, in normal processing, a composite fasta would be created
        // and that fasta would be used to align reads, create bam files, and run deep variant.
        // The deepVariantVCF is simulated output from a deepVariant run based on a composite
        // fasta created from the pathVCF file.  The  deepVariantVCF file is used
        // to create a hapltoype vcf.  The pathVCF and deepVariantVCF were manually created for testing.
        // The objective of the test below is to verify the haplotypeVCF created from CreateHaplotypeVCF
        // is correct in terms of where the SNPs fall within each haplotype.
        val pathVCF = "data/test/resequenceHaplotypeVCF/TestSample.h.vcf"
        val deepVariantVCF = "data/test/resequenceHaplotypeVCF/TestSampleDV.vcf"
        val haplotypeVCF = "$haplotypeHvcfDir/sampleDV_RESQ.vcf"
        val sampleName = "TestSampleReseqHaplotype"

        val result = CreateHaplotypeVCF().test("--path-hvcf $pathVCF --variant-vcf $deepVariantVCF --output-file $haplotypeVCF --sample-name $sampleName")
        assertEquals(0, result.statusCode )

        // Verify:  There should be as many lines in the haplotype resequenced VCF as there are in the
        // dv vcf.  Ie, that many variants.

        val dvVCFReader = VCFFileReader(File(deepVariantVCF), false)
        val hapVCFReader = VCFFileReader(File(haplotypeVCF), false)
        val dvVCFRecords = dvVCFReader.iterator().asSequence().toList()
        val hapVCFRecords = hapVCFReader.iterator().asSequence().toList()
        assertEquals(dvVCFRecords.size, hapVCFRecords.size)

        val hapSNPOffsets = mutableListOf<Position>()  // used to verify SNP positions below
        // Verify the Ref values in the output VCF match the Ref values in the input VCF
        for (idx in dvVCFRecords.indices) {
            // contigs will not match as the output VCF will have haplotype id as the contig
            // The start positions will also not match as dvVCf is composite fasta based and hapVCF
            // is haplotype sequence based.
            // But the values in the REF and ALT fields should match.
            assertEquals(dvVCFRecords[idx].reference.displayString, hapVCFRecords[idx].reference.displayString)
            assertEquals(dvVCFRecords[idx].alternateAlleles[0].displayString, hapVCFRecords[idx].alternateAlleles[0].displayString)
            // Get the haplotype id from the contig field of the hapVCF record
            val hapId = hapVCFRecords[idx].contig
            // get the start positions from the hapVCFRecord POS field
            val start = hapVCFRecords[idx].start
            hapSNPOffsets.add(Position(hapId, start))
        }

        // Verify the 5 variants on chrom 1, and the 2 variants on chrom 2 have correct haplotype-relative
        // positions.
        // Looking at the pathVCF haplotype regions, the first 3 variant positions on chrom 1
        // are in the first haplotype (12f0cec9102e84a161866e37072443b7), the 4th one falls in the 5th haplotype
        // (f50fe6d6b3d9a9d305889db977969916) and the 5th one falls in the 9th haplotype (cfa288e041fc4dbe299d20ae1920d258)
        // For chrom 2, the first variant falls in the 1st haplotype (13417ecbb38b9a159e3ca8c9dade7088) and the
        // second variant falls in the 13th variant (90248144d7173c1f1481008c43e65129)
        // The offsets for their start positions are:
        //  chr 1:  12f0cec9102e84a161866e37072443b7 = begins at position 1
        //          f50fe6d6b3d9a9d305889db977969916 = begins at 11001
        //          cfa288e041fc4dbe299d20ae1920d258 = begins at 22001
        // chr 2:  13417ecbb38b9a159e3ca8c9dade7088 = begins at position 1
        //          90248144d7173c1f1481008c43e65129 = begins at 33001
        // The SNP offsets in these haplotypes are calculated as the difference between the start position of the haplotype
        // and the start position of the variant (plus 1).  The SNP offsets are:
        val truthSNPOffsets = mutableListOf<Position>()
        truthSNPOffsets.add(Position("12f0cec9102e84a161866e37072443b7", 17))
        truthSNPOffsets.add(Position("12f0cec9102e84a161866e37072443b7", 63))
        truthSNPOffsets.add(Position("12f0cec9102e84a161866e37072443b7", 158))
        truthSNPOffsets.add(Position("f50fe6d6b3d9a9d305889db977969916", 5))
        truthSNPOffsets.add(Position("cfa288e041fc4dbe299d20ae1920d258", 527))
        truthSNPOffsets.add(Position("13417ecbb38b9a159e3ca8c9dade7088", 98))
        truthSNPOffsets.add(Position("90248144d7173c1f1481008c43e65129", 30))

        // Verify the SNP offsets in the haplotype VCF created by CreateHaplotypeVCF() match the expected values
        for (idx in hapVCFRecords.indices) {
            val hapId = hapVCFRecords[idx].contig
            val start = hapVCFRecords[idx].start
            val pos = Position(hapId, start)
            assertEquals(truthSNPOffsets[idx], pos)
        }
    }

    @Test
    fun testReseqVCFHeader() {
        // Get vcf from the data folder
        val inputVcf = "data/test/smallseq/LineA.g.vcf"
        val sampleNames = listOf("LineA")
        val reader = VCFFileReader(File(inputVcf), false)
        val hapidToLengthMap = mutableMapOf<String,Int>()
        hapidToLengthMap["hap1"] = 1010
        hapidToLengthMap["hap2"] = 1045
        val truthHeaders = reader.fileHeader

        val createResequenceVcf = CreateHaplotypeVCF()
        val vcfHeader = createResequenceVcf.reseqVCFHeader(reader, hapidToLengthMap, sampleNames)

        // Get INFO lines from vcfHeader
        val infoLines = vcfHeader.infoHeaderLines.toList()
        // verify the vcfHeader and truthHeaders have the same INFO lines
        assertEquals(truthHeaders.infoHeaderLines.size, infoLines.size)
        for (i in 0 until truthHeaders.infoHeaderLines.size) {
            assertEquals(truthHeaders.infoHeaderLines.toList()[i], infoLines[i])
        }

        // verify the vcfHeader and truthHeaders have the same FORMAT lines
        val formatLines = vcfHeader.formatHeaderLines.toList()
        assertEquals(truthHeaders.formatHeaderLines.size, formatLines.size)
        for (i in 0 until truthHeaders.formatHeaderLines.size) {
            assertEquals(truthHeaders.formatHeaderLines.toList()[i], formatLines[i])
        }

        // verify the vcfHeader contains 2 contig lines and these contig lines are
        // not the same as the contig lines in the truthHeaders.  The contig lines
        // created in createResequenceVcf.reseqVCFHeader come from the hapidToLengthMap
        // that was passed as a parmater.
        val contigLines = vcfHeader.contigLines.toList()
        assertEquals(2, contigLines.size)
        for (i in 0 until contigLines.size) {
            assert(!truthHeaders.contigLines.contains(contigLines[i]))
        }
        // verify the created VCF header contigLines contains the correct values
        val contigHap1 = contigLines.any { line -> line.toString() == "contig=<ID=hap1,length=1010>" }
        val contigHap2 = contigLines.any { line -> line.toString() == "contig=<ID=hap2,length=1045>" }
        assert(contigHap1)
        assert(contigHap2)

    }

    @Test
    fun testCreateHapPositionMap() {

        val chromToHapData = mutableMapOf<String, MutableList<HaplotypeData>>()
        val chrom1DataList = mutableListOf<HaplotypeData>()

        var hapData = HaplotypeData("hap1", "chr1",1, 1000)
        chrom1DataList.add(hapData)

        hapData = HaplotypeData("hap2", "chr1", 1001, 3500)
        chrom1DataList.add(hapData)

        hapData = HaplotypeData("hap3", "chr1", 4501, 1000)
        chrom1DataList.add(hapData)
        chromToHapData["chr1"] = chrom1DataList

        // Setup the second chromosome
        val chrom2DataList = mutableListOf<HaplotypeData>()
        hapData = HaplotypeData("hap4", "chr2",1, 1000)
        chrom2DataList.add(hapData)

        hapData = HaplotypeData("hap5", "chr2", 1001, 3500)
        chrom2DataList.add(hapData)

        hapData = HaplotypeData("hap6", "chr2",4501, 1000)
        chrom2DataList.add(hapData)
        chromToHapData["chr2"] = chrom2DataList

        val createReseqVCF = CreateHaplotypeVCF()
        val hapToLengthMap = createReseqVCF.createHapPositionMap(chromToHapData)

        assertEquals(6, hapToLengthMap.asMapOfRanges().size)
        assertEquals("hap1", hapToLengthMap.get(Position("chr1", 1)))
        assertEquals("hap1", hapToLengthMap.get(Position("chr1", 1000)))
        assertEquals("hap2", hapToLengthMap.get(Position("chr1", 2500)))
        assertEquals("hap3", hapToLengthMap.get(Position("chr1", 4501)))
        assertEquals("hap4", hapToLengthMap.get(Position("chr2", 1)))
        assertEquals("hap4", hapToLengthMap.get(Position("chr2", 1000)))
        assertEquals("hap5", hapToLengthMap.get(Position("chr2", 2303)))
        assertEquals("hap6", hapToLengthMap.get(Position("chr2", 4501)))

        println("Done with testCreateHapPositionMap")

    }

}