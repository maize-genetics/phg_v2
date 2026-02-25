package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContextBuilder
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.Position
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.fail

class Hvcf2VcfTest {

    @Test
    fun testCliktParams() {
        val hvcf2Vcf = Hvcf2Vcf()

        //check parameters
        //val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        //        .default("")
        //
        //
        //    val hvcfDir by option(help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        //        .required()
        //
        //    val donorVcfFile by option(help = "Path to the VCF file containing all the PHG SNPs.  This is typically created by running merge-gvcf on the Assembly Gvcf files.")
        //        .required()
        //
        //    val outputFile by option(help = "Output file.")
        //        .required()
        //
        //    //This is needed to create the refSeqDictionary
        //    val referenceFile by option(help = "Path to local Reference FASTA file needed for sequence dictionary")
        //        .required()

        //"--hvcf-dir /path/to/hvcf --donor-vcf-file /path/to/donor.vcf --output-file /path/to/output.vcf --reference-file /path/to/reference.fasta"

        val noHvcfDir = hvcf2Vcf.test("--donor-vcf-file /path/to/donor.vcf --output-file /path/to/output.vcf --reference-file /path/to/reference.fasta")
        assertEquals(noHvcfDir.statusCode, 1)
        assertEquals(
            "Usage: hvcf2vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --hvcf-dir\n", noHvcfDir.output
        )

        val noDonorVcf = hvcf2Vcf.test("--hvcf-dir /path/to/hvcf --output-file /path/to/output.vcf --reference-file /path/to/reference.fasta")
        assertEquals(noDonorVcf.statusCode, 1)
        assertEquals(
            "Usage: hvcf2vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --donor-vcf-file\n", noDonorVcf.output
        )

        val noOutputFile = hvcf2Vcf.test("--hvcf-dir /path/to/hvcf --donor-vcf-file /path/to/donor.vcf --reference-file /path/to/reference.fasta")
        assertEquals(noOutputFile.statusCode, 1)
        assertEquals(
            "Usage: hvcf2vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --output-file\n", noOutputFile.output
        )

        val noRefFile = hvcf2Vcf.test("--hvcf-dir /path/to/hvcf --donor-vcf-file /path/to/donor.vcf --output-file /path/to/output.vcf")
        assertEquals(noRefFile.statusCode, 1)
        assertEquals(
            "Usage: hvcf2vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --reference-file\n", noRefFile.output
        )
    }

    @Test
    fun processHVCFAndBuildVCFTest() {
        fail("Not yet implemented")
    }

    @Test
    fun createRangeHapMapToSampleGameteTest() {
        fail("Not yet implemented")
    }

    @Test
    fun processSingleHvcfFileTest() {
        fail("Not yet implemented")
    }

    @Test
    fun processSingleHvcfVariantTest() {
        //processSingleHvcfVariant(context: VariantContext): List<HvcfRangeHapIdSampleGamete>

        val hvcf2Vcf = Hvcf2Vcf()
        //Make a single sample hvcf variant context haploid and make sure the output
        val singleHaploidVCF = VariantContextBuilder()
            .chr("chr1")
            .start(100L)
            .stop(150L)
            .alleles(listOf(Allele.REF_A, Allele.create("<HAP1>", false)))
            .genotypes(GenotypeBuilder("Sample1",listOf(Allele.create("<HAP1>", false))).make())
            .make()

        val singleHapVcfResult = hvcf2Vcf.processSingleHvcfVariant(singleHaploidVCF)
        //HvcfRangeHapIdSampleGamete(val refRange: ReferenceRange, val hapId: String, val sampleGametes: List<SampleGamete>)
        assertEquals(1, singleHapVcfResult.size)
        assertEquals("HAP1", singleHapVcfResult[0].hapId)
        assertEquals(ReferenceRange("chr1", 100, 150), singleHapVcfResult[0].refRange)
        assertEquals(1, singleHapVcfResult[0].sampleGametes.size)
        assertEquals("Sample1", singleHapVcfResult[0].sampleGametes[0].name)
        assertEquals(0, singleHapVcfResult[0].sampleGametes[0].gameteId)

        //Make a single sample hvcf VariantContext diploid and make sure of the output
        val singleDiploidVCF = VariantContextBuilder()
            .chr("chr1")
            .start(100L)
            .stop(150L)
            .alleles(listOf(Allele.REF_A, Allele.create("<HAP1>", false),Allele.create("<HAP2>", false)))
            .genotypes(GenotypeBuilder("Sample1",listOf(Allele.create("<HAP1>", false),Allele.create("<HAP2>",false))).make())
            .make()

        val singleDiploidVcfResult = hvcf2Vcf.processSingleHvcfVariant(singleDiploidVCF)
        //HvcfRangeHapIdSampleGamete(val refRange: ReferenceRange, val hapId: String, val sampleGametes: List<SampleGamete>)
        assertEquals(2, singleDiploidVcfResult.size)
        assertEquals("HAP1", singleDiploidVcfResult[0].hapId)
        assertEquals(ReferenceRange("chr1", 100, 150), singleDiploidVcfResult[0].refRange)
        assertEquals(1, singleDiploidVcfResult[0].sampleGametes.size)
        assertEquals("Sample1", singleDiploidVcfResult[0].sampleGametes[0].name)
        assertEquals(0, singleDiploidVcfResult[0].sampleGametes[0].gameteId)
        assertEquals("HAP2", singleDiploidVcfResult[1].hapId)
        assertEquals(ReferenceRange("chr1", 100, 150), singleDiploidVcfResult[1].refRange)
        assertEquals(1, singleDiploidVcfResult[1].sampleGametes.size)
        assertEquals("Sample1", singleDiploidVcfResult[1].sampleGametes[0].name)
        assertEquals(1, singleDiploidVcfResult[1].sampleGametes[0].gameteId)


        //Make a multisample hvcf VariantContext haploid
        val multiHaploidVCF = VariantContextBuilder()
            .chr("chr1")
            .start(100L)
            .stop(150L)
            .alleles(listOf(Allele.REF_A, Allele.create("<HAP1>", false),Allele.create("<HAP2>", false)))
            .genotypes(GenotypeBuilder("Sample1",listOf(Allele.create("<HAP1>", false))).make(),
                GenotypeBuilder("Sample2",listOf(Allele.create("<HAP2>", false))).make(),
                GenotypeBuilder("Sample3",listOf(Allele.create("<HAP1>", false))).make())
            .make()


        val multiHaploidVcfResult = hvcf2Vcf.processSingleHvcfVariant(multiHaploidVCF).sortedBy { it.sampleGametes.first().name } //Doing a sort just in case
        //HvcfRangeHapIdSampleGamete(val refRange: ReferenceRange, val hapId: String, val sampleGametes: List<SampleGamete>)
        assertEquals(2, multiHaploidVcfResult.size)
        assertEquals("HAP1", multiHaploidVcfResult[0].hapId)
        assertEquals(ReferenceRange("chr1", 100, 150), multiHaploidVcfResult[0].refRange)
        assertEquals(2, multiHaploidVcfResult[0].sampleGametes.size)
        assertEquals("Sample1", multiHaploidVcfResult[0].sampleGametes[0].name)
        assertEquals(0, multiHaploidVcfResult[0].sampleGametes[0].gameteId)
        assertEquals("Sample3", multiHaploidVcfResult[0].sampleGametes[1].name)
        assertEquals(0, multiHaploidVcfResult[0].sampleGametes[1].gameteId)
        assertEquals("HAP2", multiHaploidVcfResult[1].hapId)
        assertEquals(ReferenceRange("chr1", 100, 150), multiHaploidVcfResult[1].refRange)
        assertEquals(1, multiHaploidVcfResult[1].sampleGametes.size)
        assertEquals("Sample2", multiHaploidVcfResult[1].sampleGametes[0].name)
        assertEquals(0, multiHaploidVcfResult[1].sampleGametes[0].gameteId)

        //Make a multisample hvcf VariantContext diploid
        val multiDiploidVCF = VariantContextBuilder()
            .chr("chr1")
            .start(100L)
            .stop(150L)
            .alleles(listOf(Allele.REF_A, Allele.create("<HAP1>", false),Allele.create("<HAP2>", false)))
            .genotypes(GenotypeBuilder("Sample1",listOf(Allele.create("<HAP1>", false),Allele.create("<HAP1>", false))).make(),
                GenotypeBuilder("Sample2",listOf(Allele.create("<HAP2>", false), Allele.create("<HAP1>", false))).make(),
                GenotypeBuilder("Sample3",listOf(Allele.create("<HAP1>", false))).make())
            .make()

        val multiDiploidVcfResult = hvcf2Vcf.processSingleHvcfVariant(multiDiploidVCF).sortedBy { it.sampleGametes.first().name } //Doing a sort just in case
                //HvcfRangeHapIdSampleGamete(val refRange: ReferenceRange, val hapId: String, val sampleGametes: List<SampleGamete>)
                assertEquals(2, multiDiploidVcfResult.size)
                assertEquals("HAP1", multiDiploidVcfResult[0].hapId)
                assertEquals(ReferenceRange("chr1", 100, 150), multiDiploidVcfResult[0].refRange)
                assertEquals(4, multiDiploidVcfResult[0].sampleGametes.size)
                assertEquals("Sample1", multiDiploidVcfResult[0].sampleGametes[0].name)
                assertEquals(0, multiDiploidVcfResult[0].sampleGametes[0].gameteId)
                assertEquals("Sample1", multiDiploidVcfResult[0].sampleGametes[1].name)
                assertEquals(1, multiDiploidVcfResult[0].sampleGametes[1].gameteId)
                assertEquals("Sample2", multiDiploidVcfResult[0].sampleGametes[2].name)
                assertEquals(1, multiDiploidVcfResult[0].sampleGametes[2].gameteId)
                assertEquals("Sample3", multiDiploidVcfResult[0].sampleGametes[3].name)
                assertEquals(0, multiDiploidVcfResult[0].sampleGametes[3].gameteId)

                assertEquals("HAP2", multiDiploidVcfResult[1].hapId)
                assertEquals(ReferenceRange("chr1", 100, 150), multiDiploidVcfResult[1].refRange)
                assertEquals(1, multiDiploidVcfResult[1].sampleGametes.size)
                assertEquals("Sample2", multiDiploidVcfResult[1].sampleGametes[0].name)
                assertEquals(0, multiDiploidVcfResult[1].sampleGametes[0].gameteId)
    }

    @Test
    fun buildPositionToRefRangeMapTest() {
        //buildPositionToRefRangeMap(ranges: Set<ReferenceRange>): TreeMap<Position, ReferenceRange> {

        val hvcf2Vcf = Hvcf2Vcf()

        //Build a simple reference range set
        val ranges = setOf(ReferenceRange("chr1", 1, 100), ReferenceRange("chr1", 101, 200), ReferenceRange("chr2", 1, 100))

        val treeMap = hvcf2Vcf.buildPositionToRefRangeMap(ranges)
        assertEquals(6, treeMap.size)
        assertEquals(ReferenceRange("chr1", 1, 100), treeMap[Position("chr1", 1)])
        assertEquals(null, treeMap[Position("chr1", 50)])
        assertEquals(ReferenceRange("chr1", 1, 100), treeMap.floorEntry(Position("chr1", 50)).value)
        assertEquals(ReferenceRange("chr1", 101, 200), treeMap[Position("chr1", 101)])
        assertEquals(null, treeMap[Position("chr1", 150)])
        assertEquals(ReferenceRange("chr1", 101, 200), treeMap.floorEntry(Position("chr1", 150)).value)
        assertEquals(ReferenceRange("chr2", 1, 100), treeMap[Position("chr2", 1)])
        assertEquals(null, treeMap[Position("chr2", 50)])
        assertEquals(ReferenceRange("chr2", 1, 100), treeMap.floorEntry(Position("chr2", 50)).value)

    }

    @Test
    fun extractVcfAndExportTest() {
        fail("Not yet implemented")
    }

    @Test
    fun extractVariantsAndWriteTest() {
        fail("Not yet implemented")
    }

}