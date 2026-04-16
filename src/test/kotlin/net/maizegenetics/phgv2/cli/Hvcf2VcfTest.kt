package net.maizegenetics.phgv2.cli

import biokotlin.genome.MAFToGVCF
import biokotlin.seqIO.NucSeqIO
import biokotlin.util.MergeGVCFUtils
import com.github.ajalt.clikt.testing.test
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.utils.HvcfVariant
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.exportVariantInfo
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertNotNull
import org.junit.jupiter.api.assertThrows
import java.io.File
import kotlin.collections.component1
import kotlin.collections.component2
import kotlin.collections.iterator
import kotlin.test.assertEquals
import kotlin.test.fail

class Hvcf2VcfTest {

    @Test
    fun testCliktParams() {
        val hvcf2Vcf = Hvcf2Vcf()

        val noHvcfDir = hvcf2Vcf.test("--pangenome-vcf-file /path/to/pangenome.vcf --output-file /path/to/output.vcf --reference-file /path/to/reference.fasta")
        assertEquals(noHvcfDir.statusCode, 1)
        assertEquals(
            "Usage: hvcf2vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --hvcf-dir\n", noHvcfDir.output
        )

        val noPangenomeVcf = hvcf2Vcf.test("--hvcf-dir /path/to/hvcf --output-file /path/to/output.vcf --reference-file /path/to/reference.fasta")
        assertEquals(noPangenomeVcf.statusCode, 1)
        assertEquals(
            "Usage: hvcf2vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --pangenome-vcf-file\n", noPangenomeVcf.output
        )

        val noOutputFile = hvcf2Vcf.test("--hvcf-dir /path/to/hvcf --pangenome-vcf-file /path/to/pangenome.vcf --reference-file /path/to/reference.fasta")
        assertEquals(noOutputFile.statusCode, 1)
        assertEquals(
            "Usage: hvcf2vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --output-file\n", noOutputFile.output
        )

        val noRefFile = hvcf2Vcf.test("--hvcf-dir /path/to/hvcf --pangenome-vcf-file /path/to/pangenome.vcf --output-file /path/to/output.vcf")
        assertEquals(noRefFile.statusCode, 1)
        assertEquals(
            "Usage: hvcf2vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --reference-file\n", noRefFile.output
        )
    }

    @Test
    fun processHVCFAndBuildVCFTest() {
        //processHVCFAndBuildVCF(dbPath: String, hvcfDir: String, pangenomeVCFFile: String, outputFile: String, referenceFile: String)
        val hvcf2Vcf = Hvcf2Vcf()
        val dbPath = "data/test/hvcf2vcf/asmHvcfs/"
        val hvcfDir = "data/test/hvcf2vcf/imputeHvcfs/"
        val pangenomeVcfFile = "data/test/hvcf2vcf/asmPangenomeVcf/merged.vcf"
        val outputFile = "data/test/hvcf2vcf/temp/imputedVariantsAB_hap_dip.vcf"
        val referenceFile = "data/test/smallseq/Ref.fa"

        hvcf2Vcf.processHVCFAndBuildVCF(dbPath, hvcfDir, pangenomeVcfFile, outputFile, referenceFile)

        //Open up the mergedVCF and make a map of Pos + Line -> Allele call
        val mergedGVCFMap = VCFFileReader(File(pangenomeVcfFile), false).flatMap {
            val position = Position(it.contig, it.start)
            it.genotypesOrderedByName.map { genotype ->
                Pair(Pair(position, genotype.sampleName), genotype.alleles.map { allele -> allele.baseString })
            }
        }.toMap()

        //For LineAB haploid
        //From
        // chr1 1 - 6500 LineA
        // chr1 6501 - 16,500 LineB
        // chr2 1 - 6500 LineB
        // chr2 6501 - 16,500 LineA
        //For diploid
        // chr1 1 - 5500 LineA|LineA
        // chr1 5501 - 11,000 LineA|LineB
        // chr1 11,001 - 16,500 LineB|LineB
        // chr2 1 - 5500 LineB|LineB
        // chr2 5501 - 11,000 LineB|LineA
        // chr2 11,001 - 16,500 LineA|LineA
        //Walk through the imputed vcf checking the variants
        VCFFileReader(File(outputFile),false).forEach {
            val position = Position(it.contig, it.start)
            val lineAsVariant = mergedGVCFMap[Pair(position,"LineA")]!!.first()
            val lineBsVariant = mergedGVCFMap[Pair(position,"LineB")]!!.first()

            if(position in Position("1",1) .. Position("1",5500)) {
                //Hap LineA, Diploid LineA|LineA
                assertEquals(lineAsVariant, it.getGenotype("LineAB_haploid").alleles.map { allele ->allele.baseString }.first())
                assertEquals(listOf(lineAsVariant,lineAsVariant), it.getGenotype("LineAB_diploid").alleles.map { allele ->allele.baseString })
            }
            else if(position in Position("1",5501) .. Position("1",6500)) {
                //Hap LineA, diploid LineA | LineB
                assertEquals(lineAsVariant, it.getGenotype("LineAB_haploid").alleles.map { allele ->allele.baseString }.first())
                assertEquals(listOf(lineAsVariant,lineBsVariant), it.getGenotype("LineAB_diploid").alleles.map { allele ->allele.baseString })
            }
            else if(position in Position("1",6500) .. Position("1",11000)) {
                //Hap LineB, diploid LineA | LineB
                assertEquals(lineBsVariant, it.getGenotype("LineAB_haploid").alleles.map { allele ->allele.baseString }.first())
                assertEquals(listOf(lineAsVariant,lineBsVariant), it.getGenotype("LineAB_diploid").alleles.map { allele ->allele.baseString })
            }
            else if(position in Position("1",11001) .. Position("1",16500)) {
                //Hap LineB, diploid LineB | LineB
                assertEquals(lineBsVariant, it.getGenotype("LineAB_haploid").alleles.map { allele ->allele.baseString }.first())
                assertEquals(listOf(lineBsVariant,lineBsVariant), it.getGenotype("LineAB_diploid").alleles.map { allele ->allele.baseString })
            }
            else if(position in Position("2",1) .. Position("2",5500)) {
                //Hap LineB, diploid LineB|LineB
                assertEquals(lineBsVariant, it.getGenotype("LineAB_haploid").alleles.map { allele ->allele.baseString }.first())
                assertEquals(listOf(lineBsVariant,lineBsVariant), it.getGenotype("LineAB_diploid").alleles.map { allele ->allele.baseString })
            }
            else if(position in Position("2",5501) .. Position("2",6500)) {
                //Hap LineB, diploid LineB | LineA
                assertEquals(lineBsVariant, it.getGenotype("LineAB_haploid").alleles.map { allele ->allele.baseString }.first())
                assertEquals(listOf(lineBsVariant,lineAsVariant), it.getGenotype("LineAB_diploid").alleles.map { allele ->allele.baseString })
            }
            else if(position in Position("2",6501) .. Position("2",11000)) {
                //Hap LineA, diploid LineB | LineA
                assertEquals(lineAsVariant, it.getGenotype("LineAB_haploid").alleles.map { allele ->allele.baseString }.first())
                assertEquals(listOf(lineBsVariant,lineAsVariant), it.getGenotype("LineAB_diploid").alleles.map { allele ->allele.baseString })
            }
            else if(position in Position("2",11001) .. Position("2",16500)) {
                //Hap LineA, diploid LineA | LineA
                assertEquals(lineAsVariant, it.getGenotype("LineAB_haploid").alleles.map { allele ->allele.baseString }.first())
                assertEquals(listOf(lineAsVariant,lineAsVariant), it.getGenotype("LineAB_diploid").alleles.map { allele ->allele.baseString })
            }
            else {
                fail("Position ${position.contig}:${position.position} is out of the expected range")
            }
        }

        File(outputFile).deleteRecursively()
    }

    @Test
    fun createRangeHapMapToSampleGameteTest() {
        val hvcf2Vcf = Hvcf2Vcf()

        val hvcfDir = "data/test/hvcf2vcf/asmHvcfs/hvcf_files/"

        val hvcfRangeValues = hvcf2Vcf.createRangeHapMapToSampleGamete(hvcfDir)

        val truthOutput = buildTruthHVCFRecords().associate { Pair(it.refRange, it.hapId) to it.sampleGametes }

        for((key,rangeValue) in hvcfRangeValues) {
            assertTrue(truthOutput.containsKey(key))
            assertEquals(truthOutput.getValue(key).sorted(), rangeValue.sorted())
        }
    }

    @Test
    fun createRangeHapMapToSampleGameteWithMissingTest() {
        val hvcf2Vcf = Hvcf2Vcf()

        val hvcfDir = "data/test/hvcf2vcf/asmHvcfs/hvcf_files_with_missing/"

        val hvcfRangeValues = hvcf2Vcf.createRangeHapMapToSampleGamete(hvcfDir)

        for((key,rangeValue) in hvcfRangeValues) {
            println("${key}\t${rangeValue}")
        }
    }


    @Test
    fun processSingleHvcfFileTest() {
        val hvcf2Vcf = Hvcf2Vcf()

        //using this input hvcf as it is short
        val inputHvcf = "data/test/hvcf2vcf/asmHvcfs/hvcf_files/LineA.h.vcf"

        val hvcfRangeValues = hvcf2Vcf.processSingleHvcfFile(File(inputHvcf))

        val truthOutput = buildTruthHVCFRecords()

        for(range in hvcfRangeValues) {
            assertTrue(truthOutput.contains(range))
        }
    }

    fun buildTruthHVCFRecords(): Set<HvcfRangeHapIdSampleGamete> {
        return setOf(
            HvcfRangeHapIdSampleGamete(ReferenceRange("1",1,1000),"12f0cec9102e84a161866e37072443b7",  listOf(SampleGamete("LineA", 0), SampleGamete("LineA2", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("1",1001,5500),"3149b3144f93134eb29661bade697fc6",  listOf(SampleGamete("LineA", 0), SampleGamete("LineA2", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("1",5501,6500),"1b568197f6f329ec5b71f66e49a732fb",  listOf(SampleGamete("LineA", 0), SampleGamete("LineA2", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("1",6501,11000),"369464a8743d2e40ad83d1375c196bdd",  listOf(SampleGamete("LineA", 0), SampleGamete("LineA2", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("1",11001,12000),"f50fe6d6b3d9a9d305889db977969916",  listOf(SampleGamete("LineA", 0), SampleGamete("LineA2", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("1",12001,16500),"d4c8b5505d7046b41d7f69b246063ebb",  listOf(SampleGamete("LineA", 0), SampleGamete("LineA2", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("2",1,1000),"13417ecbb38b9a159e3ca8c9dade7088",  listOf(SampleGamete("LineA", 0), SampleGamete("LineA2", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("2",1001,5500),"c16ac825052c0456069f6408652aadf8",  listOf(SampleGamete("LineA", 0), SampleGamete("LineA2", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("2",5501,6500),"50044914d5111c5b5ec58c9d720e3b2d",  listOf(SampleGamete("LineA", 0), SampleGamete("LineA2", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("2",6501,11000),"c472bb8d63e19218a4089821ae666db3",  listOf(SampleGamete("LineA", 0), SampleGamete("LineA2", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("2",11001,12000),"3ec680649615da0685b8c245e0f196e2",  listOf(SampleGamete("LineA", 0), SampleGamete("LineA2", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("1",1,1000), "4fc7b8af32ddd74e07cb49d147ef1938", listOf(SampleGamete("LineB", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("1",1001,5500), "8967fabf10e55d881caa6fe192e7d4ca", listOf(SampleGamete("LineB", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("1",5501,6500), "05efe15d97db33185b64821791b01b0f", listOf(SampleGamete("LineB", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("1",6501,11000), "8f7de1a693aa15fb8fb7b85e7a8b5e95", listOf(SampleGamete("LineB", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("1",11001,12000), "6b5f46bd5c31917af3ab6c3ccc8668cd", listOf(SampleGamete("LineB", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("1",12001,16500), "aff71f19de448514a6d9208b1fcb4e8a", listOf(SampleGamete("LineB", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("2",1,1000), "180417a01edbfed525d7c238910e0ff4", listOf(SampleGamete("LineB", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("2",1001,5500), "8bcf66e8c49da2d9ad8cbafa0bb7a93d", listOf(SampleGamete("LineB", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("2",5501,6500), "45b121547c7ae517a181fdd2621495c4", listOf(SampleGamete("LineB", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("2",6501,11000), "bc94073196b0b2c13e62b5fa47c76b51", listOf(SampleGamete("LineB", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("2",11001,12000), "b787382b1337fd694dd8d77de0141da4", listOf(SampleGamete("LineB", 0))),
            HvcfRangeHapIdSampleGamete(ReferenceRange("2",12001,16500), "6fb2de47c835bd9ab026c02d62f49807", listOf(SampleGamete("LineB", 0)))
        )
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
        //Note do not need to unit test extractVariantsAndWrite as it needs a writer and is a very simple function that is just called by this.
        //It was more to make the code more modular/smaller but likely does not need its own test.

        val hvcf2Vcf = Hvcf2Vcf()
        val refFasta = "data/test/smallseq/Ref.fa"
        val pangenomeVCF = "data/test/hvcf2vcf/asmPangenomeVcf/merged.vcf"
        val outputFile = "data/test/hvcf2vcf/temp/imputed.vcf"

        val refSeq = NucSeqIO(refFasta).readAll()
        //Make 2 reference ranges, 1-100 and 101-200.  The first refRange will have variants from LineA and the 2nd will be from LineB
        //The hapids don't even need to be real just need to match
        val asmHapIdMap: Map<Pair<ReferenceRange, String>, List<HvcfVariant>> = mapOf(
            Pair(ReferenceRange("1", 1, 100), "LineA") to
                listOf(HvcfVariant(ReferenceRange("1",1, 100), "LineA", "HAPA1")),
            Pair(ReferenceRange("1",1,100), "LineB") to listOf(HvcfVariant(ReferenceRange("1",1, 100), "LineB", "HAPB1")),
            Pair(ReferenceRange("1",101,200), "LineA") to listOf(HvcfVariant(ReferenceRange("1",101, 200), "LineA", "HAPA2")),
            Pair(ReferenceRange("1",101,200), "LineB") to listOf(HvcfVariant(ReferenceRange("1",101, 200), "LineB", "HAPB2"))
            )

        val refRangeAndHapIdMap = mapOf<Pair<ReferenceRange,String>, List<SampleGamete>>(
            Pair(ReferenceRange("1",1,100), "HAPA1") to listOf(SampleGamete("OutputAB", 0)),
            Pair(ReferenceRange("1",101,200), "HAPB2") to listOf(SampleGamete("OutputAB", 0)),
            )


        val positionToRangeMap = hvcf2Vcf.buildPositionToRefRangeMap(setOf(ReferenceRange("1", 1, 100), ReferenceRange("1", 101, 200)))

        hvcf2Vcf.extractVcfAndExport(asmHapIdMap,refRangeAndHapIdMap, positionToRangeMap, pangenomeVCF, outputFile, refSeq)

        //check the results
        val expectedFile = "data/test/hvcf2vcf/truth/expectedSimpleFile.vcf"

        compareVCFs(expectedFile, outputFile)
        //flip it just to be sure
        compareVCFs(outputFile, expectedFile)


        //delete outputFile
        //Eventually use temp dir but requires some changes as it needs to pass a file around...
        File(outputFile).deleteRecursively()
    }

    fun compareVCFs(expectedVCF: String, actualVCF: String) {
        //Load in expectedFile and outputFile and compare results
        val truthContextMap = VCFFileReader(File(expectedVCF),false).associateBy { Position(it.contig, it.start) }

        val actualContextMap = VCFFileReader(File(actualVCF),false).associateBy { Position(it.contig, it.start) }

        //Loop through each map and compare the position, alleles and genotype calls
        for((truthContextKey, truthContext) in truthContextMap) {
            assertTrue(actualContextMap.containsKey(truthContextKey))
            val actualContext = actualContextMap[truthContextKey]!!
            assertEquals(truthContext.contig, actualContext.contig)
            assertEquals(truthContext.start, actualContext.start)
            assertEquals(truthContext.alleles, actualContext.alleles)
            val truthGenotypes = truthContext.genotypes
            val actualGenotypes = actualContext.genotypes
            for(sample in truthGenotypes.sampleNamesOrderedByName) {
                val actualGenotype = actualGenotypes.get(sample)
                assertEquals(truthGenotypes.get(sample).alleles, actualGenotype.alleles)

            }
        }
    }

    /**
     * Code to build the LineBGVCF and merged VCF
     */
    //@Test
    fun buildLineBAndSharedVCF() {
        val lineBMaf = "data/test/hvcf2vcf/temp/LineB.maf"
        val outputLineB = "data/test/hvcf2vcf/temp/LineB.g.vcf"
        val lineAMaf = "data/test/hvcf2vcf/temp/LineA.maf"
        val outputLineA = "data/test/hvcf2vcf/temp/LineA.g.vcf"
        val ref = "data/test/smallSeq/Ref.fa"
        val mafToGVCF = MAFToGVCF()

        val refSeq = NucSeqIO(ref).readAll()

        val variantsLineB = mafToGVCF.getVariantContextsfromMAF(lineBMaf, refSeq, "LineB", false, false)
        exportVariantInfo("LineB", variantsLineB["LineB"]!!.sortedBy { variant -> Position(variant.chr, variant.startPos) }, outputLineB, refSeq, setOf(),mafToGVCF)

        val variantsLineA = mafToGVCF.getVariantContextsfromMAF(lineAMaf, refSeq, "LineA", false, false)
        exportVariantInfo("LineA", variantsLineA["LineA"]!!.sortedBy { variant -> Position(variant.chr, variant.startPos) }, outputLineA, refSeq, setOf(),mafToGVCF)

        val mergedVCF = "data/test/hvcf2vcf/temp/merged.vcf"
        //Now we need to merge this pluse LineAs.  I moved LineA.g.vcf to temp so we can merge
        MergeGVCFUtils.mergeGVCFs("data/test/hvcf2vcf/temp/", mergedVCF, null)
    }

    @Test
    fun buildVariantContextTest() {
        val hvcf2Vcf = Hvcf2Vcf()

        //This is an input variant from multiple haploids so it should only have one allele
        val context = VariantContextBuilder(".","1",10L, 10L, listOf(Allele.REF_A, Allele.ALT_G))
            .genotypes(GenotypeBuilder("Sample1", listOf(Allele.ALT_G)).make())
            .make()

        val positionToRefRangeMissing = hvcf2Vcf.buildPositionToRefRangeMap(setOf(ReferenceRange("1", 11, 20)))

        val nullOutputVC = hvcf2Vcf.buildVariantContext(context, positionToRefRangeMissing, emptyMap(), emptyMap())
        assertEquals(null, nullOutputVC)

        //Check chroms not matching
        val positionToRefRangeChrWrong = hvcf2Vcf.buildPositionToRefRangeMap(setOf(ReferenceRange("chr1", 11, 20)))
        val nullOutputChrWrong = hvcf2Vcf.buildVariantContext(context, positionToRefRangeChrWrong, emptyMap(), emptyMap())
        assertEquals(null, nullOutputChrWrong)

        //Test first that the position is not in positionToRangeMap
        val positionToRangeMap = hvcf2Vcf.buildPositionToRefRangeMap(setOf(ReferenceRange("1", 1, 100), ReferenceRange("1", 101, 200)))
        val asmHapIdMap: Map<Pair<ReferenceRange, String>, List<HvcfVariant>> = mapOf(Pair(ReferenceRange("1", 1, 100), "Sample1") to
                listOf(HvcfVariant(ReferenceRange("1",1, 100), "Sample1", "HAP1")))

        val refRangeAndHapIdMap: Map<Pair<ReferenceRange, String>, List<SampleGamete>> = mapOf(Pair(ReferenceRange("1", 1, 100), "HAP1") to listOf(SampleGamete("OutputSample1", 0)))


        val outputVCF = hvcf2Vcf.buildVariantContext(context, positionToRangeMap, asmHapIdMap, refRangeAndHapIdMap)
        assertNotNull(outputVCF)
        assertEquals("OutputSample1", outputVCF.sampleNames.first())
        //check the GT field should have alt G calls
        assertEquals(listOf(Allele.REF_A, Allele.ALT_G), outputVCF.alleles)
        assertEquals("1", outputVCF.contig)
        assertEquals(10, outputVCF.start)
        assertEquals(10, outputVCF.end)

        val firstGT = outputVCF.genotypes.first()
        assertEquals("OutputSample1", firstGT.sampleName)
        assertEquals(Allele.ALT_G, firstGT.alleles[0])

    }

    @Test
    fun buildOutputGenotypesTest() {
        //buildOutputGenotypes(sampleGameteAndAllelePairs: List<Pair<SampleGamete, Allele>>): List<Genotype>
        val hvcf2Vcf = Hvcf2Vcf()

        //Make up both haploid and diploid only and mixed lists
        val haploidList = listOf(Pair(SampleGamete("Sample1", 0), Allele.create("<HAP1>", false)),
            Pair(SampleGamete("Sample2", 0), Allele.create("<HAP2>", false)),
            Pair(SampleGamete("Sample3", 0), Allele.create("<HAP1>", false)))

        val haploidGenotypes = hvcf2Vcf.buildOutputGenotypes(haploidList)
        //check to make sure that we have the right genotypes out
        val truthHaploidGenotypes = mapOf<String, Genotype>(
            "Sample1" to GenotypeBuilder("Sample1", listOf(Allele.create("<HAP1>", false))).make(),
            "Sample2" to GenotypeBuilder("Sample2", listOf(Allele.create("<HAP2>", false))).make(),
            "Sample3" to GenotypeBuilder("Sample3", listOf(Allele.create("<HAP1>", false))).make()
        )

        for(hapGenotype in haploidGenotypes) {
            assertTrue(truthHaploidGenotypes.containsKey(hapGenotype.sampleName))
            val truthGenotype = truthHaploidGenotypes[hapGenotype.sampleName]!!
            //make sure the actual allele calls match
            assertEquals(truthGenotype.alleles, hapGenotype.alleles)
        }


        //Check diploid only
        val diploidList = listOf(Pair(SampleGamete("Sample1", 0), Allele.create("<HAP1>", false)),
            Pair(SampleGamete("Sample1", 1), Allele.create("<HAP1>", false)),
            Pair(SampleGamete("Sample2", 0), Allele.create("<HAP2>", false)),
            Pair(SampleGamete("Sample2", 1), Allele.create("<HAP1>", false)),
            Pair(SampleGamete("Sample3", 0), Allele.create("<HAP1>", false)),
            Pair(SampleGamete("Sample3", 1), Allele.create("<HAP2>", false)))

        val diploidGenotypes = hvcf2Vcf.buildOutputGenotypes(diploidList)

        val truthDiploidGenotypes = mapOf<String, Genotype>(
            "Sample1" to GenotypeBuilder("Sample1", listOf(Allele.create("<HAP1>", false), Allele.create("<HAP1>", false))).make(),
            "Sample2" to GenotypeBuilder("Sample2", listOf(Allele.create("<HAP2>", false), Allele.create("<HAP1>", false))).make(),
            "Sample3" to GenotypeBuilder("Sample3", listOf(Allele.create("<HAP1>", false), Allele.create("<HAP2>", false))).make()
        )

        for(diploidGenotype in diploidGenotypes) {
            assertTrue(truthDiploidGenotypes.contains(diploidGenotype.sampleName))
            val truthGenotype = truthDiploidGenotypes[diploidGenotype.sampleName]!!
            //make sure the actual allele calls match
            assertEquals(truthGenotype.alleles, diploidGenotype.alleles)
        }

        //Do a mixed sample
        val mixedList = listOf(Pair(SampleGamete("Sample1", 0), Allele.create("<HAP1>", false)),
            Pair(SampleGamete("Sample1", 1), Allele.create("<HAP1>", false)),
            Pair(SampleGamete("Sample2", 0), Allele.create("<HAP2>", false)),
            Pair(SampleGamete("Sample3", 0), Allele.create("<HAP1>", false)),
            Pair(SampleGamete("Sample3", 1), Allele.create("<HAP2>", false)))

        val mixedGenotypes = hvcf2Vcf.buildOutputGenotypes(mixedList)

        val truthMixedGenotypes = mapOf<String, Genotype>(
            "Sample1" to GenotypeBuilder("Sample1", listOf(Allele.create("<HAP1>", false), Allele.create("<HAP1>", false))).make(),
            "Sample2" to GenotypeBuilder("Sample2", listOf(Allele.create("<HAP2>", false))).make(),
            "Sample3" to GenotypeBuilder("Sample3", listOf(Allele.create("<HAP1>", false), Allele.create("<HAP2>", false))).make()
        )

        for(mixedGenotype in mixedGenotypes) {
            assertTrue(truthMixedGenotypes.containsKey(mixedGenotype.sampleName))

            val truthGenotype = truthMixedGenotypes[mixedGenotype.sampleName]!!
            assertEquals(truthGenotype.alleles, mixedGenotype.alleles)
        }

    }

    @Test
    fun extractAllelesForEachSampleGameteTest() {
        val hvcf2vcf = Hvcf2Vcf()

        val variantContext = VariantContextBuilder(".","1",10L,10L, listOf(Allele.REF_A, Allele.ALT_G))
            .genotypes(GenotypeBuilder("Sample1", listOf(Allele.ALT_G)).make())
            .make()

        //First check for errors being thrown
        //First make asmHapIdMap not have a valid refRange.  Should return empty list.
        assertEquals(0,hvcf2vcf.extractAllelesForEachSampleGamete(variantContext,
            mapOf(Pair(ReferenceRange("1", 1, 5), "Sample1") to listOf(
                HvcfVariant(
                    ReferenceRange("1", 1, 5),
                    "Sample1",
                    "Hap1"
                )
            )),
            ReferenceRange("1", 10, 10),
            mapOf(Pair(ReferenceRange("1", 1, 5), "HAP1") to listOf(SampleGamete("Sample1", 0)))
        ).size )

        //This needs to have an empty value.  This is effectively saying that there is no haplotype information for this variant's sample so we should not impute anything and just return an empty list
        val emptyList = hvcf2vcf.extractAllelesForEachSampleGamete(variantContext,
            mapOf(Pair(ReferenceRange("1", 1, 5), "Sample1") to listOf(
                HvcfVariant(
                    ReferenceRange("1", 1, 5),
                    "Sample1",
                    "HAP1"
                )
            ),
                Pair(ReferenceRange("1",6,20),"Sample1") to listOf(
                    HvcfVariant(
                        ReferenceRange("1", 6, 20),
                        "Sample1",
                        "HAP1"
                    )
                )

            ),
            ReferenceRange("1", 6, 20),
            mapOf(Pair(ReferenceRange("1", 1, 5), "HAP1") to listOf(SampleGamete("Sample1", 0)))
        )

        assertEquals(0, emptyList.size)

        //Check that it actually works if the maps are consistent with what is requested
        val allelesForSampleGamete = hvcf2vcf.extractAllelesForEachSampleGamete(variantContext,
            mapOf(Pair(ReferenceRange("1", 1, 5), "Sample1") to listOf(
                HvcfVariant(
                    ReferenceRange("1", 1, 5),
                    "Sample1",
                    "HAP1"
                )
            ),
                Pair(ReferenceRange("1",6,20),"Sample1") to listOf(
                    HvcfVariant(
                        ReferenceRange("1", 6, 20),
                        "Sample1",
                        "HAP1"
                    )
                )

            ),
            ReferenceRange("1", 6, 20),
            mapOf(Pair(ReferenceRange("1", 1, 5), "HAP1") to listOf(SampleGamete("Sample1", 0)),
                Pair(ReferenceRange("1", 6, 20), "HAP1") to listOf(SampleGamete("Sample1", 0)))
        )

        //List<Pair<SampleGamete, Allele>>
        assertEquals(1, allelesForSampleGamete.size)
        assertEquals("Sample1", allelesForSampleGamete[0].first.name)
        assertEquals(0, allelesForSampleGamete[0].first.gameteId)
        assertEquals(Allele.ALT_G, allelesForSampleGamete[0].second)
    }

}