package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFHeader
import htsjdk.variant.vcf.VCFHeaderLine
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.createGenericHeader
import org.junit.Ignore
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.junit5.JUnit5Asserter.fail

class ConvertVcf2Ps4gFileTest {

    @Test
    fun testCliktParams() {
        val convertVcf2Ps4gFile = ConvertVcf2Ps4gFile()

        val noSampleVCF = convertVcf2Ps4gFile.test("--gamete-vcf testDir.vcf --output-dir testDir")
        assertEquals(1, noSampleVCF.statusCode)
        assertEquals("Usage: convert-vcf2ps4g-file [<options>]\n\n" +
                "Error: missing option --sample-vcf\n", noSampleVCF.stderr)

        val noGameteVCF = convertVcf2Ps4gFile.test("--sample-vcf testDir.vcf --output-dir testDir")
        assertEquals(1, noGameteVCF.statusCode)
        assertEquals("Usage: convert-vcf2ps4g-file [<options>]\n\n" +
                "Error: missing option --gamete-vcf\n", noGameteVCF.stderr)

        val noOutputDir = convertVcf2Ps4gFile.test("--sample-vcf testDir.vcf --gamete-vcf testDir.vcf")
        assertEquals(1, noOutputDir.statusCode)
        assertEquals("Usage: convert-vcf2ps4g-file [<options>]\n\n" +
                "Error: missing option --output-dir\n", noOutputDir.stderr)
    }


    @Test
    fun testCreatePositionSampleGameteLookup() {
        val inputVCFFile = File("data/test/ps4gTests/refPanel.vcf")
        val convertVcf2Ps4gFile = ConvertVcf2Ps4gFile()
        val positionSampleGameteLookup = convertVcf2Ps4gFile.createPositionSampleGameteLookup(inputVCFFile.path)

        val sample1Gamete = SampleGamete("sample1", 0)
        val sample2Gamete = SampleGamete("sample2", 0)
        val sample3Gamete = SampleGamete("sample3", 0)
        val sample4Gamete = SampleGamete("sample4", 0)


        assertEquals(40, positionSampleGameteLookup.size)
        for(i in 0 until 40) {
            val position = Position("chr1", i * 10 + 1)
            val sampleGameteMap = positionSampleGameteLookup[position]
            assertEquals(2, sampleGameteMap?.size)

            assertEquals(2, sampleGameteMap?.get("A")?.size)
            assertEquals(2, sampleGameteMap?.get("C")?.size)

            val sampleGameteA = sampleGameteMap?.get("A")
            val sampleGameteC = sampleGameteMap?.get("C")

            when {
                i % 2 == 0 -> {
                    // sample1 and sample2 are homozygous for A
                    assertTrue(sampleGameteA?.contains(sample1Gamete) ?: false)
                    assertTrue(sampleGameteA?.contains(sample2Gamete) ?: false)
                    assertTrue(sampleGameteC?.contains(sample3Gamete) ?: false)
                    assertTrue(sampleGameteC?.contains(sample4Gamete) ?: false)
                }

                i % 3 == 0 -> {
                    // sample1 and sample4 are homozygous for A
                    assertTrue(sampleGameteA?.contains(sample1Gamete) ?: false)
                    assertTrue(sampleGameteA?.contains(sample4Gamete) ?: false)
                    assertTrue(sampleGameteC?.contains(sample2Gamete) ?: false)
                    assertTrue(sampleGameteC?.contains(sample3Gamete) ?: false)
                }
                else -> {
                    // sample1 and sample3 are homozygous for A
                    assertTrue(sampleGameteA?.contains(sample1Gamete) ?: false)
                    assertTrue(sampleGameteA?.contains(sample3Gamete) ?: false)
                    assertTrue(sampleGameteC?.contains(sample2Gamete) ?: false)
                    assertTrue(sampleGameteC?.contains(sample4Gamete) ?: false)
                }
            }
        }
    }

    @Test
    fun testCreatePS4GData() {
        fail("Not yet implemented")
    }

    @Test
    fun testProcessVariantPosition() {
        // fun processVariantPosition(
        //        position: Position,
        //        sampleNameToIdxMap: Map<String, Int>,
        //        positionSampleGameteLookup: Map<Position, Map<String, List<SampleGamete>>>,
        //        record: VariantContext,
        //        sampleGameteCount: MutableMap<SampleGamete, MutableMap<SampleGamete, Int>>,
        //        gameteToCountMap: MutableMap<SampleGamete, SampleGameteCountMaps>,
        //        gameteToIdxMap: Map<SampleGamete, Int>
        //    )
        


    }

    @Test
    fun testProcessSingleGenotype() {
        //create the allele map
        val alleleMap = mapOf(
            "A" to listOf(SampleGamete("LineA",0), SampleGamete("LineB",0)),
            "C" to listOf(SampleGamete("LineC",0), SampleGamete("LineD",0)),
            "G" to listOf(SampleGamete("LineD",1), SampleGamete("LineF",0))
        )

        val gameteToIdxMap = mapOf<SampleGamete,Int>(SampleGamete("LineA",0) to 1,
            SampleGamete("LineB",0) to 2,
            SampleGamete("LineC",0) to 3,
            SampleGamete("LineD",0) to 4,
            SampleGamete("LineD",1) to 5,
            SampleGamete("LineF",0) to 6
        )

        val genotype1 = GenotypeBuilder().name("genotype1")
            .alleles(listOf(Allele.REF_A, Allele.ALT_C)).make()

        val encodedPosition1 = PS4GUtils.encodePositionFromIdxAndPos(1,1000)

        val convertVcf2Ps4gFile = ConvertVcf2Ps4gFile()

        val sampleGameteCountMap = mutableMapOf<SampleGamete, MutableMap<SampleGamete,Int>>()
        val gameteToCountMap = mutableMapOf<SampleGamete, SampleGameteCountMaps>()

        assertEquals(0, sampleGameteCountMap.size)
        assertEquals(0, gameteToCountMap.size)

        convertVcf2Ps4gFile.processSingleGenotype(genotype1, 0, Allele.REF_A, alleleMap,
            sampleGameteCountMap, gameteToCountMap, gameteToIdxMap, encodedPosition1)

        assertEquals(1, sampleGameteCountMap.size)
        assertEquals(1, gameteToCountMap.size)
        assertEquals(2, sampleGameteCountMap[SampleGamete("genotype1",0)]?.size)
        val genotype1CountMap = gameteToCountMap[SampleGamete("genotype1",0)]?.countMap!!

        assertEquals(1, genotype1CountMap.get(Pair(encodedPosition1, listOf(1,2))))

        convertVcf2Ps4gFile.processSingleGenotype(genotype1, 1, Allele.ALT_C, alleleMap,
            sampleGameteCountMap, gameteToCountMap, gameteToIdxMap, encodedPosition1)

        assertEquals(2, sampleGameteCountMap.size)
        assertEquals(2, gameteToCountMap.size)
        assertEquals(2, sampleGameteCountMap[SampleGamete("genotype1",1)]?.size)

    }

    @Test
    fun testCreateGameteToIdxMap() {
        val convertVcf2Ps4gFile = ConvertVcf2Ps4gFile()

        val positionSampleGameteLookup = mapOf<Position, Map<String, List<SampleGamete>>>(
            Position("chr1", 1) to mapOf(
                "imputeSample1" to listOf(SampleGamete("sample1", 0)),
                "imputeSample2" to listOf(SampleGamete("sample2", 0))),
            Position("chr1", 2) to mapOf(
                "imputeSample1" to listOf(SampleGamete("sample1", 0), SampleGamete("sample1", 1)),
                "imputeSample2" to listOf(SampleGamete("sample1", 1))),
            Position("chr1", 3) to mapOf(
                "imputeSample1" to listOf(SampleGamete("sample1", 0), SampleGamete("sample1", 1)),
                "imputeSample2" to listOf(SampleGamete("sample4", 0), SampleGamete("sample3", 0))),
        )

        val truthMap = mapOf<SampleGamete, Int>(
            SampleGamete("sample1", 0) to 0,
            SampleGamete("sample1", 1) to 1,
            SampleGamete("sample2", 0) to 2,
            SampleGamete("sample4", 0) to 4,
            SampleGamete("sample3", 0) to 3
        )



        val gameteToIdxMap = convertVcf2Ps4gFile.createGameteToIdxMap(positionSampleGameteLookup)
        assertEquals(truthMap.size, gameteToIdxMap.size)
        truthMap.forEach { (sampleGamete, idx) ->
            assertEquals(idx, gameteToIdxMap[sampleGamete])
        }

    }

    @Ignore
    @Test
    fun buildSimpleReferencePanelVCF() {
        //This test was what was run to get the refPanel.vcf.  I am leaving it in so it is clear what was done to make the file
        val outputFile = "data/test/ps4gTests/refPanel.vcf"
        val header = createGenericHeader(listOf("sample1", "sample2", "sample3", "sample4"),emptySet())
        val writer = VariantContextWriterBuilder()
            .unsetOption(Options.INDEX_ON_THE_FLY)
            .setOutputFile(File(outputFile))
            .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
            .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
            .build()

        writer.writeHeader(header)

        (0 until 400 step 10).forEach { idx ->
            //build a variant context for the 4 samples
            val vcb = VariantContextBuilder("src", "chr1", (idx+1).toLong(), (idx + 1.toLong()), listOf(Allele.REF_A, Allele.ALT_C),)

            val genotypes = if(idx % 20 == 0) {
                listOf<Genotype>(
                    GenotypeBuilder("sample1", listOf(Allele.REF_A)).make(),
                    GenotypeBuilder("sample2", listOf(Allele.REF_A)).make(),
                    GenotypeBuilder("sample3", listOf(Allele.ALT_C)).make(),
                    GenotypeBuilder("sample4", listOf(Allele.ALT_C)).make()
                )
            }
            else if(idx % 30 == 0) {
                listOf<Genotype>(
                    GenotypeBuilder("sample1", listOf(Allele.REF_A)).make(),
                    GenotypeBuilder("sample2", listOf(Allele.ALT_C)).make(),
                    GenotypeBuilder("sample3", listOf(Allele.ALT_C)).make(),
                    GenotypeBuilder("sample4", listOf(Allele.REF_A)).make()
                )
            }
            else {
                listOf<Genotype>(
                    GenotypeBuilder("sample1", listOf(Allele.REF_A)).make(),
                    GenotypeBuilder("sample2", listOf(Allele.ALT_C)).make(),
                    GenotypeBuilder("sample3", listOf(Allele.REF_A)).make(),
                    GenotypeBuilder("sample4", listOf(Allele.ALT_C)).make()
                )
            }

            vcb.genotypes(genotypes)
            writer.add(vcb.make())


        }

        writer.close()

    }

    @Test
    fun testConvertCountMapsToData() {
        fail("Not yet implemented")
    }



}