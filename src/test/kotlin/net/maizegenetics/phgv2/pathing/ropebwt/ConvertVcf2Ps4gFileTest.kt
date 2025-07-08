package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.createGenericHeader
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.Disabled
import org.junit.jupiter.api.Test
import java.io.File

class ConvertVcf2Ps4gFileTest {

    @Test
    fun testCliktParams() {
        val convertVcf2Ps4gFile = ConvertVcf2Ps4gFile()

        val noSampleVCF = convertVcf2Ps4gFile.test("--ref-panel-vcf testDir.vcf --output-dir testDir")
        assertEquals(1, noSampleVCF.statusCode)
        assertEquals("Usage: convert-vcf2ps4g-file [<options>]\n\n" +
                "Error: missing option --to-impute-vcf\n", noSampleVCF.stderr)

        val noGameteVCF = convertVcf2Ps4gFile.test("--to-impute-vcf testDir.vcf --output-dir testDir")
        assertEquals(1, noGameteVCF.statusCode)
        assertEquals("Usage: convert-vcf2ps4g-file [<options>]\n\n" +
                "Error: missing option --ref-panel-vcf\n", noGameteVCF.stderr)

        val noOutputDir = convertVcf2Ps4gFile.test("--to-impute-vcf testDir.vcf --ref-panel-vcf testDir.vcf")
        assertEquals(1, noOutputDir.statusCode)
        assertEquals("Usage: convert-vcf2ps4g-file [<options>]\n\n" +
                "Error: missing option --output-dir\n", noOutputDir.stderr)
    }


    @Test
    fun testCreateRefPanelPositionSampleGameteLookup() {
        val inputVCFFile = File("data/test/ps4gTests/refPanel.vcf")
        val convertVcf2Ps4gFile = ConvertVcf2Ps4gFile()
        val positionSampleGameteLookup = convertVcf2Ps4gFile.createRefPanelPositionSampleGameteLookup(inputVCFFile.path)

        val sample1Gamete = SampleGamete("sample1", 0)
        val sample2Gamete = SampleGamete("sample2", 0)
        val sample3Gamete = SampleGamete("sample3", 0)
        val sample4Gamete = SampleGamete("sample4", 0)


        assertEquals(400, positionSampleGameteLookup.size)
        for(i in 0 until 400) {
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
        //Load in the first variant from the refPanel
        val refVCFFile = "data/test/ps4gTests/refPanel.vcf"
        val inputVCFFile = "data/test/ps4gTests/imputeFile.vcf"
        val convertVcf2Ps4gFile = ConvertVcf2Ps4gFile()

        //fun createPS4GData(sampleVcf: String, positionSampleGameteLookup: Map<Position, Map<String, List<SampleGamete>>>): Map<SampleGamete, PS4GDataWithMaps> {
        val positionSampleGameteLookup = convertVcf2Ps4gFile.createRefPanelPositionSampleGameteLookup(refVCFFile)
        val ps4GData = convertVcf2Ps4gFile.createPS4GData(inputVCFFile, positionSampleGameteLookup)

        assertEquals(4, ps4GData.size)

        //check the counts
        //Each imputeSample should have 3 ps4gData entries
        val impute1PS4GData = ps4GData[SampleGamete("impute1",0)]!!
        assertEquals(3, impute1PS4GData.ps4gData.size)

        val impute2PS4GData = ps4GData[SampleGamete("impute2",0)]!!
        assertEquals(3, impute2PS4GData.ps4gData.size)

        val impute3PS4GData = ps4GData[SampleGamete("impute3",0)]!!
        assertEquals(3, impute3PS4GData.ps4gData.size)

        val impute4PS4GData = ps4GData[SampleGamete("impute4",0)]!!
        assertEquals(3, impute4PS4GData.ps4gData.size)


        //Check the ps4g entries
        val impute1PS4GList = impute1PS4GData.ps4gData
        assertEquals(listOf(0,1) ,impute1PS4GList[0].gameteList)
        assertEquals(0, impute1PS4GList[0].pos)
        assertEquals(3, impute1PS4GList[0].count)

        assertEquals(listOf(0,3) ,impute1PS4GList[1].gameteList)
        assertEquals(0, impute1PS4GList[1].pos)
        assertEquals(1, impute1PS4GList[1].count)

        assertEquals(listOf(0,2) ,impute1PS4GList[2].gameteList)
        assertEquals(2, impute1PS4GList[2].pos)
        assertEquals(1, impute1PS4GList[2].count)

        //check impute2
        val impute2PS4GList = impute2PS4GData.ps4gData
        assertEquals(listOf(0,1) ,impute2PS4GList[0].gameteList)
        assertEquals(0, impute2PS4GList[0].pos)
        assertEquals(3, impute2PS4GList[0].count)

        assertEquals(listOf(1,2) ,impute2PS4GList[1].gameteList)
        assertEquals(0, impute2PS4GList[1].pos)
        assertEquals(1, impute2PS4GList[1].count)

        assertEquals(listOf(1,3) ,impute2PS4GList[2].gameteList)
        assertEquals(2, impute2PS4GList[2].pos)
        assertEquals(1, impute2PS4GList[2].count)

        //check impute3
        val impute3PS4GList = impute3PS4GData.ps4gData
        assertEquals(listOf(2,3) ,impute3PS4GList[0].gameteList)
        assertEquals(0, impute3PS4GList[0].pos)
        assertEquals(3, impute3PS4GList[0].count)

        assertEquals(listOf(1,2) ,impute3PS4GList[1].gameteList)
        assertEquals(0, impute3PS4GList[1].pos)
        assertEquals(1, impute3PS4GList[1].count)

        assertEquals(listOf(0,2) ,impute3PS4GList[2].gameteList)
        assertEquals(2, impute3PS4GList[2].pos)
        assertEquals(1, impute3PS4GList[2].count)


        //check impute4
        val impute4PS4GList = impute4PS4GData.ps4gData
        assertEquals(listOf(2,3) ,impute4PS4GList[0].gameteList)
        assertEquals(0, impute4PS4GList[0].pos)
        assertEquals(3, impute4PS4GList[0].count)

        assertEquals(listOf(0,3) ,impute4PS4GList[1].gameteList)
        assertEquals(0, impute4PS4GList[1].pos)
        assertEquals(1, impute4PS4GList[1].count)

        assertEquals(listOf(1,3) ,impute4PS4GList[2].gameteList)
        assertEquals(2, impute4PS4GList[2].pos)
        assertEquals(1, impute4PS4GList[2].count)



    }

    @Test
    fun testProcessVariantPosition() {
        //Load in the first variant from the refPanel
        val refVCFFile = "data/test/ps4gTests/refPanel.vcf"
        val inputVCFFile = "data/test/ps4gTests/imputeFile.vcf"
        val convertVcf2Ps4gFile = ConvertVcf2Ps4gFile()

        VCFFileReader(File(inputVCFFile), false).use { vcfReader ->
            val firstVC = vcfReader.iterator().next()

            val position = Position(firstVC.contig, firstVC.start)

            val contigNameToIdxMap = mapOf("chr1" to 0, "chr2" to 1, "chr3" to 2)
            val sampleGameteCountMap = mutableMapOf<SampleGamete, MutableMap<SampleGamete, Int>>()
            val gameteToCountMap = mutableMapOf<SampleGamete, SampleGameteCountMaps>()
            val gameteToIdxMap = mutableMapOf<SampleGamete, Int>(SampleGamete("sample1", 0) to 0,
                SampleGamete("sample2", 0) to 1,
                SampleGamete("sample3", 0) to 2,
                SampleGamete("sample4", 0) to 3,
            )

            convertVcf2Ps4gFile.processVariantPosition(position, contigNameToIdxMap,
                convertVcf2Ps4gFile.createRefPanelPositionSampleGameteLookup(refVCFFile),
                firstVC, sampleGameteCountMap, gameteToCountMap, gameteToIdxMap)

            val encodedPos = PS4GUtils.encodePositionFromIdxAndPos(0,firstVC.start)
            //check that the gameteToCountMap is correct
            assertEquals(4, gameteToCountMap.size)
            assertEquals(4, sampleGameteCountMap.size)
            assertEquals(1, gameteToCountMap[SampleGamete("impute1",0)]?.countMap?.size)
            assertTrue(gameteToCountMap[SampleGamete("impute1",0)]?.countMap?.containsKey(Pair(encodedPos, listOf(0,1))) ?: false)
            assertEquals(1, gameteToCountMap[SampleGamete("impute1",0)]?.countMap?.get(Pair(encodedPos, listOf(0,1))))

            assertEquals(1, gameteToCountMap[SampleGamete("impute2",0)]?.countMap?.size)
            assertTrue(gameteToCountMap[SampleGamete("impute2",0)]?.countMap?.containsKey(Pair(encodedPos, listOf(0,1))) ?: false)
            assertEquals(1, gameteToCountMap[SampleGamete("impute2",0)]?.countMap?.get(Pair(encodedPos, listOf(0,1))) ?: 0)

            assertEquals(1, gameteToCountMap[SampleGamete("impute3",0)]?.countMap?.size)
            assertTrue(gameteToCountMap[SampleGamete("impute3",0)]?.countMap?.containsKey(Pair(encodedPos, listOf(2,3))) ?: false)
            assertEquals(1, gameteToCountMap[SampleGamete("impute3",0)]?.countMap?.get(Pair(encodedPos, listOf(2,3))) ?: 0)

            assertEquals(1, gameteToCountMap[SampleGamete("impute4",0)]?.countMap?.size)
            assertTrue(gameteToCountMap[SampleGamete("impute4",0)]?.countMap?.containsKey(Pair(encodedPos, listOf(2,3))) ?: false)
            assertEquals(1, gameteToCountMap[SampleGamete("impute4",0)]?.countMap?.get(Pair(encodedPos, listOf(2,3))) ?: 0)

            assertEquals(2, sampleGameteCountMap[SampleGamete("impute1",0)]?.size)
            assertEquals(2, sampleGameteCountMap[SampleGamete("impute2",0)]?.size)
            assertEquals(2, sampleGameteCountMap[SampleGamete("impute3",0)]?.size)
            assertEquals(2, sampleGameteCountMap[SampleGamete("impute4",0)]?.size)
        }
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

    @Disabled
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

        (0 until 4000 step 10).forEach { idx ->
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
        //fun convertCountMapsToData(countMaps: Map<SampleGamete, SampleGameteCountMaps>, gameteToIdxMap: Map<SampleGamete, Int>)
        // : Map<SampleGamete, PS4GDataWithMaps>

        val convertVcf2Ps4gFile = ConvertVcf2Ps4gFile()

        val gameteToIdxMap = mapOf<SampleGamete, Int>(
            SampleGamete("sample1", 0) to 0,
            SampleGamete("sample2", 0) to 1,
            SampleGamete("sample3", 0) to 2,
            SampleGamete("sample4", 0) to 3,
        )
        val imputeSample1CountMap = mutableMapOf<Pair<Int, List<Int>>, Int>(
            Pair(0, listOf(0, 1)) to 20,
            Pair(1, listOf(0)) to 10,
            Pair(2, listOf(0, 2)) to 5,
        )
        val imputeSample1RefCount = mutableMapOf<SampleGamete, Int>(
            SampleGamete("sample1", 0) to 35,
            SampleGamete("sample2", 0) to 20,
            SampleGamete("sample3", 0) to 5
        )

        val imputeSample2CountMap = mutableMapOf<Pair<Int, List<Int>>, Int>(
            Pair(0, listOf(0,)) to 10,
            Pair(1, listOf(0,1)) to 5,
            Pair(2, listOf(1, 2)) to 15,
        )

        val imputeSample2RefCount = mutableMapOf<SampleGamete, Int>(
            SampleGamete("sample1", 0) to 15,
            SampleGamete("sample2", 0) to 20,
            SampleGamete("sample3", 0) to 15
        )

        val sampleGameteCountMap = mapOf<SampleGamete, SampleGameteCountMaps>(
            SampleGamete("imputeSample1", 0) to SampleGameteCountMaps(
                imputeSample1CountMap, imputeSample1RefCount
            ),
            SampleGamete("imputeSample2", 0) to SampleGameteCountMaps(
                imputeSample2CountMap, imputeSample2RefCount
            )
        )

        val ps4gData = convertVcf2Ps4gFile.convertCountMapsToData(sampleGameteCountMap, gameteToIdxMap)

        assertEquals(2, ps4gData.size)

        //(val ps4gData: List<PS4GData>, val sampleGameteCount: Map<SampleGamete,Int>, val gameteToIdxMap: Map<SampleGamete,Int>)
        val impute1PS4GData = ps4gData[SampleGamete("imputeSample1",0)]!!
        val ps4gList = impute1PS4GData.ps4gData
        assertEquals(3, ps4gList.size)
        assertEquals(listOf(0,1), ps4gList[0].gameteList)
        assertEquals(0, ps4gList[0].pos)
        assertEquals(20, ps4gList[0].count)
        assertEquals(listOf(0), ps4gList[1].gameteList)
        assertEquals(1, ps4gList[1].pos)
        assertEquals(10, ps4gList[1].count)
        assertEquals(listOf(0,2), ps4gList[2].gameteList)
        assertEquals(2, ps4gList[2].pos)
        assertEquals(5, ps4gList[2].count)
        assertEquals(35, impute1PS4GData.sampleGameteCount[SampleGamete("sample1",0)])
        assertEquals(20, impute1PS4GData.sampleGameteCount[SampleGamete("sample2",0)])
        assertEquals(5, impute1PS4GData.sampleGameteCount[SampleGamete("sample3",0)])

        val impute2PS4GData = ps4gData[SampleGamete("imputeSample2",0)]!!
        val ps4gList2 = impute2PS4GData.ps4gData
        assertEquals(3, ps4gList2.size)
        assertEquals(listOf(0), ps4gList2[0].gameteList)
        assertEquals(0, ps4gList2[0].pos)
        assertEquals(10, ps4gList2[0].count)
        assertEquals(listOf(0,1), ps4gList2[1].gameteList)
        assertEquals(1, ps4gList2[1].pos)
        assertEquals(5, ps4gList2[1].count)
        assertEquals(listOf(1,2), ps4gList2[2].gameteList)
        assertEquals(2, ps4gList2[2].pos)
        assertEquals(15, impute2PS4GData.sampleGameteCount[SampleGamete("sample1",0)])
        assertEquals(20, impute2PS4GData.sampleGameteCount[SampleGamete("sample2",0)])
        assertEquals(15, impute2PS4GData.sampleGameteCount[SampleGamete("sample3",0)])




    }



}