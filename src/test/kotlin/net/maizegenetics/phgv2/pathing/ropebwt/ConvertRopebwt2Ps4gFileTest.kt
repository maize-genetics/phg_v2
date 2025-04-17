package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.utils.Position
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction
import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test
import java.io.File

class ConvertRopebwt2Ps4gFileTest {

    @Test
    fun testCliktParams() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()

        val noBedFile = convertRopebwt2Ps4gFile.test("--output-dir testDir --vcf-dir testDir")
        assertEquals(1, noBedFile.statusCode)
        assertEquals("Usage: convert-ropebwt2ps4g-file [<options>]\n\n" +
                "Error: missing option --ropebwt-bed\n", noBedFile.stderr)

        val noOutputDir = convertRopebwt2Ps4gFile.test("--ropebwt-bed testDir --vcf-dir testDir")
        assertEquals(1, noOutputDir.statusCode)
        assertEquals("Usage: convert-ropebwt2ps4g-file [<options>]\n\n" +
                "Error: missing option --output-dir\n", noOutputDir.stderr)

        val noHvcfDir = convertRopebwt2Ps4gFile.test("--ropebwt-bed testDir --output-dir testDir")
        assertEquals(1, noHvcfDir.statusCode)
        assertEquals("Usage: convert-ropebwt2ps4g-file [<options>]\n\n" +
                "Error: missing option --vcf-dir\n", noHvcfDir.stderr)

    }

    @Test
    fun testConvertCountMapToPS4GData() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()
        val countMap = mapOf(Pair(1, listOf(1, 2)) to 3, Pair(2, listOf(3, 4)) to 5)
        val ps4gData = convertRopebwt2Ps4gFile.convertCountMapToPS4GData(countMap)
        assertEquals(2, ps4gData.size)
        assertEquals(PS4GData(listOf(1, 2), 1, 3), ps4gData[0])
        assertEquals(PS4GData(listOf(3, 4), 2, 5), ps4gData[1])
    }

    @Test
    fun testFindBestMems() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()
        val tempMems = listOf(
            MEM("read1", 1, 20, 2, listOf(MEMHit("contig1", "+", 1), MEMHit("contig2", "+", 2))),
            MEM("read1", 0, 19, 2, listOf(MEMHit("contig1", "+", 3), MEMHit("contig2", "+", 4)))
        )
        val bestMems = convertRopebwt2Ps4gFile.findBestMems(tempMems, 1, 5)
        assertEquals(4, bestMems.size)
        assertEquals(MEMHit("contig1", "+", 1), bestMems[0])
        assertEquals(MEMHit("contig2", "+", 2), bestMems[1])
        assertEquals(MEMHit("contig1", "+", 3), bestMems[2])
        assertEquals(MEMHit("contig2", "+", 4), bestMems[3])

        val bestMems2 = convertRopebwt2Ps4gFile.findBestMems(tempMems, 1, 2)
        assertEquals(0, bestMems2.size)

        val tempMemsShort = listOf(
            MEM("read1", 1, 20, 2, listOf(MEMHit("contig1", "+", 1), MEMHit("contig2", "+", 2))),
            MEM("read1", 0, 15, 2, listOf(MEMHit("contig1", "+", 3), MEMHit("contig2", "+", 4)))
        )

        val bestMemsNoShort = convertRopebwt2Ps4gFile.findBestMems(tempMemsShort, 17, 5)
        assertEquals(2, bestMemsNoShort.size)
        assertEquals(MEMHit("contig1", "+", 1), bestMemsNoShort[0])
        assertEquals(MEMHit("contig2", "+", 2), bestMemsNoShort[1])


        val emptyTempMems = listOf<MEM>()
        val bestMemsEmpty = convertRopebwt2Ps4gFile.findBestMems(emptyTempMems, 1, 5)
        assertEquals(0, bestMemsEmpty.size)

        val bestMemsAllTooShort = convertRopebwt2Ps4gFile.findBestMems(tempMemsShort, 25, 2)
        assertEquals(0, bestMemsAllTooShort.size)
    }

    @Test
    fun testCreateConsensusPositionAndGametes() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()
        val chrIndexMap = mapOf(Pair("chr1", 1), Pair("chr2", 2))
        val encodedPositions = listOf(Pair("chr1_gamete1", PS4GUtils.encodePosition(Position("chr1",20),chrIndexMap)), Pair("chr1_gamete2", PS4GUtils.encodePosition(Position("chr1",30),chrIndexMap)))
        val gameteToIdxMap = mapOf(Pair("gamete1", 1), Pair("gamete2", 2))
        val consensus = convertRopebwt2Ps4gFile.createConsensusPositionAndGametes(encodedPositions, chrIndexMap, gameteToIdxMap)

        assertEquals(PS4GUtils.encodePositionFromIdxAndPos(1,25), consensus.first)
        assertEquals(listOf(1, 2), consensus.second)

        //check if the encoded positions are the exact same
        val encodedPositionsSamePos = listOf(Pair("chr1_gamete1", PS4GUtils.encodePosition(Position("chr1",20),chrIndexMap)), Pair("chr1_gamete2", PS4GUtils.encodePosition(Position("chr1",20),chrIndexMap)))
        val consensusSamePos = convertRopebwt2Ps4gFile.createConsensusPositionAndGametes(encodedPositionsSamePos, chrIndexMap, gameteToIdxMap)

        assertEquals(PS4GUtils.encodePositionFromIdxAndPos(1,20), consensusSamePos.first)
        assertEquals(listOf(1, 2), consensusSamePos.second)

        //Try empty list
        val consensusEmpty = convertRopebwt2Ps4gFile.createConsensusPositionAndGametes(listOf(), chrIndexMap, gameteToIdxMap)
        assertEquals(-1, consensusEmpty.first)
        assertEquals(0, consensusEmpty.second.size)

    }

    @Test
    fun testEncodeHitsToPosition() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()
        val splineLookup = mutableMapOf<String, PolynomialSplineFunction>()
        val splineBuilder = AkimaSplineInterpolator()
        val listOfPoints = mutableListOf(
            Pair(1.0, 1.0),
            Pair(3.0, 3.0),
            Pair(5.0, 5.0),
            Pair(7.0, 7.0),
            Pair(9.0, 9.0)
        )
        convertRopebwt2Ps4gFile.buildSpline(listOfPoints, splineBuilder, splineLookup, "chr1", "sample1")

        val listOfPoints2 = mutableListOf(
            Pair(1.0, 2.0),
            Pair(3.0, 4.0),
            Pair(5.0, 6.0),
            Pair(7.0, 8.0),
            Pair(9.0, 10.0)
        )
        convertRopebwt2Ps4gFile.buildSpline(listOfPoints2, splineBuilder, splineLookup, "chr1", "sample2")

        val listOfPoints3 = mutableListOf(
            Pair(1.0, 3.0),
            Pair(3.0, 5.0),
            Pair(5.0, 7.0),
            Pair(7.0, 9.0),
            Pair(11.0, 13.0)
        )
        convertRopebwt2Ps4gFile.buildSpline(listOfPoints3, splineBuilder, splineLookup, "chr1", "sample3")

        //single hit should return the position
        val singleHit = listOf(MEMHit("chr1_sample1", "+", 1))
        val encodedSingleHit = convertRopebwt2Ps4gFile.encodeHitsToPosition(singleHit, splineLookup)
        assertEquals(1, encodedSingleHit.size)
        assertEquals(Pair("chr1_sample1", 1), encodedSingleHit[0])

        val singleHitPos4 = listOf(MEMHit("chr1_sample1", "+", 4))
        val encodedSingleHitPos4 = convertRopebwt2Ps4gFile.encodeHitsToPosition(singleHitPos4, splineLookup)
        assertEquals(1, encodedSingleHitPos4.size)
        assertEquals(Pair("chr1_sample1", 4), encodedSingleHitPos4[0])


        //multiple hits should return the position from the spline
        // in this simple case they are 3,4,5
        val multipleHits = listOf(MEMHit("chr1_sample1", "+", 3), MEMHit("chr1_sample2", "+", 3), MEMHit("chr1_sample3", "+", 3))
        val encodedMultipleHits = convertRopebwt2Ps4gFile.encodeHitsToPosition(multipleHits, splineLookup)
        assertEquals(3, encodedMultipleHits.size)
        assertEquals(Pair("chr1_sample1", 3), encodedMultipleHits[0])
        assertEquals(Pair("chr1_sample2", 4), encodedMultipleHits[1])
        assertEquals(Pair("chr1_sample3", 5), encodedMultipleHits[2])

        //multiple hits with a missing spline should return -1
        val missingSpline = listOf(MEMHit("chr1_sample1", "+", 3), MEMHit("chr1_sample2", "+", 3), MEMHit("chr1_sample3", "+", 10))
        val encodedMissingSpline = convertRopebwt2Ps4gFile.encodeHitsToPosition(missingSpline, splineLookup)
        assertEquals(3, encodedMissingSpline.size)
        assertEquals(Pair("chr1_sample1", 3), encodedMissingSpline[0])
        assertEquals(Pair("chr1_sample2", 4), encodedMissingSpline[1])

        val missingChrSpline = listOf(MEMHit("chr1_sample1", "+", 3), MEMHit("chr1_sample2", "+", 3), MEMHit("chr1_sample10", "+", 5))
        val encodedMissingChrSpline = convertRopebwt2Ps4gFile.encodeHitsToPosition(missingChrSpline, splineLookup)
        assertEquals(2, encodedMissingChrSpline.size)
        assertEquals(Pair("chr1_sample1", 3), encodedMissingChrSpline[0])
        assertEquals(Pair("chr1_sample2", 4), encodedMissingChrSpline[1])
    }

    @Test
    fun testProcessMemsForRead() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()
        //Start with having no passing hits
        val memList = listOf(
            MEM("read1", 1, 20, 2, listOf(MEMHit("chr1_sample1", "+", 1), MEMHit("chr1_sample2", "+", 2))),
            MEM("read1", 0, 18, 2, listOf(MEMHit("chr1_sample1", "+", 3), MEMHit("chr1_sample2", "+", 4)))
        )
        val chrIndexMap = mapOf(Pair("chr1", 0), Pair("chr2", 1))
        val gameteToIdxMap = mapOf(Pair("sample1", 0), Pair("sample2", 1))

        val noPassingHits = convertRopebwt2Ps4gFile.processMemsForRead(memList, emptyMap(), chrIndexMap, 5, 1, gameteToIdxMap)
        //Pair(-1, listOf())
        assertEquals(-1, noPassingHits.first)
        assertEquals(0, noPassingHits.second.size)

        val noPassingHits2 = convertRopebwt2Ps4gFile.processMemsForRead(memList, mapOf(), chrIndexMap, 50, 30, gameteToIdxMap)
        //Also should be Pair(-1, listOf())
        assertEquals(-1, noPassingHits2.first)
        assertEquals(0, noPassingHits2.second.size)


        //Now we test with making the hits pass
        val splineLookup = mutableMapOf<String, PolynomialSplineFunction>()
        val splineBuilder = AkimaSplineInterpolator()
        val listOfPoints = mutableListOf(
            Pair(1.0, 1.0),
            Pair(3.0, 3.0),
            Pair(5.0, 5.0),
            Pair(7.0, 7.0),
            Pair(9.0, 9.0)
        )
        convertRopebwt2Ps4gFile.buildSpline(listOfPoints, splineBuilder, splineLookup, "chr1", "sample1")

        val listOfPoints2 = mutableListOf(
            Pair(1.0, 2.0),
            Pair(3.0, 4.0),
            Pair(5.0, 6.0),
            Pair(7.0, 8.0),
            Pair(9.0, 10.0)
        )
        convertRopebwt2Ps4gFile.buildSpline(listOfPoints2, splineBuilder, splineLookup, "chr1", "sample2")

        val processedMems = convertRopebwt2Ps4gFile.processMemsForRead(memList, splineLookup, chrIndexMap, 19, 10, gameteToIdxMap)

        assertEquals(2, processedMems.first) // 1 + 3 = 4 /2 = 2
        assertEquals(2, processedMems.second.size)
        assertEquals(0, processedMems.second[0])
        assertEquals(1, processedMems.second[1])
    }

    @Test
    fun testProcessTempMEMs() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()
        val tempMems1 = mutableListOf(
            MEM("read1", 1, 20, 2, listOf(MEMHit("chr1_sample1", "+", 1), MEMHit("chr1_sample2", "+", 2)))
        )

        val tempMems2 = mutableListOf(
            MEM("read1", 1, 20, 2, listOf(MEMHit("chr1_sample1", "+", 1)))
        )

        val tempMems3 = mutableListOf(
            MEM("read1", 1, 20, 2, listOf(MEMHit("chr1_sample1", "+", 1), MEMHit("chr1_sample2", "+", 2))),
            MEM("read1", 0, 18, 2, listOf(MEMHit("chr1_sample1", "+", 3), MEMHit("chr1_sample2", "+", 4)))
        )

        // create a spline lookup (NOTE - need at least 5 observations so AkimaSplineInterpolator does not throw exception)
        val splineLookup = mutableMapOf<String, PolynomialSplineFunction>()
        val splineBuilder = AkimaSplineInterpolator()
        val listOfPoints = mutableListOf(
            Pair(1.0, 1.0),
            Pair(3.0, 3.0),
            Pair(5.0, 5.0),
            Pair(7.0, 7.0),
            Pair(9.0, 9.0)
        )
        convertRopebwt2Ps4gFile.buildSpline(listOfPoints, splineBuilder, splineLookup, "chr1", "sample1")

        val listOfPoints2 = mutableListOf(
            Pair(1.0, 2.0),
            Pair(3.0, 4.0),
            Pair(5.0, 6.0),
            Pair(7.0, 8.0),
            Pair(9.0, 10.0)
        )
        convertRopebwt2Ps4gFile.buildSpline(listOfPoints2, splineBuilder, splineLookup, "chr1", "sample2")

        val chrIndexMap = mapOf(Pair("chr1", 0), Pair("chr2", 1))
        val gameteToIdxMap = mapOf(Pair("sample1", 0), Pair("sample2", 1))

        val countMap = mutableMapOf<Pair<Int, List<Int>>, Int>()
        val sampleGameteCountMap = mutableMapOf<SampleGamete, Int>()
        val gameteIdxToSampleGameteMap = mapOf(Pair(0, SampleGamete("sample1", 0)), Pair(1, SampleGamete("sample2", 1)))

        convertRopebwt2Ps4gFile.processTempMEMs(tempMems1, splineLookup, chrIndexMap, 5, 10, gameteToIdxMap, countMap, sampleGameteCountMap, gameteIdxToSampleGameteMap)
        assertEquals(1, countMap.size)
        assertEquals(1, countMap[Pair(2, listOf(0,1))])

        convertRopebwt2Ps4gFile.processTempMEMs(tempMems2, splineLookup, chrIndexMap, 5, 10, gameteToIdxMap, countMap, sampleGameteCountMap, gameteIdxToSampleGameteMap)
        assertEquals(2, countMap.size)
        assertEquals(1, countMap[Pair(1, listOf(0))])
        assertEquals(1, countMap[Pair(2, listOf(0,1))])

        convertRopebwt2Ps4gFile.processTempMEMs(tempMems3, splineLookup, chrIndexMap, 19, 10, gameteToIdxMap, countMap, sampleGameteCountMap, gameteIdxToSampleGameteMap)
        assertEquals(2, countMap.size)
        assertEquals(1, countMap[Pair(1, listOf(0))])
        assertEquals(2, countMap[Pair(2, listOf(0,1))])
    }

    //buildPS4GData(ropebwtBed: String,  splineLookup: Map<String, PolynomialSplineFunction>, chrIndexMap:Map<String,Int>,
    //                      gameteToIdxMap: Map<String,Int>,
    //                      minMEMLength: Int, maxNumHits: Int) : Pair<List<PS4GData>, Map<SampleGamete,Int>>
    @Test
    fun testBuildPS4GData() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()
        val ropebwtBed = "data/test/ropebwt/LineA_FullChr.bed"
        val hvcfDir = "data/test/ropebwt/testHVCFs"

        val truthData=  setOf(PS4GData(listOf(0),6, 2), PS4GData(listOf(0),4,1),
            PS4GData(listOf(0),8, 1), PS4GData(listOf(0),12, 2)
        )

        val (splineLookup, chrIndexMap, gameteToIdxMap) = convertRopebwt2Ps4gFile.buildSplineLookup(hvcfDir, "hvcf")

        val ps4gData = convertRopebwt2Ps4gFile.buildPS4GData(ropebwtBed, splineLookup, chrIndexMap, gameteToIdxMap, 148, 10)



        assertEquals(4, ps4gData.first.size)
        assertEquals(1, ps4gData.second.size)
        assertEquals(6, ps4gData.second[SampleGamete("LineA", 0)])
        for(data in ps4gData.first) {
            assertTrue(truthData.contains(data))
        }

    }

    @Test
    fun testBuildSpline() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()
        //make a simple linear spline
        val listOfPoints = mutableListOf(
            Pair(1.0, 1.0),
            Pair(3.0, 3.0),
            Pair(5.0, 5.0),
            Pair(7.0, 7.0),
            Pair(9.0, 9.0)
        )
        val splineBuilder = AkimaSplineInterpolator()
        val splineMap = mutableMapOf<String, PolynomialSplineFunction>()

        convertRopebwt2Ps4gFile.buildSpline(listOfPoints, splineBuilder, splineMap, "chr1", "sample1")

        assertEquals(1, splineMap.size)
        assertTrue(splineMap.containsKey("chr1_sample1"))
        val spline = splineMap["chr1_sample1"]!!
        assertEquals(1.0, spline.value(1.0))
        assertEquals(3.0, spline.value(3.0))
        assertEquals(5.0, spline.value(5.0))
        assertEquals(2.0, spline.value(2.0), 0.0001)
        assertEquals(4.0, spline.value(4.0), 0.0001)
        assertFalse(spline.isValidPoint(10.0))


        //make a spline with multiple values for the same x
        //The code should remove one of them based on how it sees things
        //This list of points should be the same as the previous one
        val listOfPoints2 = mutableListOf(
            Pair(1.0, 1.0),
            Pair(3.0, 3.0),
            Pair(3.0, 4.0),
            Pair(5.0, 5.0),
            Pair(5.0, 6.0),
            Pair(7.0, 7.0),
            Pair(7.0, 9.0),
            Pair(9.0, 9.0),
            Pair(9.0, 11.0)
        )
        convertRopebwt2Ps4gFile.buildSpline(listOfPoints2, splineBuilder, splineMap, "chr1", "sample2")

        assertEquals(2, splineMap.size)
        assertTrue(splineMap.containsKey("chr1_sample2"))
        val spline2 = splineMap["chr1_sample2"]!!
        assertEquals(1.0, spline2.value(1.0))
        assertEquals(3.0, spline2.value(3.0))
        assertEquals(5.0, spline2.value(5.0))
        assertEquals(2.0, spline2.value(2.0), 0.0001)
        assertEquals(4.0, spline2.value(4.0), 0.0001)
        assertFalse(spline2.isValidPoint(10.0))


    }

    @Test
    fun testCheckMapAndAddToIndex() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()
        val stringToIndexMap = mutableMapOf(Pair("test1", 0), Pair("test2", 1))
        convertRopebwt2Ps4gFile.checkMapAndAddToIndex(stringToIndexMap, "test3")
        assertEquals(3, stringToIndexMap.size)
        assertEquals(2, stringToIndexMap["test3"])

        convertRopebwt2Ps4gFile.checkMapAndAddToIndex(stringToIndexMap, "test1")
        assertEquals(3, stringToIndexMap.size)
        assertEquals(0, stringToIndexMap["test1"])
    }

    @Test
    fun testProcessHvcfFileIntoSplines() {
        val inputFile = "data/test/ropebwt/testHVCFs/LineA.h.vcf"
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()
        val splineMap = mutableMapOf<String, PolynomialSplineFunction>()
        val chrIndexMap = mutableMapOf("1" to 0, "2" to 1)
        val gameteIndexMap = mutableMapOf("LineA" to 0, "LineB" to 1)

        convertRopebwt2Ps4gFile.processHvcfFileIntoSplines(File(inputFile), splineMap, chrIndexMap, gameteIndexMap)


        //Test some of the values in the spline
        val chr1Spline = splineMap["1_LineA"]!!

        println(chr1Spline.value(1500.0))
        println(chr1Spline.value(1500.0).toInt())
        println(PS4GUtils.decodePosition(chr1Spline.value(1500.0).toInt()))

        assertEquals(0, PS4GUtils.decodePosition(chr1Spline.value(1.0).toInt()).position)
        assertEquals(256, PS4GUtils.decodePosition(chr1Spline.value(256.0).toInt()).position)
        assertEquals(1024, PS4GUtils.decodePosition(chr1Spline.value(1500.0).toInt()).position) // 1500/256 = 5
        assertEquals(2560, PS4GUtils.decodePosition(chr1Spline.value(3000.0).toInt()).position) // 3000/256 = 11

        assertFalse(chr1Spline.isValidPoint(30000.0))
    }

    //buildSplineLookup(hvcfDir: String) : Triple<Map<String, PolynomialSplineFunction>, Map<String,Int>, Map<String,Int>>
    @Test
    fun testBuildSplineLookup() {
        val hvcfDir = "data/test/ropebwt/testHVCFs"
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()
        val (splineMap, chrIndexMap, gameteIndexMap) = convertRopebwt2Ps4gFile.buildSplineLookup(hvcfDir,"hvcf")

        assertEquals(2, splineMap.size)
        assertEquals(2, chrIndexMap.size)
        assertEquals(1, gameteIndexMap.size)

        val chr1Spline = splineMap["1_LineA"]!!
        assertEquals(0, PS4GUtils.decodePosition(chr1Spline.value(1.0).toInt()).position)
        assertEquals(256, PS4GUtils.decodePosition(chr1Spline.value(256.0).toInt()).position)
        assertEquals(1024, PS4GUtils.decodePosition(chr1Spline.value(1500.0).toInt()).position) // 1500/256 = 5
        assertEquals(2560, PS4GUtils.decodePosition(chr1Spline.value(3000.0).toInt()).position) // 3000/256 = 11

        assertFalse(chr1Spline.isValidPoint(30000.0))

    }


}