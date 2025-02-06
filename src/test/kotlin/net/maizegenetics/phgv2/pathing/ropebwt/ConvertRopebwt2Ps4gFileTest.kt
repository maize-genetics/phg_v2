package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.utils.Position
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test
import kotlin.test.fail

class ConvertRopebwt2Ps4gFileTest {

    @Test
    fun testCliktParams() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()

        val noBedFile = convertRopebwt2Ps4gFile.test("--output-dir testDir --hvcf-dir testDir")
        assertEquals(1, noBedFile.statusCode)
        assertEquals("Usage: convert-ropebwt2ps4g-file [<options>]\n\n" +
                "Error: missing option --ropebwt-bed\n", noBedFile.stderr)

        val noOutputDir = convertRopebwt2Ps4gFile.test("--ropebwt-bed testDir --hvcf-dir testDir")
        assertEquals(1, noOutputDir.statusCode)
        assertEquals("Usage: convert-ropebwt2ps4g-file [<options>]\n\n" +
                "Error: missing option --output-dir\n", noOutputDir.stderr)

        val noHvcfDir = convertRopebwt2Ps4gFile.test("--ropebwt-bed testDir --output-dir testDir")
        assertEquals(1, noHvcfDir.statusCode)
        assertEquals("Usage: convert-ropebwt2ps4g-file [<options>]\n\n" +
                "Error: missing option --hvcf-dir\n", noHvcfDir.stderr)

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
        fail("Not implemented")
    }

    //fun processMemsForRead(tempMems: List<MEM>, splineLookup: Map<String, PolynomialSplineFunction>,
    //                           chrIndexMap: Map<String, Int>, minMEMLength: Int, maxNumHits: Int,
    //                           gameteToIdxMap: Map<String, Int>): Pair<Int, List<Int>>
    @Test
    fun testProcessMemsForRead() {
        fail("Not implemented")
    }

    //processTempMEMs(
    //        tempMems: MutableList<MEM>,
    //        splineLookup: Map<String, PolynomialSplineFunction>,
    //        chrIndexMap: Map<String, Int>,
    //        minMEMLength: Int,
    //        maxNumHits: Int,
    //        gameteToIdxMap: Map<String, Int>,
    //        countMap: MutableMap<Pair<Int, List<Int>>, Int>,
    //        sampleGameteCountMap: MutableMap<SampleGamete, Int>,
    //        gameteIdxToSampleGameteMap: Map<Int, SampleGamete>
    //    )
    @Test
    fun testProcessTempMEMs() {
        fail("Not implemented")
    }

    //buildPS4GData(ropebwtBed: String,  splineLookup: Map<String, PolynomialSplineFunction>, chrIndexMap:Map<String,Int>,
    //                      gameteToIdxMap: Map<String,Int>,
    //                      minMEMLength: Int, maxNumHits: Int) : Pair<List<PS4GData>, Map<SampleGamete,Int>>
    @Test
    fun testBuildPS4GData() {
        fail("Not implemented")
    }

    //fun buildSpline(
    //        listOfPoints: MutableList<Pair<Double, Double>>,
    //        splineBuilder: SplineInterpolator,
    //        splineMap: MutableMap<String, PolynomialSplineFunction>,
    //        currentChrom: String,
    //        sampleName: String?
    //    )
    @Test
    fun testBuildSpline() {
        fail("Not implemented")
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

    //fun processHvcfFileIntoSplines(
    //        hvcfFile: File?,
    //        splineMap: MutableMap<String, PolynomialSplineFunction>,
    //        chrIndexMap : MutableMap<String,Int>,
    //        gameteIndexMap: MutableMap<String, Int>
    //    )
    @Test
    fun testProcessHvcfFileIntoSplines() {
        fail("Not implemented")
    }

    //buildSplineLookup(hvcfDir: String) : Triple<Map<String, PolynomialSplineFunction>, Map<String,Int>, Map<String,Int>>
    @Test
    fun testBuildSplineLookup() {
        fail("Not implemented")
    }


}