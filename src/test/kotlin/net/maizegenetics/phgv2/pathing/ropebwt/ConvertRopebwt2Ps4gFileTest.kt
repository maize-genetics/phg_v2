package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File

class ConvertRopebwt2Ps4gFileTest {

    companion object {
        val tempTestDir = "${TestExtension.tempDir}ConvertRopeBwt2PS4GTestDir/"


        //Setup/download  files
        //Resetting on both setup and teardown just to be safe.
        @JvmStatic
        @BeforeAll
        fun setup() {
            resetDirs()
            setupDebugLogging()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            resetDirs()
        }

        private fun resetDirs() {

            File(TestExtension.tempDir).deleteRecursively()
            File(tempTestDir).deleteRecursively()

            File(TestExtension.tempDir).mkdirs()
            File(tempTestDir).mkdirs()
        }
    }

    @Test
    fun testCliktParams() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()

        val noBedFile = convertRopebwt2Ps4gFile.test("--output-dir testDir --spline-knot-dir ./knotFiles/")
        assertEquals(1, noBedFile.statusCode)
        assertEquals("Usage: convert-ropebwt2ps4g-file [<options>]\n\n" +
                "Error: must provide one of --ropebwt-bed, --ropebwt-bed-files, --ropebwt-bed-list-file, --ropebwt-bed-dir\n", noBedFile.stderr)

        val noOutputDir = convertRopebwt2Ps4gFile.test("--ropebwt-bed testDir --spline-knot-dir ./knotFiles/")
        assertEquals(1, noOutputDir.statusCode)
        assertEquals("Usage: convert-ropebwt2ps4g-file [<options>]\n\n" +
                "Error: missing option --output-dir\n", noOutputDir.stderr)

        val noHvcfDir = convertRopebwt2Ps4gFile.test("--ropebwt-bed testDir --output-dir testDir")
        assertEquals(1, noHvcfDir.statusCode)
        assertEquals("Usage: convert-ropebwt2ps4g-file [<options>]\n\n" +
                "Error: missing option --spline-knot-dir\n", noHvcfDir.stderr)

        val negativeRange = convertRopebwt2Ps4gFile.test("--ropebwt-bed testDir --output-dir testDir --spline-knot-dir ./knotFiles/ --max-range -2")
        assertEquals(1, negativeRange.statusCode)
        assertEquals("Usage: convert-ropebwt2ps4g-file [<options>]\n\n" +
                "Error: invalid value for --max-range: max range must be non-negative\n", negativeRange.stderr)
    }

    @Test
    fun testWriteOutput() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()
        val ropebwtBed = "data/test/ropebwt/LineA_FullChr.bed"
        val hvcfDir = "data/test/ropebwt/testHVCFs"

        SplineUtils.buildSplineKnots(hvcfDir, "hvcf", tempTestDir)

        convertRopebwt2Ps4gFile.test("--output-dir $tempTestDir --spline-knot-dir $tempTestDir --ropebwt-bed $ropebwtBed")
        assertEquals(File("$tempTestDir/LineA_FullChr.ps4g").exists(), true)

        convertRopebwt2Ps4gFile.test("--output-dir $tempTestDir/test_A.ps4g --spline-knot-dir $tempTestDir --ropebwt-bed $ropebwtBed")
        assertEquals(File("$tempTestDir/test_A.ps4g").exists(), true)
    }
    @Test
    fun testConvertCountMapToPS4GData() {
        val countMap = mapOf(Pair(Position("1",1), listOf(1, 2)) to 3, Pair(Position("1",2), listOf(3, 4)) to 5)
        val ps4gData = PS4GUtils.convertCountMapToPS4GData(countMap)
        assertEquals(2, ps4gData.size)
        assertEquals(PS4GData(listOf(1, 2), Position("1",1), 3), ps4gData[0])
        assertEquals(PS4GData(listOf(3, 4), Position("1",2), 5), ps4gData[1])
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

        val encodedPositions = listOf(Pair("chr1_gamete1", Position("1",5120)),Pair("chr1_gamete2", Position("1",7680)))


        val gameteToIdxMap = mapOf(Pair("gamete1", 1), Pair("gamete2", 2))
        val consensus = convertRopebwt2Ps4gFile.createConsensusPositionAndGametes(encodedPositions, gameteToIdxMap)

        assertEquals(6400, consensus.first.position)
        assertEquals(listOf(1, 2), consensus.second)

        //check if the encoded positions are the exact same
        val encodedPositionsSamePos = listOf(Pair("chr1_gamete1", Position("1",5120)), Pair("chr1_gamete2", Position("1",5120)))
        val consensusSamePos = convertRopebwt2Ps4gFile.createConsensusPositionAndGametes(encodedPositionsSamePos, gameteToIdxMap)

        assertEquals(5120, consensusSamePos.first.position)
        assertEquals(listOf(1, 2), consensusSamePos.second)

        //Try empty list
        val consensusEmpty = convertRopebwt2Ps4gFile.createConsensusPositionAndGametes(listOf(), gameteToIdxMap)
        assertEquals(-1, consensusEmpty.first.position)
        assertEquals(0, consensusEmpty.second.size)

    }

    @Test
    fun testLookupHitsToRefPosition() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()

        val splineKnotLookup = mutableMapOf<String,List<Triple<Int,String,Int>>>()

        splineKnotLookup["chr1_sample1"] = listOf(
            Triple(1, "chr1", 1),
            Triple(3, "chr1", 3),
            Triple(5, "chr1", 5),
            Triple(7, "chr1", 7),
            Triple(9, "chr1", 9)
        )
        splineKnotLookup["chr1_sample2"] = listOf(
            Triple(1, "chr1", 2),
            Triple(3, "chr1", 4),
            Triple(5, "chr1", 6),
            Triple(7, "chr1", 8),
            Triple(9, "chr1", 10)
        )
        splineKnotLookup["chr1_sample3"] = listOf(
            Triple(1, "chr1", 3),
            Triple(3, "chr1", 5),
            Triple(5, "chr1", 7),
            Triple(7, "chr1", 9),
            Triple(11, "chr1", 13)
        )

        val splineLookup = LinearLookupFunction(splineKnotLookup, mapOf("chr1" to 0))

        //single hit should return the position
        val singleHit = listOf(MEMHit("chr1_sample1", "+", 1))
        val encodedSingleHit = convertRopebwt2Ps4gFile.lookupHitsToRefPosition(singleHit, splineLookup)
        assertEquals(1, encodedSingleHit.size)
        assertEquals(Pair("chr1_sample1", Position("chr1",1)), encodedSingleHit[0])

        val singleHitPos4 = listOf(MEMHit("chr1_sample1", "+", 4))
        val encodedSingleHitPos4 = convertRopebwt2Ps4gFile.lookupHitsToRefPosition(singleHitPos4, splineLookup)
        assertEquals(1, encodedSingleHitPos4.size)
        assertEquals(Pair("chr1_sample1", Position("chr1",4)), encodedSingleHitPos4[0])

        //multiple hits should return the position from the spline
        // in this simple case they are 3,4,5
        val multipleHits = listOf(MEMHit("chr1_sample1", "+", 3), MEMHit("chr1_sample2", "+", 3), MEMHit("chr1_sample3", "+", 3))
        val encodedMultipleHits = convertRopebwt2Ps4gFile.lookupHitsToRefPosition(multipleHits, splineLookup)
        assertEquals(3, encodedMultipleHits.size)
        assertEquals(Pair("chr1_sample1", Position("chr1",3)), encodedMultipleHits[0])
        assertEquals(Pair("chr1_sample2", Position("chr1",4)), encodedMultipleHits[1])
        assertEquals(Pair("chr1_sample3", Position("chr1",5)), encodedMultipleHits[2])

        //multiple hits with a missing spline should return -1
        val missingSpline = listOf(MEMHit("chr1_sample1", "+", 3), MEMHit("chr1_sample2", "+", 3), MEMHit("chr1_sample3", "+", 10))
        val encodedMissingSpline = convertRopebwt2Ps4gFile.lookupHitsToRefPosition(missingSpline, splineLookup)
        assertEquals(3, encodedMissingSpline.size)
        assertEquals(Pair("chr1_sample1", Position("chr1",3)), encodedMissingSpline[0])
        assertEquals(Pair("chr1_sample2", Position("chr1",4)), encodedMissingSpline[1])

        val missingChrSpline = listOf(MEMHit("chr1_sample1", "+", 3), MEMHit("chr1_sample2", "+", 3), MEMHit("chr1_sample10", "+", 5))
        val encodedMissingChrSpline = convertRopebwt2Ps4gFile.lookupHitsToRefPosition(missingChrSpline, splineLookup)
        assertEquals(2, encodedMissingChrSpline.size)
        assertEquals(Pair("chr1_sample1", Position("chr1",3)), encodedMissingChrSpline[0])
        assertEquals(Pair("chr1_sample2", Position("chr1",4)), encodedMissingChrSpline[1])
    }

    @Test
    fun testProcessMemsForRead() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()
        //Start with having no passing hits
        val memList = listOf(
            MEM("read1", 1, 20, 2, listOf(MEMHit("chr1_sample1", "+", 1), MEMHit("chr1_sample2", "+", 2))),
            MEM("read1", 0, 18, 2, listOf(MEMHit("chr1_sample1", "+", 3), MEMHit("chr1_sample2", "+", 4)))
        )
        val chrIndexMap = mapOf(Pair("1", 0), Pair("2", 1))
        val gameteToIdxMap = mapOf(Pair("sample1", 0), Pair("sample2", 1))

        val emptySplines = LinearLookupFunction(emptyMap(), emptyMap())

        val noPassingHits = convertRopebwt2Ps4gFile.processMemsForRead(memList, emptySplines, 5, 1, 10, gameteToIdxMap)
        //Pair(-1, listOf())
        assertEquals(-1, noPassingHits.first.position)
        assertEquals(0, noPassingHits.second.size)

        val noPassingHits2 = convertRopebwt2Ps4gFile.processMemsForRead(memList, emptySplines, 50, 30, 10, gameteToIdxMap)
        //Also should be Pair(-1, listOf())
        assertEquals(-1, noPassingHits2.first.position)
        assertEquals(0, noPassingHits2.second.size)

        val noPassingHits3 = convertRopebwt2Ps4gFile.processMemsForRead(memList, emptySplines, 5, 30, 0, gameteToIdxMap)
        //Also should be Pair(-1, listOf())
        assertEquals(-1, noPassingHits3.first.position)
        assertEquals(0, noPassingHits3.second.size)

        val knots = buildSimpleKnotMap()

        val splineLookup =LinearLookupFunction(knots,chrIndexMap)

        val processedMems = convertRopebwt2Ps4gFile.processMemsForRead(memList, splineLookup, 19, 10, 10, gameteToIdxMap)

        assertEquals(2, processedMems.first.position) // 1 + 3 = 4 /2 = 2
        assertEquals(2, processedMems.second.size)
        assertEquals(0, processedMems.second[0])
        assertEquals(1, processedMems.second[1])

    }

    private fun buildSimpleKnotMap(): MutableMap<String, List<Triple<Int, String, Int>>> {
        val knots = mutableMapOf<String, List<Triple<Int, String, Int>>>()

        knots["chr1_sample1"] = listOf(
            Triple(1, "1", 1),
            Triple(3, "1", 3),
            Triple(5, "1", 5),
            Triple(7, "1", 7),
            Triple(9, "1", 9)
        )
        knots["chr1_sample2"] = listOf(
            Triple(1, "1", 2),
            Triple(3, "1", 4),
            Triple(5, "1", 6),
            Triple(7, "1", 8),
            Triple(9, "1", 10)
        )
        return knots
    }

    private fun buildSimpleKnotMapOffset(): MutableMap<String, List<Triple<Int, String, Int>>> {
        val knots = mutableMapOf<String, List<Triple<Int, String, Int>>>()

        knots["chr1_sample1"] = listOf(
            Triple(11, "1", 1),
            Triple(13, "1", 3),
            Triple(15, "1", 5),
            Triple(17, "1", 7),
            Triple(19, "1", 9)
        )
        knots["chr1_sample2"] = listOf(
            Triple(11, "1", 2),
            Triple(13, "1", 4),
            Triple(15, "1", 6),
            Triple(17, "1", 8),
            Triple(19, "1", 10)
        )
        return knots
    }

    private fun buildSimpleKnotMapReversed(): MutableMap<String, List<Triple<Int, String, Int>>> {
        val knots = mutableMapOf<String, List<Triple<Int, String, Int>>>()

        knots["chr1_sample1"] = listOf(
            Triple(1, "1", 9),
            Triple(3, "1", 7),
            Triple(5, "1", 5),
            Triple(7, "1", 3),
            Triple(9, "1", 1)
        )
        knots["chr1_sample2"] = listOf(
            Triple(1, "1", 2),
            Triple(3, "1", 4),
            Triple(5, "1", 6),
            Triple(7, "1", 8),
            Triple(9, "1", 10)
        )
        return knots
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

        val knots = buildSimpleKnotMap()
        val chrIndexMap = mapOf(Pair("1", 0), Pair("2", 1))
        val splineLookup = LinearLookupFunction(knots, chrIndexMap)


        val gameteToIdxMap = mapOf(Pair("sample1", 0), Pair("sample2", 1))

        val countMap = mutableMapOf<Pair<Position, List<Int>>, Int>()
        val sampleGameteCountMap = mutableMapOf<SampleGamete, Int>()
        val gameteIdxToSampleGameteMap = mapOf(Pair(0, SampleGamete("sample1", 0)), Pair(1, SampleGamete("sample2", 1)))

        convertRopebwt2Ps4gFile.processTempMEMs(tempMems1, splineLookup,  5, 10, 10, gameteToIdxMap, countMap, sampleGameteCountMap, gameteIdxToSampleGameteMap)
        assertEquals(1, countMap.size)
        assertEquals(1, countMap[Pair(Position("1",2), listOf(0,1))])

        convertRopebwt2Ps4gFile.processTempMEMs(tempMems2, splineLookup,  5, 10, 10, gameteToIdxMap, countMap, sampleGameteCountMap, gameteIdxToSampleGameteMap)
        assertEquals(2, countMap.size)
        assertEquals(1, countMap[Pair(Position( "1",1), listOf(0))])
        assertEquals(1, countMap[Pair(Position("1",2), listOf(0,1))])

        convertRopebwt2Ps4gFile.processTempMEMs(tempMems3, splineLookup,  19, 10, 10, gameteToIdxMap, countMap, sampleGameteCountMap, gameteIdxToSampleGameteMap)
        assertEquals(2, countMap.size)
        assertEquals(1, countMap[Pair(Position("1",1), listOf(0))])
        assertEquals(2, countMap[Pair(Position("1",2), listOf(0,1))])
    }

    @Test
    fun testBuildPS4GData() {
        val convertRopebwt2Ps4gFile = ConvertRopebwt2Ps4gFile()
        val ropebwtBed = "data/test/ropebwt/LineA_FullChr.bed"
        val hvcfDir = "data/test/ropebwt/testHVCFs"

        val truthDataWideRange=  setOf(PS4GData(listOf(0),Position("1",6), 2), PS4GData(listOf(0),Position("1",4),1),
            PS4GData(listOf(0),Position("1",8), 1), PS4GData(listOf(0),Position("1",12), 2)
        )

        val truthDataNarrowRange=  setOf(PS4GData(listOf(0),Position("1",6), 2), PS4GData(listOf(0),Position("1",4),1),
            PS4GData(listOf(0),Position("1",8), 1)
        )

        SplineUtils.buildSplineKnots(hvcfDir, "hvcf", tempTestDir)

        val (splineKnots, chrIndexMap, gameteToIdxMap) = SplineUtils.loadSplineKnotLookupFromDirectory(tempTestDir)

        val splineLookup = LinearLookupFunction(splineKnots, chrIndexMap)


        val ps4gData = convertRopebwt2Ps4gFile.buildPS4GData(ropebwtBed, splineLookup, gameteToIdxMap, 148, 10, 50)

        assertEquals(4, ps4gData.first.size)
        assertEquals(1, ps4gData.second.size)
        assertEquals(6, ps4gData.second[SampleGamete("LineA", 0)])
        for(data in ps4gData.first) {
            assertTrue(truthDataWideRange.contains(data))
        }

        val ps4gData2 = convertRopebwt2Ps4gFile.buildPS4GData(ropebwtBed, splineLookup, gameteToIdxMap, 148, 10, 15)

        assertEquals(3, ps4gData2.first.size)
        assertEquals(1, ps4gData2.second.size)
        assertEquals(4, ps4gData2.second[SampleGamete("LineA", 0)])
        for(data in ps4gData2.first) {
            assertTrue(truthDataNarrowRange.contains(data))
        }

        resetDirs()
    }

    @Test
    fun testLinearLookup() {
        val knots = buildSimpleKnotMap()
        val splineLookup = LinearLookupFunction(knots, mapOf("1" to 0, "2" to 1))

        //MEMHit("chr1_sample1", "+", 1), MEMHit("chr1_sample2", "+", 2))
        val lookup1 = splineLookup.value(Position("chr1_sample1", 1))
        assertNotNull(lookup1)

        assertEquals("1", lookup1.contig)
        assertEquals(1, lookup1.position)


        val lookup2 = splineLookup.value(Position("chr1_sample1", 2))
        assertNotNull(lookup2)
        assertEquals("1", lookup2.contig)
        assertEquals(2, lookup2.position)


        val lookupUnknownContig = splineLookup.value(Position("chr1_sample3", 1))
        assertEquals("unknown", lookupUnknownContig.contig)
        assertEquals(0, lookupUnknownContig.position)

        val lookupUnknownPosition = splineLookup.value(Position("chr1_sample1", 100))
        assertEquals("unknown", lookupUnknownPosition.contig)
        assertEquals(0, lookupUnknownPosition.position)

        val negativeKnots = buildSimpleKnotMapReversed()

        val negativeSplineLookup = LinearLookupFunction(negativeKnots, mapOf("1" to 0, "2" to 1))

        val lookup1Negative = negativeSplineLookup.value(Position("chr1_sample1", 1))
        assertNotNull(lookup1Negative)

        assertEquals("1", lookup1Negative.contig)
        assertEquals(9, lookup1Negative.position)


        val offsetMap = buildSimpleKnotMapOffset()
        val offsetSplineLookup = LinearLookupFunction(offsetMap, mapOf("1" to 0, "2" to 1))
        val lookup1Offset = offsetSplineLookup.value(Position("chr1_sample1", 11))
        assertNotNull(lookup1Offset)
        assertEquals("1", lookup1Offset.contig)
        assertEquals(1, lookup1Offset.position)
        //Check positions before the knots should be unknown
        val lookupBeforeKnots = offsetSplineLookup.value(Position("chr1_sample1", 10))
        assertEquals("unknown", lookupBeforeKnots.contig)
        assertEquals(0, lookupBeforeKnots.position)
        //Check positions after the knots should be unknown
        val lookupAfterKnots = offsetSplineLookup.value(Position("chr1_sample1", 20))
        assertEquals("unknown", lookupAfterKnots.contig)
        assertEquals(0, lookupAfterKnots.position)
    }


}