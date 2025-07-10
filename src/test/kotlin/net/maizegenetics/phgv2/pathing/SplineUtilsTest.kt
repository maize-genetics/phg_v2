package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.pathing.ropebwt.IndexMaps
import net.maizegenetics.phgv2.pathing.ropebwt.PS4GUtils
import net.maizegenetics.phgv2.pathing.ropebwt.SplineKnotLookup
import net.maizegenetics.phgv2.pathing.ropebwt.SplineUtils
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.BeforeEach
import org.junit.jupiter.api.Test
import java.io.File

class SplineUtilsTest {

    companion object {
        val tempTestDir = "${TestExtension.tempDir}splineTest/"


        //Setup/download  files
        //Resetting on both setup and teardown just to be safe.
        @JvmStatic
        @BeforeAll
        fun setupBeforeAll() {
            resetDirs()
            setupDebugLogging()
        }


//        @JvmStatic
        @BeforeEach
        fun setupBeforeEach() {
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
    fun testBuildSpline() {
        //make a simple linear spline
        val listOfPoints = mutableListOf(
            Pair(1.0, 1.0),
            Pair(3.0, 3.0),
            Pair(5.0, 5.0),
            Pair(7.0, 7.0),
            Pair(9.0, 9.0)
        )

        val splineKnotMap = mutableMapOf<String, Pair<DoubleArray, DoubleArray>>()

        SplineUtils.buildSplineKnotsForASMChrom(listOfPoints, splineKnotMap, "chr1", "sample1")

        val splineMap = SplineUtils.convertKnotsToSpline(splineKnotMap)

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
        //Add more points to the existing knots and build new splines
        SplineUtils.buildSplineKnotsForASMChrom(listOfPoints2, splineKnotMap, "chr1", "sample2")

        val splineMap2 = SplineUtils.convertKnotsToSpline(splineKnotMap)


        assertEquals(2, splineMap2.size)
        assertTrue(splineMap2.containsKey("chr1_sample2"))
        val spline2 = splineMap2["chr1_sample2"]!!
        assertEquals(1.0, spline2.value(1.0))
        assertEquals(3.0, spline2.value(3.0))
        assertEquals(5.0, spline2.value(5.0))
        assertEquals(2.0, spline2.value(2.0), 0.0001)
        assertEquals(4.0, spline2.value(4.0), 0.0001)
        assertFalse(spline2.isValidPoint(10.0))


    }

    @Test
    fun testCheckMapAndAddToIndex() {
        val stringToIndexMap = mutableMapOf(Pair("test1", 0), Pair("test2", 1))
        SplineUtils.checkMapAndAddToIndex(stringToIndexMap, "test3")
        assertEquals(3, stringToIndexMap.size)
        assertEquals(2, stringToIndexMap["test3"])

        SplineUtils.checkMapAndAddToIndex(stringToIndexMap, "test1")
        assertEquals(3, stringToIndexMap.size)
        assertEquals(0, stringToIndexMap["test1"])
    }


    @Test
    fun testProcessHvcfFileIntoSplines() {
        val inputFile = "data/test/ropebwt/testHVCFs/LineA.h.vcf"
        val chrIndexMap = mutableMapOf("1" to 0, "2" to 1)
        val gameteIndexMap = mutableMapOf("LineA" to 0, "LineB" to 1)

        val splineKnotLookup = SplineUtils.processHvcfFileIntoSplineKnots(File(inputFile), chrIndexMap, gameteIndexMap)

        val splineMap = SplineUtils.convertKnotsToSpline(splineKnotLookup.splineKnotMap)


        //Test some of the values in the spline
        val chr1Spline = splineMap["1_LineA"]!!

        assertEquals(0, PS4GUtils.decodePosition(chr1Spline.value(1.0).toInt()).position)
        assertEquals(256, PS4GUtils.decodePosition(chr1Spline.value(256.0).toInt()).position)
        assertEquals(1024, PS4GUtils.decodePosition(chr1Spline.value(1500.0).toInt()).position) // 1500/256 = 5
        assertEquals(2560, PS4GUtils.decodePosition(chr1Spline.value(3000.0).toInt()).position) // 3000/256 = 11

        assertFalse(chr1Spline.isValidPoint(30000.0))
    }

    @Test
    fun testProcessGvcfFileIntoSplines() {
        val inputFile = "data/test/smallseq/LineA.g.vcf"
        val chrIndexMap = mutableMapOf("1" to 0, "2" to 1)
        val gameteIndexMap = mutableMapOf("LineA" to 0, "LineB" to 1)

        val splineKnotLookup = SplineUtils.processGvcfFileIntoSplineKnots(File(inputFile), chrIndexMap, gameteIndexMap)

        val splineMap = SplineUtils.convertKnotsToSpline(splineKnotLookup.splineKnotMap)


        //Test some of the values in the spline
        val chr1Spline = splineMap["1_LineA"]!!

        assertEquals(0, PS4GUtils.decodePosition(chr1Spline.value(1.0).toInt()).position)
        assertEquals(256, PS4GUtils.decodePosition(chr1Spline.value(256.0).toInt()).position)
        assertEquals(2048, PS4GUtils.decodePosition(chr1Spline.value(1500.0).toInt()).position) // 1500/256 = 5
        assertEquals(2816, PS4GUtils.decodePosition(chr1Spline.value(3000.0).toInt()).position) // 3000/256 = 11

        assertFalse(chr1Spline.isValidPoint(3000000.0))
        resetDirs()
    }

    @Test
    fun testBuildSplineLookup() {
        val hvcfDir = "data/test/ropebwt/testHVCFs"
        //(vcfDir: String, vcfType: String, outputDir: String ,minIndelLength: Int = 10, maxNumPointsPerChrom: Int = 250_000, contigSet : Set<String> = emptySet(), randomSeed: Long = 12345)
        SplineUtils.buildSplineKnots(hvcfDir,"hvcf", tempTestDir)

        val (splineKnotMap, chrIndexMap, gameteIndexMap) = SplineUtils.loadSplineKnotLookupFromDirectory(tempTestDir)

        val splineMap = SplineUtils.convertKnotsToSpline(splineKnotMap)

        assertEquals(2, splineMap.size)
        assertEquals(2, chrIndexMap.size)
        assertEquals(1, gameteIndexMap.size)

        val chr1Spline = splineMap["1_LineA"]!!
        assertEquals(0, PS4GUtils.decodePosition(chr1Spline.value(1.0).toInt()).position)
        assertEquals(256, PS4GUtils.decodePosition(chr1Spline.value(256.0).toInt()).position)
        assertEquals(1024, PS4GUtils.decodePosition(chr1Spline.value(1500.0).toInt()).position) // 1500/256 = 5
        assertEquals(2560, PS4GUtils.decodePosition(chr1Spline.value(3000.0).toInt()).position) // 3000/256 = 11

        assertFalse(chr1Spline.isValidPoint(30000.0))

        resetDirs()
    }

    @Test
    fun testSerializingSplineLookup() {
        val hvcfDir = "data/test/ropebwt/testHVCFs"

        SplineUtils.buildSplineKnots(hvcfDir,"hvcf", tempTestDir)

        val (splineKnotMap, chrIndexMap, gameteIndexMap) = SplineUtils.loadSplineKnotLookupFromDirectory(tempTestDir)

        val splineMap = SplineUtils.convertKnotsToSpline(splineKnotMap)

        val splineArrays = mutableMapOf<String, Pair<DoubleArray,DoubleArray>>()
        val rand = java.util.Random(123345)
        for ((key, spline) in splineMap) {
            val x = DoubleArray(200_000)
            val y = DoubleArray(200_000)
            for (i in 0 until 200_000) {
                x[i] = rand.nextDouble()
                y[i] = rand.nextDouble()
            }
            splineArrays[key] = Pair(x, y)
        }

        //Need to do a reset otherwise we have too many files and the map gets too big
        resetDirs()

        val outputSplineFile = "${tempTestDir}Sample1_spline_knots.json.gz"
        val outputIndexFile = "${tempTestDir}index_maps.json.gz"

        SplineUtils.writeSplineKnotsToFile(splineArrays, outputSplineFile)
        SplineUtils.writeIndexMapsToFile(IndexMaps(chrIndexMap, gameteIndexMap), outputIndexFile)



        //Read the file back in and check that the values are the same
        val (splineMap2, chrIndexMap2, gameteIndexMap2) = SplineUtils.loadSplineKnotLookupFromDirectory(tempTestDir)

        assertEquals(splineArrays.size, splineMap2.size)
        //check the entries of the arrays are the same
        for (key in splineArrays.keys) {
            assertTrue(splineMap2.containsKey(key))
            val (x, y) = splineArrays[key]!!
            val (x2, y2) = splineMap2[key]!!

            assertEquals(x.size, x2.size)
            assertEquals(y.size, y2.size)

            //Check x vs y sizes as well
            assertEquals(x.size, y.size)
            assertEquals(x2.size, y2.size)

            for (i in x.indices) {
                assertEquals(x[i], x2[i])
                assertEquals(y[i], y2[i])
            }
        }
        assertEquals(chrIndexMap.size, chrIndexMap2.size)
        assertEquals(gameteIndexMap.size, gameteIndexMap2.size)

        resetDirs()
    }

    @Test
    fun testDownsamplePoints() {
        val points = mutableListOf<Pair<Double, Double>>()
        for (i in 0 until 1000) {
            points.add(Pair(i.toDouble(), i.toDouble()))
        }

        val points2 = mutableListOf<Pair<Double, Double>>()
        for (i in 0 until 99) {
            points2.add(Pair(i.toDouble(), i.toDouble()))
        }

        val splineKnotLookup = mutableMapOf<String, MutableList<Pair<Double,Double>>>("chr1" to points, "chr2" to points2)

        assertEquals(1000, splineKnotLookup["chr1"]!!.size)
        assertEquals(99, splineKnotLookup["chr2"]!!.size)

        SplineUtils.downsamplePoints(splineKnotLookup, 100)

        assertEquals(100, splineKnotLookup["chr1"]!!.size)
        assertEquals(99, splineKnotLookup["chr2"]!!.size)

    }

    @Test
    fun testSortAndSplitPoints() {
        val firstPoint = Pair(1.0, PS4GUtils.encodePositionNoLookup(Position("0",1)).toDouble())
        val secondPoint = Pair(10001.0, PS4GUtils.encodePositionNoLookup(Position("0",10001)).toDouble())

        val points = mutableListOf(firstPoint, secondPoint)

        val sortedPoints = SplineUtils.sortAndSplitPoints(points)

        assertEquals(5, sortedPoints.size)
        assertEquals(0, PS4GUtils.decodePosition(sortedPoints[0].second.toInt()).position)
        assertEquals(0, PS4GUtils.decodePosition(sortedPoints[0].second.toInt()).contig.toInt())

        for(i in 1 until sortedPoints.size) {
            assertEquals(0, PS4GUtils.decodePosition(sortedPoints[i].second.toInt()).contig.toInt())
            // We ask for 4 points so they should be 2500, 5000, 7500, 10000.
            // We bin them by dividing by 256 and multiplying by 256 to get the binned position
            assertEquals(((i * 2500)/256) * 256, PS4GUtils.decodePosition(sortedPoints[i].second.toInt()).position)
        }

        //Test the negative strand too
        val firstPointNegative = Pair(1.0, PS4GUtils.encodePositionNoLookup(Position("1",10001)).toDouble())
        val secondPointNegative = Pair(10001.0, PS4GUtils.encodePositionNoLookup(Position("1",1)).toDouble())
        val pointsNegative = mutableListOf(firstPointNegative, secondPointNegative )

        val sortedPointsNegative = SplineUtils.sortAndSplitPoints(pointsNegative)


        assertEquals(5, sortedPointsNegative.size)
        assertEquals(1, PS4GUtils.decodePosition(sortedPointsNegative[0].second.toInt()).contig.toInt())
        assertEquals((10001/256) * 256, PS4GUtils.decodePosition(sortedPointsNegative[0].second.toInt()).position)


        for(i in 1 until sortedPointsNegative.size) {
            assertEquals(1, PS4GUtils.decodePosition(sortedPointsNegative[i].second.toInt()).contig.toInt())
            // We ask for 4 points so they should be 10000. 7500, 5000, 2500, 0
            // We bin them by dividing by 256 and multiplying by 256 to get the binned position
            assertEquals(((10001 - (i * 2500))/256) * 256, PS4GUtils.decodePosition(sortedPointsNegative[i].second.toInt()).position)
        }

    }

    @Test
    fun testDownsamplePointsByChrLength() {
        //downsamplePointsByChrLength(splineKnotMap:MutableMap<String, MutableList<Pair<Double,Double>>>, numBpsPerKnot: Int = 50_000, randomSeed : Long = 12345)
        //make a spline map with random increasing values
        val splineKnotMap = mutableMapOf<String, MutableList<Pair<Double, Double>>>()
        val rand = java.util.Random(12345)
        for (i in 1..5) {
            val points = mutableListOf<Pair<Double, Double>>()
            for (j in 0 until 10_000) {
                points.add(Pair(j.toDouble() * i, rand.nextDouble() * 1000 + i * 1000))
            }
            splineKnotMap["$i"] = points
        }
        //Add one for 0 that has 999 points
        splineKnotMap["0"] = mutableListOf<Pair<Double, Double>>()
        for (j in 0 until 999) {
            splineKnotMap["0"]!!.add(Pair(j.toDouble(), rand.nextDouble() * 1000))
        }

        assertEquals(6, splineKnotMap.size)
        for (key in splineKnotMap.keys) {
            if(key == "0") {
                //The 0 chromosome should have 999 points
                assertEquals(999, splineKnotMap[key]!!.size)
            } else {
                //The other chromosomes should have 10_000 points
                assertEquals(10_000, splineKnotMap[key]!!.size)
            }
        }
        SplineUtils.downsamplePointsByChrLength(splineKnotMap, 1000, 12345)
        //Check that the number of points is reduced
        assertEquals(6, splineKnotMap.size)
        for (key in splineKnotMap.keys) {
            if(key == "0") {
                //The 0 chromosome should have 999 points
                assertEquals(999, splineKnotMap[key]!!.size)
            } else {

                //The number of points should be reduced to 1000
                assertTrue(
                    splineKnotMap[key]!!.size <= (key.toInt() * 10),
                    "Spline map for $key has more than 1000 points: ${splineKnotMap[key]!!.size}"
                )
//            //Check that the points are still increasing
//            for (j in 1 until splineKnotMap[key]!!.size) {
//                assertTrue(splineKnotMap[key]!![j].first > splineKnotMap[key]!![j - 1].first, "Points are not increasing for $key at index $j")
//            }
            }
        }

    }
}