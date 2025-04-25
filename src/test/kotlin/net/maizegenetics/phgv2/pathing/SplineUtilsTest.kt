package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.pathing.ropebwt.PS4GUtils
import net.maizegenetics.phgv2.pathing.ropebwt.SplineUtils
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction
import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test
import java.io.File

class SplineUtilsTest {

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
        val splineBuilder = AkimaSplineInterpolator()
        val splineMap = mutableMapOf<String, PolynomialSplineFunction>()

        SplineUtils.buildSpline(listOfPoints, splineBuilder, splineMap, "chr1", "sample1")

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
        SplineUtils.buildSpline(listOfPoints2, splineBuilder, splineMap, "chr1", "sample2")

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
        val splineMap = mutableMapOf<String, PolynomialSplineFunction>()
        val chrIndexMap = mutableMapOf("1" to 0, "2" to 1)
        val gameteIndexMap = mutableMapOf("LineA" to 0, "LineB" to 1)

        SplineUtils.processHvcfFileIntoSplines(File(inputFile), splineMap, chrIndexMap, gameteIndexMap)


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

    @Test
    fun testBuildSplineLookup() {
        val hvcfDir = "data/test/ropebwt/testHVCFs"
        val (splineMap, chrIndexMap, gameteIndexMap) = SplineUtils.buildSplineLookup(hvcfDir,"hvcf")

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