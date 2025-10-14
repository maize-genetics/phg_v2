package net.maizegenetics.phgv2.pathing

import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.pathing.ropebwt.IndexMaps
import net.maizegenetics.phgv2.pathing.ropebwt.LinearLookupFunction
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

        //check to see if the maps have exactly what we expect
        //We do not need to test convertKnotsToSpline as it needs to be tested as part of the LinearLookupFunction

        assertEquals(2, splineKnotLookup.splineKnotMap.size)
        assertEquals(2, chrIndexMap.size)
        assertEquals(chrIndexMap["1"]!!, splineKnotLookup.chrIndexMap["1"])
        assertEquals(chrIndexMap["2"]!!, splineKnotLookup.chrIndexMap["2"])
        assertEquals(2, gameteIndexMap.size)
        assertEquals(gameteIndexMap["LineA"]!!, splineKnotLookup.gameteIndexMap["LineA"])
        assertEquals(gameteIndexMap["LineB"]!!, splineKnotLookup.gameteIndexMap["LineB"])

        for((key, value) in splineKnotLookup.splineKnotMap) {
            for((asmPos,refChr,refPos) in value) {
                assertEquals(key, refChr)
                assertEquals(asmPos/256, refPos)
            }
        }
    }

    @Test
    fun testProcessGvcfFileIntoSplines() {
        val inputFile = "data/test/smallseq/LineA.g.vcf"
        val chrIndexMap = mutableMapOf("1" to 0, "2" to 1)
        val gameteIndexMap = mutableMapOf("LineA" to 0, "LineB" to 1)

        val splineKnotLookup = SplineUtils.processGvcfFileIntoSplineKnots(File(inputFile), chrIndexMap, gameteIndexMap)

        assertEquals(2, splineKnotLookup.splineKnotMap.size)
        assertEquals(2, chrIndexMap.size)
        assertEquals(chrIndexMap["1"]!!, splineKnotLookup.chrIndexMap["1"])
        assertEquals(chrIndexMap["2"]!!, splineKnotLookup.chrIndexMap["2"])
        assertEquals(2, gameteIndexMap.size)
        assertEquals(gameteIndexMap["LineA"]!!, splineKnotLookup.gameteIndexMap["LineA"])
        assertEquals(gameteIndexMap["LineB"]!!, splineKnotLookup.gameteIndexMap["LineB"])


        //These have been verified manually
        val chr1Knots = splineKnotLookup.splineKnotMap["1"]!!
        assertEquals(chr1Knots[0],Triple(1, "1", 0))
        assertEquals(chr1Knots[1],Triple(3106, "1", 12))
        assertEquals(chr1Knots[2],Triple(3357, "1", 12))
        assertEquals(chr1Knots[4],Triple(3893, "1", 13))
    }

    @Test
    fun testBuildSplineLookup() {
        val hvcfDir = "data/test/ropebwt/testHVCFs"
        //(vcfDir: String, vcfType: String, outputDir: String ,minIndelLength: Int = 10, maxNumPointsPerChrom: Int = 250_000, contigSet : Set<String> = emptySet(), randomSeed: Long = 12345)
        SplineUtils.buildSplineKnots(hvcfDir,"hvcf", tempTestDir)

        val (splineKnotMap, chrIndexMap, gameteIndexMap) = SplineUtils.loadSplineKnotLookupFromDirectory(tempTestDir)

        assertEquals(2, splineKnotMap.size)
        assertEquals(2, chrIndexMap.size)
        assertEquals(1, gameteIndexMap.size)

        val lookup = LinearLookupFunction(splineKnotMap)


        //Check some values
        val pos1 = lookup.value(Position("1",1))
        assertEquals("1", pos1.contig)
        assertEquals(0, pos1.position)
        val pos2 = lookup.value(Position("1",256))
        assertEquals("1", pos2.contig)
        assertEquals(0, pos2.position)
        val pos3 = lookup.value(Position("1",3000))
        assertEquals("1", pos3.contig)
        assertEquals(10, pos3.position) // 3000/256 = 11
        val pos4 = lookup.value(Position("1",5000))
        assertEquals("1", pos4.contig)
        assertEquals(18, pos4.position) // 5000/256 = 19

        val unknown = lookup.value(Position("1",30000))
        assertEquals("unknown", unknown.contig)
        assertEquals(0, unknown.position)

        resetDirs()
    }

    @Test
    fun testSerializingSplineLookup() {
        val hvcfDir = "data/test/ropebwt/testHVCFs"

        SplineUtils.buildSplineKnots(hvcfDir,"hvcf", tempTestDir)

        val (splineKnotMap, chrIndexMap, gameteIndexMap) = SplineUtils.loadSplineKnotLookupFromDirectory(tempTestDir)

        //Need to do a reset otherwise we have too many files and the map gets too big
        resetDirs()

        val outputSplineFile = "${tempTestDir}Sample1_spline_knots.json.gz"
        val outputIndexFile = "${tempTestDir}index_maps.json.gz"

        SplineUtils.writeSplineKnotsToFile(splineKnotMap, outputSplineFile)
        SplineUtils.writeIndexMapsToFile(IndexMaps(chrIndexMap, gameteIndexMap), outputIndexFile)

        //Read the file back in and check that the values are the same
        val (splineMap2, chrIndexMap2, gameteIndexMap2) = SplineUtils.loadSplineKnotLookupFromDirectory(tempTestDir)

        assertEquals(splineKnotMap.size, splineMap2.size)

        //check that the entries of the maps are the same
        for(key in splineKnotMap.keys) {
            assertTrue(splineMap2.containsKey(key))
            val list1 = splineKnotMap[key]!!
            val list2 = splineMap2[key]!!
            assertEquals(list1.size, list2.size)
            for(i in list1.indices) {
                assertEquals(list1[i], list2[i])
            }
        }

        assertEquals(chrIndexMap.size, chrIndexMap2.size)
        assertEquals(gameteIndexMap.size, gameteIndexMap2.size)

        resetDirs()
    }

    @Test
    fun testDownsamplePoints() {
        val points = mutableListOf<Triple<Int,String,Int>>()
        for (i in 0 until 1000) {
            points.add(Triple(i,"chr1", i))
        }

        val points2 = mutableListOf<Triple<Int,String,Int>>()
        for (i in 0 until 99) {
            points2.add(Triple(i,"chr2", i))
        }

        val splineKnotLookup = mutableMapOf("chr1" to points, "chr2" to points2)

        assertEquals(1000, splineKnotLookup["chr1"]!!.size)
        assertEquals(99, splineKnotLookup["chr2"]!!.size)

        SplineUtils.downsamplePointsByChrLength(splineKnotLookup, 100)

        assertEquals(9, splineKnotLookup["chr1"]!!.size)
        assertEquals(99, splineKnotLookup["chr2"]!!.size)

    }

    @Test
    fun testDownsamplePointsByChrLength() {
        //downsamplePointsByChrLength(splineKnotMap:MutableMap<String, MutableList<Pair<Double,Double>>>, numBpsPerKnot: Int = 50_000, randomSeed : Long = 12345)
        //make a spline map with random increasing values
        val splineKnotMap = mutableMapOf<String, MutableList<Triple<Int,String,Int>>>()
        val rand = java.util.Random(12345)
        for (i in 1..5) {
            val points = mutableListOf<Triple<Int,String,Int>>()
            for (j in 0 until 10_000) {
                points.add(Triple(j * i, "$i",rand.nextInt() * 1000 + i * 1000))
            }
            splineKnotMap["$i"] = points
        }
        //Add one for 0 that has 999 points
        splineKnotMap["0"] = mutableListOf<Triple<Int,String,Int>>()
        for (j in 0 until 999) {
            splineKnotMap["0"]!!.add(Triple(j,"0", rand.nextInt() * 1000))
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
            }
        }

        //Make a spline map with no knots
        val emptySplineKnotMap = mutableMapOf<String, MutableList<Triple<Int,String,Int>>>()
        emptySplineKnotMap["0"] = mutableListOf<Triple<Int,String,Int>>()
        SplineUtils.downsamplePointsByChrLength(emptySplineKnotMap, 1000, 12345)
        //Check that the empty spline map is still empty
        assertEquals(1, emptySplineKnotMap.size)
        assertEquals(0, emptySplineKnotMap["0"]!!.size)


    }
}