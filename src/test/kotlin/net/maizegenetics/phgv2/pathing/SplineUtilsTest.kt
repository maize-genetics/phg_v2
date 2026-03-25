package net.maizegenetics.phgv2.pathing

import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContextBuilder
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.pathing.ropebwt.IndexMaps
import net.maizegenetics.phgv2.pathing.ropebwt.LinearLookupFunction
import net.maizegenetics.phgv2.pathing.ropebwt.SplineUtils
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.setupDebugLogging
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
                assertEquals(key.split("_").first(), refChr)
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
        val chr1Knots = splineKnotLookup.splineKnotMap["1_LineA"]!!
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

        val lookup = LinearLookupFunction(splineKnotMap, chrIndexMap)


        //Check some values
        val pos1 = lookup.value(Position("1_LineA",1))
        assertEquals("1", pos1.contig)
        assertEquals(0, pos1.position)
        val pos2 = lookup.value(Position("1_LineA",256))
        assertEquals("1", pos2.contig)
        assertEquals(0, pos2.position)
        val pos3 = lookup.value(Position("1_LineA",3000))
        assertEquals("1", pos3.contig)
        assertEquals(10, pos3.position) // 3000/256 = 11
        val pos4 = lookup.value(Position("1_LineA",5000))
        assertEquals("1", pos4.contig)
        assertEquals(18, pos4.position) // 5000/256 = 19

        val unknown = lookup.value(Position("1_LineA",30000))
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

    @Test
    fun testFindAsmCoords() {
        //Need to test to make sure we handle each case of gvcf for both forward and reverse strand and for each version of the disableASMCoordinates
        //Testing refBlock
        //VariantContextBuilder(final String source, final String contig, final long start, final long stop, final Collection<Allele> alleles)
        val genotypeRefA = GenotypeBuilder("sample1").alleles(listOf(Allele.REF_A)).make()
        val genotypeAltG = GenotypeBuilder("sample1").alleles(listOf(Allele.ALT_G)).make()
        val genotypeAltGGG = GenotypeBuilder("sample1").alleles(listOf(Allele.create("GGG",false))).make()
        val refBlockVariantPos = VariantContextBuilder("src","chr1",100,200, listOf(Allele.REF_A, Allele.NON_REF_ALLELE))
            .attribute("End", 200)
            .attribute("ASM_Start", 300)
            .attribute("ASM_End", 400)
            .genotypes(genotypeRefA)
            .make()

        val asmBasedCoordsForRefBlockPos = SplineUtils.findAsmCoords(refBlockVariantPos,400)
        assertEquals(Triple(300,400,-1), asmBasedCoordsForRefBlockPos)

        val countBasedCoordsForRefBlockPos = SplineUtils.findAsmCoords(refBlockVariantPos, 400,true)
        assertEquals(Triple(400,500,501), countBasedCoordsForRefBlockPos)

        val refBlockVariantNeg = VariantContextBuilder("src","chr1",100,200, listOf(Allele.REF_A, Allele.NON_REF_ALLELE))
            .attribute("End", 200)
            .attribute("ASM_Start",400)
            .attribute("ASM_End",300)
            .genotypes(genotypeRefA)
            .make()

        val asmBasedCoordsForRefBlockNeg = SplineUtils.findAsmCoords(refBlockVariantNeg,400)
        assertEquals(Triple(400,300,-1), asmBasedCoordsForRefBlockNeg)

        val countBasedCoordsForRefBlockNeg = SplineUtils.findAsmCoords(refBlockVariantNeg, 400,true)
        assertEquals(Triple(400,500,501), countBasedCoordsForRefBlockNeg)

        //Testing SNP  SNP does not need to be tested for strand as it doesnt matter as its only a single bp
        val snpVariant = VariantContextBuilder("src","chr1",100,100, listOf(Allele.REF_A, Allele.ALT_G))
            .attribute("End", 100)
            .attribute("ASM_Start",300)
            .attribute("ASM_End",300)
            .genotypes(genotypeAltG)
            .make()

        val asmBasedCoordsForSNP = SplineUtils.findAsmCoords(snpVariant,400)
        assertEquals(Triple(300,300,-1), asmBasedCoordsForSNP)

        val countBasedCoordsForSNP = SplineUtils.findAsmCoords(snpVariant,400,true)
        assertEquals(Triple(400,400,401), countBasedCoordsForSNP)

        //Testing Insertion
        val insertionVariantPos = VariantContextBuilder("src","chr1",100,100, listOf(Allele.REF_G, Allele.create("GGG",false)))
            .attribute("End", 100)
            .attribute("ASM_Start",300)
            .attribute("ASM_End",302)
            .genotypes(genotypeAltGGG)
            .make()

        val asmBasedCoordsForInsertionPos = SplineUtils.findAsmCoords(insertionVariantPos,400)
        assertEquals(Triple(300,302,-1), asmBasedCoordsForInsertionPos)

        val countBasedCoordsForInsertionPos = SplineUtils.findAsmCoords(insertionVariantPos,400,true)
        assertEquals(Triple(400,402,403), countBasedCoordsForInsertionPos)

        val insertionVariantNeg = VariantContextBuilder("src","chr1",100,100, listOf(Allele.REF_G, Allele.create("GGG",false)))
            .attribute("End", 100)
            .attribute("ASM_Start",302)
            .attribute("ASM_End",300)
            .genotypes(genotypeAltGGG)
            .make()

        val asmBasedCoordsForInsertionNeg = SplineUtils.findAsmCoords(insertionVariantNeg,400)
        assertEquals(Triple(302,300,-1), asmBasedCoordsForInsertionNeg)

        val countBasedCoordsForInsertionNeg = SplineUtils.findAsmCoords(insertionVariantNeg,400,true)
        assertEquals(Triple(400,402,403), countBasedCoordsForInsertionNeg)

        //Testing Deletion Pos and Neg have the same alt size so we only have to check one.
        val deletionVar = VariantContextBuilder("src","chr1",100,102, listOf(Allele.create("AAA",true), Allele.ALT_G))
            .attribute("End", 103)
            .attribute("ASM_Start",300)
            .attribute("ASM_End",300)
            .genotypes(genotypeAltG)
            .make()

        val asmBasedCoordsForDeletion = SplineUtils.findAsmCoords(deletionVar,400)
        assertEquals(Triple(300,300,-1), asmBasedCoordsForDeletion)

        val countBasedCoordsForDeletion = SplineUtils.findAsmCoords(deletionVar,400, true)
        assertEquals(Triple(400,400,401), countBasedCoordsForDeletion)

        //Test missing  This will not have ASM coords as it doesnt make sense to have them if its missing...
        val refBlockVariantMissing = VariantContextBuilder("src","chr1",100,200, listOf(Allele.REF_A, Allele.NON_REF_ALLELE))
            .attribute("End", 200)
            .genotypes(GenotypeBuilder("sample1").alleles(listOf(Allele.NO_CALL)).make())
            .make()

        val asmBasedCoordsForMissing = SplineUtils.findAsmCoords(refBlockVariantMissing,400)
        assertEquals(Triple(null,null,-1), asmBasedCoordsForMissing)

        val countBasedCoordsForMissing = SplineUtils.findAsmCoords(refBlockVariantMissing,400, true)
        assertEquals(Triple(null,null,400), countBasedCoordsForMissing)
    }

    @Test
    fun testDisableASMCoords() {
        val lineAFile = File("data/test/buildSplineKnots/LineA.g.vcf")
        val chrIndexMap = mutableMapOf<String, Int>()
        val gameteIndexMap = mutableMapOf<String, Int>()
        val splineKnotLookup = SplineUtils.processGvcfFileIntoSplineKnots(lineAFile, chrIndexMap, gameteIndexMap, minIndelLength = 1 ,disableASMCoordinates = true, binSize = 1)

        //Verify the output.
        val knots = splineKnotLookup.splineKnotMap["1_LineA"]!!
        assertEquals(8, knots.size)
        assertEquals(Triple(1,"1",1), knots[0])
        assertEquals(Triple(102,"1",102), knots[1])
        assertEquals(Triple(103,"1",112), knots[2]) //This seems odd but we need to midpoint the indel
        assertEquals(Triple(104,"1",123), knots[3])
        assertEquals(Triple(111,"1",130), knots[4])
        assertEquals(Triple(121,"1",131), knots[5]) //This seems odd but we need to midpoint the indel
        assertEquals(Triple(132,"1",132), knots[6])
        assertEquals(Triple(135,"1",135), knots[7])

        //Test unsorted errors
        val lineAFileUnsortedPos = File("data/test/buildSplineKnots/LineA_unsortedPos.g.vcf")
        //this should throw an illegal state exception
        assertThrows(IllegalStateException::class.java) {
            SplineUtils.processGvcfFileIntoSplineKnots(lineAFileUnsortedPos, chrIndexMap, gameteIndexMap, minIndelLength = 1 ,disableASMCoordinates = true, binSize = 1)
        }

        val lineAFileUnsortedChrom = File("data/test/buildSplineKnots/LineA_unsortedChr.g.vcf")
        assertThrows(IllegalStateException::class.java) {
            SplineUtils.processGvcfFileIntoSplineKnots(lineAFileUnsortedChrom, chrIndexMap, gameteIndexMap, minIndelLength = 1 ,disableASMCoordinates = true, binSize = 1)
        }
    }
}