package net.maizegenetics.phgv2.utils

import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.AgcCompress
import net.maizegenetics.phgv2.cli.CreateMafVcf
import net.maizegenetics.phgv2.cli.Initdb
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.String
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class CreateMafVcfUtilsTest {

    @Test
    fun testMergeConsecutiveRegions() {
        val regions = mutableListOf(Pair(Position("chr1",1),Position("chr1",10)),
            Pair(Position("chr1",11),Position("chr1",20)),
            Pair(Position("chr1",21),Position("chr1",30)),
            Pair(Position("chr1",31),Position("chr1",40)),
            Pair(Position("chr1",41),Position("chr1",50)),
            Pair(Position("chr1",51),Position("chr1",60)),
            Pair(Position("chr1",61),Position("chr1",70)),
            Pair(Position("chr1",71),Position("chr1",80)),
            Pair(Position("chr1",81),Position("chr1",90)),
            Pair(Position("chr1",91),Position("chr1",100)),
            Pair(Position("chr1",201),Position("chr1",210)),
            Pair(Position("chr1",211),Position("chr1",220)),
            Pair(Position("chr1",221),Position("chr1",230)),
            Pair(Position("chr1",231),Position("chr1",240)),
            Pair(Position("chr1",241),Position("chr1",250)),
            Pair(Position("chr1",251),Position("chr1",260)),
            Pair(Position("chr1",261),Position("chr1",270)),
            Pair(Position("chr1",271),Position("chr1",280)))

        val mergedRegions = CreateMafVcfUtils.mergeConsecutiveRegions(regions)

        assertEquals(2, mergedRegions.size)
        assertEquals(Pair(Position("chr1",1),Position("chr1",100)), mergedRegions[0])
        assertEquals(Pair(Position("chr1",201),Position("chr1",280)), mergedRegions[1])


        //make inverted regions
        val invertedRegions = mutableListOf(Pair(Position("chr1",280),Position("chr1",271)),
            Pair(Position("chr1",270),Position("chr1",261)),
            Pair(Position("chr1",260),Position("chr1",251)),
            Pair(Position("chr1",250),Position("chr1",241)),
            Pair(Position("chr1",240),Position("chr1",231)),
            Pair(Position("chr1",230),Position("chr1",221)),
            Pair(Position("chr1",220),Position("chr1",211)),
            Pair(Position("chr1",210),Position("chr1",201)),
            Pair(Position("chr1",100),Position("chr1",91)),
            Pair(Position("chr1",90),Position("chr1",81)),
            Pair(Position("chr1",80),Position("chr1",71)),
            Pair(Position("chr1",70),Position("chr1",61)),
            Pair(Position("chr1",60),Position("chr1",51)),
            Pair(Position("chr1",50),Position("chr1",41)),
            Pair(Position("chr1",40),Position("chr1",31)),
            Pair(Position("chr1",30),Position("chr1",21)),
            Pair(Position("chr1",20),Position("chr1",11)),
            Pair(Position("chr1",10),Position("chr1",1)))

        val mergedInvertedRegions = CreateMafVcfUtils.mergeConsecutiveRegions(invertedRegions)
        assertEquals(2, mergedInvertedRegions.size)
        assertEquals(Pair(Position("chr1",280),Position("chr1",201)), mergedInvertedRegions[0])
        assertEquals(Pair(Position("chr1",100),Position("chr1",1)), mergedInvertedRegions[1])


        //test mixed positive then inverted regions
        val mixedRegions = mutableListOf(Pair(Position("chr1",1),Position("chr1",10)),
            Pair(Position("chr1",11),Position("chr1",20)),
            Pair(Position("chr1",21),Position("chr1",30)),
            Pair(Position("chr1",31),Position("chr1",40)),
            Pair(Position("chr1",41),Position("chr1",50)),
            Pair(Position("chr1",51),Position("chr1",60)),
            Pair(Position("chr1",61),Position("chr1",70)),
            Pair(Position("chr1",71),Position("chr1",80)),
            Pair(Position("chr1",81),Position("chr1",90)),
            Pair(Position("chr1",91),Position("chr1",100)),
            Pair(Position("chr1",280),Position("chr1",271)),
            Pair(Position("chr1",270),Position("chr1",261)),
            Pair(Position("chr1",260),Position("chr1",251)),
            Pair(Position("chr1",250),Position("chr1",241)),
            Pair(Position("chr1",240),Position("chr1",231)),
            Pair(Position("chr1",230),Position("chr1",221)),
            Pair(Position("chr1",220),Position("chr1",211)),
            Pair(Position("chr1",210),Position("chr1",201)))

        val mergedMixedRegions = CreateMafVcfUtils.mergeConsecutiveRegions(mixedRegions)
        assertEquals(2, mergedMixedRegions.size)
        assertEquals(Pair(Position("chr1",1),Position("chr1",100)), mergedMixedRegions[0])
        assertEquals(Pair(Position("chr1",280),Position("chr1",201)), mergedMixedRegions[1])


        //test mixed inverted then positive regions
        val mixedRegions2 = mutableListOf(Pair(Position("chr1",280),Position("chr1",271)),
            Pair(Position("chr1",270),Position("chr1",261)),
            Pair(Position("chr1",260),Position("chr1",251)),
            Pair(Position("chr1",250),Position("chr1",241)),
            Pair(Position("chr1",240),Position("chr1",231)),
            Pair(Position("chr1",230),Position("chr1",221)),
            Pair(Position("chr1",220),Position("chr1",211)),
            Pair(Position("chr1",210),Position("chr1",201)),
            Pair(Position("chr1",1),Position("chr1",10)),
            Pair(Position("chr1",11),Position("chr1",20)),
            Pair(Position("chr1",21),Position("chr1",30)),
            Pair(Position("chr1",31),Position("chr1",40)),
            Pair(Position("chr1",41),Position("chr1",50)),
            Pair(Position("chr1",51),Position("chr1",60)),
            Pair(Position("chr1",61),Position("chr1",70)),
            Pair(Position("chr1",71),Position("chr1",80)),
            Pair(Position("chr1",81),Position("chr1",90)),
            Pair(Position("chr1",91),Position("chr1",100)))

        val mergedMixedRegions2 = CreateMafVcfUtils.mergeConsecutiveRegions(mixedRegions2)
        assertEquals(2, mergedMixedRegions2.size)
        assertEquals(Pair(Position("chr1",280),Position("chr1",201)), mergedMixedRegions2[0])
        assertEquals(Pair(Position("chr1",1),Position("chr1",100)), mergedMixedRegions2[1])


        //Test merging when they are consecutive but switch strands
        val mixedRegions3 = mutableListOf(Pair(Position("chr1",1),Position("chr1",10)),
            Pair(Position("chr1",11),Position("chr1",20)),
            Pair(Position("chr1",21),Position("chr1",30)),
            Pair(Position("chr1",31),Position("chr1",40)),
            Pair(Position("chr1",41),Position("chr1",50)),
            Pair(Position("chr1",51),Position("chr1",60)),
            Pair(Position("chr1",61),Position("chr1",70)),
            Pair(Position("chr1",71),Position("chr1",80)),
            Pair(Position("chr1",81),Position("chr1",90)),
            Pair(Position("chr1",91),Position("chr1",100)),
            Pair(Position("chr1",150),Position("chr1",141)),
            Pair(Position("chr1",140),Position("chr1",121)),
            Pair(Position("chr1",120),Position("chr1",101)))
        val mergedMixedRegions3 = CreateMafVcfUtils.mergeConsecutiveRegions(mixedRegions3)
        assertEquals(2, mergedMixedRegions3.size)
        assertEquals(Pair(Position("chr1",1),Position("chr1",100)), mergedMixedRegions3[0])
        assertEquals(Pair(Position("chr1",150),Position("chr1",101)), mergedMixedRegions3[1])


        //have inverted regions that are single bp regions
        val singleRegions = mutableListOf(Pair(Position("chr1",1),Position("chr1",30)),
            Pair(Position("chr1",31),Position("chr1",31)),
            Pair(Position("chr1",32),Position("chr1",60)))

        val mergedSingleRegions = CreateMafVcfUtils.mergeConsecutiveRegions(singleRegions)
        assertEquals(1, mergedSingleRegions.size)
        assertEquals(Pair(Position("chr1",1),Position("chr1",60)), mergedSingleRegions[0])

        val singleRegionsInverted = mutableListOf(Pair(Position("chr1",60),Position("chr1",32)),
            Pair(Position("chr1",31),Position("chr1",31)),
            Pair(Position("chr1",30),Position("chr1",1)))

        val mergedSingleRegionsInverted = CreateMafVcfUtils.mergeConsecutiveRegions(singleRegionsInverted)
        assertEquals(1, mergedSingleRegionsInverted.size)
        assertEquals(Pair(Position("chr1",60),Position("chr1",1)), mergedSingleRegionsInverted[0])
    }

    @Test
    fun testResizePositionRange() {
        val range = Pair(Position("chr1",5), Position("chr1",20))
        val resizedRange = CreateMafVcfUtils.resizePositionRange(range, 15, true)

        assertEquals(Position("chr1",15), resizedRange.first)
        assertEquals(Position("chr1",20), resizedRange.second)

        val resizedRange2 = CreateMafVcfUtils.resizePositionRange(range, 15, false)

        assertEquals(Position("chr1",5), resizedRange2.first)
        assertEquals(Position("chr1",15), resizedRange2.second)

    }

    @Test
    fun testRetrieveSequenceForRanges() {
        //Need to build an AGC archive
        //write the assemblies list
        val assembliesList = "data/test/smallseq/assembliesList.txt"
        getBufferedWriter(assembliesList).use { myWriter ->
            myWriter.write("data/test/smallseq/LineA.fa\n")
            myWriter.write("data/test/smallseq/LineB.fa\n")
        }

        File(TestExtension.testTileDBURI).mkdirs()
        val initdb = Initdb()
        initdb.test("--db-path ${TestExtension.testTileDBURI}")
        //Create the agc record:
        val agcCompress = AgcCompress()

        val agcResult = agcCompress.test("--fasta-list $assembliesList --db-path ${TestExtension.testTileDBURI} --reference-file ${TestExtension.smallseqRefFile}")
        println(agcResult.output)


//        retrieveSequenceForRanges(dbPath: String, ranges: List<String>, condaEnvPrefix:String = "", numRangesPerAgcQuery: Int = -1)
        //Build 10 ranges to pull out
        val ranges = listOf<String>( "1@LineA:1-10", "1@LineA:11-20", "1@LineA:21-30", "1@LineA:31-40", "1@LineA:41-50", "1@LineA:51-60", "1@LineA:61-70", "1@LineA:71-80", "1@LineA:81-90", "1@LineA:91-100" )




        val fullChrom = CreateMafVcfUtils.retrieveSequenceForRanges(TestExtension.testTileDBURI, ranges)

        val splitChrom = CreateMafVcfUtils.retrieveSequenceForRanges(TestExtension.testTileDBURI, ranges, numRangesPerAgcQuery = 3) //Doing 3 so its uneven


        //Check to see if the results are the same
        assertEquals(fullChrom, splitChrom)


        //Delete the agc archive
        File(assembliesList).delete()
        File("${TestExtension.testTileDBURI}assemblies.agc").delete()
    }

}