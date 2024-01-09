package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.random.Random
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class MostLikelyParentsTest {

    @Test
    fun testLikelyParents() {
        //build a haplotype groph from Ref, lineA, and lineB
        //for that need an hvcfdir with the files in it
        //make sure the test output directory exists
        File(TestExtension.testOutputDir).mkdirs()

        val myGraph = HaplotypeGraph(listOf(TestExtension.smallseqRefHvcfFile, TestExtension.smallseqLineAHvcfFile,
            TestExtension.smallseqLineBHvcfFile))

        //some readMapping counts
        val refrangeToReadCounts = createReadMappings(myGraph)

        val parentList = MostLikelyParents(myGraph).findMostLikelyParents(refrangeToReadCounts, 3, 1.0)
        assertEquals(2, parentList.size)
        assertEquals("LineA", parentList[0].first.name)
        assertEquals("LineB", parentList[1].first.name)
        assertEquals(190, parentList[0].second)
        assertEquals(38, parentList[1].second)

    }

    private fun createReadMappings(graph: HaplotypeGraph): Map<ReferenceRange, Map<List<String>, Int>> {
        val probMapToA = 0.7
        val probMapToB = 0.5
        val probMapToRef = 0.1

        val sampleGameteA = SampleGamete("LineA")
        val sampleGameteB = SampleGamete("LineB")
        val sampleGameteRef = SampleGamete("Ref")

        val readCounts = mutableMapOf<ReferenceRange, Map<List<String>, Int>>()
        for (range in graph.ranges()) {
            val readMap = mutableMapOf<List<String>, Int>()

            //generate some mappings
            //for each range add 2 reads that maps only to A
            //add 1 read that maps only to B
            //add 2 reads that maps to A and B
            //add 1 read that maps to A and B and Ref

            //since 38 ranges have haplotypes for all samples,
            //hapA should have a count of 38 * 5 = 190
            //hapB count is B only = 38
            //hapRef only is 0
            val hapA = graph.sampleToHapId(range, sampleGameteA)
            val hapB = graph.sampleToHapId(range, sampleGameteB)
            val hapRef = graph.sampleToHapId(range, sampleGameteRef)
            if (hapA != null && hapB != null && hapRef != null) {
                readMap[listOf(hapA)] = 2
                readMap[listOf(hapB)] = 1
                readMap[listOf(hapA,hapB)] = 2
                readMap[listOf(hapA, hapB, hapRef)] = 1
            }

            readCounts[range] = readMap

        }

        return readCounts
    }
}