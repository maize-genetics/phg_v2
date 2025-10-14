package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.junit.jupiter.api.*
import java.io.File
import kotlin.test.assertEquals

class ConvertRm2Ps4gFileTest {

    companion object {
        val tempTestDir = "${TestExtension.tempDir}ropebwtTest/"

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
            File(TestExtension.testOutputFastaDir).deleteRecursively()
            File(TestExtension.testOutputDir).deleteRecursively()
            File(tempTestDir).deleteRecursively()

            File(TestExtension.tempDir).mkdirs()
            File(TestExtension.testOutputFastaDir).mkdirs()
            File(TestExtension.testOutputDir).mkdirs()
            File(tempTestDir).mkdirs()
        }
    }

    @Test
    fun testCliktParams() {
        val convertRm2Ps4gFile = ConvertRm2Ps4gFile()

        val noReadMappingFile = convertRm2Ps4gFile.test("--output-dir outputDir --hvcf-dir hvcfDir")
        assertEquals(1, noReadMappingFile.statusCode)
        Assertions.assertEquals(
            "Usage: convert-rm2ps4g-file [<options>]\n\n" +
                    "Error: missing option --read-mapping-file\n", noReadMappingFile.stderr
        )

        val noOutputDir = convertRm2Ps4gFile.test("--read-mapping-file readMappingFile --hvcf-dir hvcfDir")
        assertEquals(1, noOutputDir.statusCode)
        Assertions.assertEquals(
            "Usage: convert-rm2ps4g-file [<options>]\n\n" +
                    "Error: missing option --output-dir\n", noOutputDir.stderr
        )

        val noHVCFDir = convertRm2Ps4gFile.test("--read-mapping-file readMappingFile --output-dir outputDir")
        assertEquals(1, noHVCFDir.statusCode)
        Assertions.assertEquals(
            "Usage: convert-rm2ps4g-file [<options>]\n\n" +
                    "Error: missing option --hvcf-dir\n", noHVCFDir.stderr
        )
    }

    @Test
    fun testConvertReadMappingFile() {
        val convertRm2Ps4gFile = ConvertRm2Ps4gFile()
        val hvcfFiles = listOf("data/test/smallseq/LineA.h.vcf", "data/test/smallseq/LineB.h.vcf", "data/test/smallseq/Ref.h.vcf")
        val graph = HaplotypeGraph(hvcfFiles)
        val readMappingFile = "data/test/ropebwt/LineA_1_readMapping.txt"
        val cliCommand = "dummy Command"
        val outputDir = tempTestDir

        val readMappings = bufferedReader(readMappingFile).readLines().filter { !it.startsWith("#") }
            .filter { !it.startsWith("HapIds") }.map { it.split("\t") }.map { Pair(it[0].split(","), it[1].toInt()) }.toMap()

        convertRm2Ps4gFile.convertReadMappingFile(readMappingFile, outputDir, graph, cliCommand)

        //check the output file
        val outputFile = File("$outputDir/LineA_1_readMapping_ps4g.txt")
        val outputLines = outputFile.readLines()
        val header = outputLines.filter { it.startsWith("#") }
        val sampleGameteCountOutput = header.filter { it.startsWith("#Line") || it.startsWith("#Ref") }.map { it.split("\t") }.map { Triple(it.first(), it[1].toInt(), it[2].toInt()) }

        val hapIdToSampleGamete = graph.hapIdsToSampleGametes()

        val chrToIndexMap = graph.contigs.mapIndexed { index, s -> Pair(s, index) }.toMap()


        val gameteToIdxMap = sampleGameteCountOutput.associate { it.first.replace("#","") to it.second }
        val gameteIdxToTotalCount = sampleGameteCountOutput.associate { it.first to it.third}

        val ps4GData  = outputLines.filter { !it.startsWith("#") }.filter { !it.startsWith("gameteSet") }.map { it.split("\t") }

        assertEquals(93, ps4GData.size)

        val ps4gMap = ps4GData.associate {
            Triple(
                it[0].split(",").map { it.toInt() },
                chrToIndexMap[it[1]],
                it[2].toInt()
            ) to it[3].toInt()
        }

        for(readMapping in readMappings) {
            val expectedHapIdSet = readMapping.key

            val expectedGameteIdx = expectedHapIdSet.flatMap { hapIdToSampleGamete[it]!! }.map { gameteToIdxMap["$it"]!! }.sorted()

            val expectedData = readMapping.value
            val refRange = graph.hapIdToRefRangeMap()[readMapping.key[0]]!!.first()

            val observedDataCount = ps4gMap[Triple(expectedGameteIdx,chrToIndexMap[refRange.contig],refRange.start/256)]!!

            assertEquals(expectedData, observedDataCount)
        }
    }

    @Test
    fun testReadInReadMappingFile() {
        val convertRm2Ps4gFile = ConvertRm2Ps4gFile()
        val (header, readMappings) = convertRm2Ps4gFile.readInReadMappingFile("data/test/ropebwt/LineA_1_readMapping.txt")

        assertEquals(2, header.size)
        //#sampleName=LineA
        //#filename1=data/test/kmerReadMapping/simulatedReads/LineA_1.fq
        assertEquals("#sampleName=LineA", header[0])
        assertEquals("#filename1=data/test/kmerReadMapping/simulatedReads/LineA_1.fq", header[1])

        assertEquals(93, readMappings.size)
        // manually check a few of them
        //12f0cec9102e84a161866e37072443b7	853
        assertEquals(853, readMappings[listOf("12f0cec9102e84a161866e37072443b7")])
        //25413a93cd7622fa44d015d6914f344d,9650638ee16a4853495610120e1323f8,bb4b6391a180e13f4e2448674cba6d43	5
        assertEquals(5, readMappings[listOf("25413a93cd7622fa44d015d6914f344d", "9650638ee16a4853495610120e1323f8", "bb4b6391a180e13f4e2448674cba6d43")])
        //8f2731708b30e1c402c6d6a69a983fe4,c052e7e5a036bddc9f364324b203f1b3	24
        assertEquals(24, readMappings[listOf("8f2731708b30e1c402c6d6a69a983fe4", "c052e7e5a036bddc9f364324b203f1b3")])
    }

    @Test
    fun testConvertReadMappingDataToPS4G() {
        val convertRm2Ps4gFile = ConvertRm2Ps4gFile()
        val hvcfFiles = listOf("data/test/smallseq/LineA.h.vcf", "data/test/smallseq/LineB.h.vcf", "data/test/smallseq/Ref.h.vcf")
        val graph = HaplotypeGraph(hvcfFiles)

        val readMappings = mapOf(
            listOf("12f0cec9102e84a161866e37072443b7") to 853,
            listOf("25413a93cd7622fa44d015d6914f344d", "9650638ee16a4853495610120e1323f8", "bb4b6391a180e13f4e2448674cba6d43") to 5,
            listOf("8f2731708b30e1c402c6d6a69a983fe4", "c052e7e5a036bddc9f364324b203f1b3") to 24,
            listOf("8f2731708b30e1c402c6d6a69a983fe4") to 15
        )
        val (ps4GData, sampleGameteCount, gameteToIdxMap) = convertRm2Ps4gFile.convertReadMappingDataToPS4G(readMappings, graph)
        //check the ps4G Data
        assertEquals(4, ps4GData.size)
        //12f0cec9102e84a161866e37072443b7	853
        assertEquals(1, ps4GData[0].gameteList.size)
        assertEquals(0, ps4GData[0].gameteList[0])
        assertEquals(853, ps4GData[0].count)
        assertEquals(Position("1", 0), ps4GData[0].refPos)

        //25413a93cd7622fa44d015d6914f344d,9650638ee16a4853495610120e1323f8,bb4b6391a180e13f4e2448674cba6d43	5
        assertEquals(3, ps4GData[1].gameteList.size)
        assertEquals(0, ps4GData[1].gameteList[0])
        assertEquals(1, ps4GData[1].gameteList[1])
        assertEquals(2, ps4GData[1].gameteList[2])
        assertEquals(5, ps4GData[1].count)
        //34001 /256 = 132
        assertEquals(Position("1", 132), ps4GData[1].refPos)

        //8f2731708b30e1c402c6d6a69a983fe4,c052e7e5a036bddc9f364324b203f1b3	24
        assertEquals(2, ps4GData[2].gameteList.size)
        //this is swapped 8f2731 is from LineB and c052e7 is from LineA.  This is correct as we sort by id and LineA comes first
        assertEquals(0, ps4GData[2].gameteList[0])
        assertEquals(1, ps4GData[2].gameteList[1])
        assertEquals(24, ps4GData[2].count)
        //28501 /256 = 111
        assertEquals(Position("2", 111), ps4GData[2].refPos)

        //8f2731708b30e1c402c6d6a69a983fe4	15
        assertEquals(1, ps4GData[3].gameteList.size)
        assertEquals(1, ps4GData[3].gameteList[0])
        assertEquals(15, ps4GData[3].count)
        //28501 /256 = 111
        assertEquals(Position("2", 111), ps4GData[3].refPos)


        //check the sampleGameteCount
        assertEquals(3, sampleGameteCount.size)
        assertEquals(882, sampleGameteCount[SampleGamete("LineA")])
        assertEquals(44, sampleGameteCount[SampleGamete("LineB")])
        assertEquals(5, sampleGameteCount[SampleGamete("Ref")])

        //check the gameteToIdxMap
        assertEquals(3, gameteToIdxMap.size)
        assertEquals(0, gameteToIdxMap[SampleGamete("LineA")])
        assertEquals(1, gameteToIdxMap[SampleGamete("LineB")])
        assertEquals(2, gameteToIdxMap[SampleGamete("Ref")])


    }

    @Test
    fun testCreatePS4GFileForSingleMapping() {
        val convertRm2Ps4gFile = ConvertRm2Ps4gFile()
        val hvcfFiles = listOf("data/test/smallseq/LineA.h.vcf", "data/test/smallseq/LineB.h.vcf", "data/test/smallseq/Ref.h.vcf")

        val graph = HaplotypeGraph(hvcfFiles)
        val hapIdToRanges = graph.hapIdToRefRangeMap()
        val contigToIdxMap = graph.contigs.mapIndexed { index, s -> Pair(s, index) }.toMap()
        val gameteCountMap = mutableMapOf<SampleGamete, Int>()
        val gameteToIdxMap = graph.sampleGametesInGraph().mapIndexed { index, s -> Pair(s, index) }.toMap()

        //12f0cec9102e84a161866e37072443b7	853
        val hapIdSet = listOf("12f0cec9102e84a161866e37072443b7")
        val ps4GDataFirst = convertRm2Ps4gFile.createPS4GFileForSingleMapping(
            hapIdSet, hapIdToRanges, contigToIdxMap, graph, gameteCountMap, 853, gameteToIdxMap
        )
        assertEquals(1, ps4GDataFirst.gameteList.size)
        assertEquals(0, ps4GDataFirst.gameteList[0])
        assertEquals(853, ps4GDataFirst.count)
        assertEquals(Position("1", 0), ps4GDataFirst.refPos)


        //3149b3144f93134eb29661bade697fc6,8967fabf10e55d881caa6fe192e7d4ca	40
        val hapIdSet2 = listOf("3149b3144f93134eb29661bade697fc6", "8967fabf10e55d881caa6fe192e7d4ca")
        val ps4GDataSecond = convertRm2Ps4gFile.createPS4GFileForSingleMapping(
            hapIdSet2, hapIdToRanges, contigToIdxMap, graph, gameteCountMap, 40, gameteToIdxMap
        )
        assertEquals(2, ps4GDataSecond.gameteList.size)
        assertEquals(0, ps4GDataSecond.gameteList[0])
        assertEquals(1, ps4GDataSecond.gameteList[1])
        assertEquals(40, ps4GDataSecond.count)
        //1001/256 = 3.9 -> 3
        assertEquals(Position("1", 3), ps4GDataSecond.refPos)

        //5be0e52ccfe573a6e42e4dd1e8658105,a7214fe07512b511ddf13edade461b39,fafafee7c250c76b7e0571fde286022e	33
        val hapIdSet3 = listOf("5be0e52ccfe573a6e42e4dd1e8658105", "a7214fe07512b511ddf13edade461b39", "fafafee7c250c76b7e0571fde286022e")
        val ps4GDataThird = convertRm2Ps4gFile.createPS4GFileForSingleMapping(
            hapIdSet3, hapIdToRanges, contigToIdxMap, graph, gameteCountMap, 33, gameteToIdxMap
        )
        assertEquals(3, ps4GDataThird.gameteList.size)
        assertEquals(0, ps4GDataThird.gameteList[0])
        assertEquals(1, ps4GDataThird.gameteList[1])
        assertEquals(2, ps4GDataThird.gameteList[2])
        assertEquals(33, ps4GDataThird.count)
        //49501/256 = 193
        assertEquals(Position("1", 193), ps4GDataThird.refPos)
    }
}