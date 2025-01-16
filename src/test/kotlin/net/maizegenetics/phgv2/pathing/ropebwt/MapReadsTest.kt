package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.pathing.KeyFileData
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import java.io.File
import kotlin.test.assertEquals
import kotlin.test.assertTrue
import kotlin.test.fail

class MapReadsTest {

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
        val mapReads = MapReads()

        val noIndex = mapReads.test("--read-files test1.fq --output-dir testDir")
        assertEquals(1, noIndex.statusCode)
        assertEquals("Usage: map-reads [<options>]\n\n" +
                "Error: missing option --index\n", noIndex.stderr)

        val noReads = mapReads.test("--index testIndex --output-dir testDir")
        assertEquals(1, noReads.statusCode)
        assertEquals("Usage: map-reads [<options>]\n\n" +
                "Error: must provide one of --key-file, --read-files\n", noReads.stderr)

        val bothReadInputs = mapReads.test("--index testIndex --key-file testKeyFile --read-files test1.fq --output-dir testDir")
        assertEquals(1, bothReadInputs.statusCode)
        assertEquals("Usage: map-reads [<options>]\n\n" +
                "Error: option --key-file cannot be used with --read-files\n", bothReadInputs.stderr)

        val noOutputDir = mapReads.test("--index testIndex --read-files test1.fq")
        assertEquals(1, noOutputDir.statusCode)
        assertEquals("Usage: map-reads [<options>]\n\n" +
                "Error: missing option --output-dir\n", noOutputDir.stderr)
    }

    @Test
    fun testParseMem() {
        val mapReads = MapReads()
        val alignmentString = "ST-E00317:129:HVMFTCCXX:7:1101:5944:1309\t1\t150\t22\t22\t" +
                "598e490a361fdaf13f1f7e4fd6aeeac4:+:22036\t2bb15f8e74f01c33164895f6aa0c1955:+:25669\te0eec1384d9dfa4e62fcca27e4caadaa:-:35325\t" +
                "6d6966023182c5a182d1a0b650e4f620:+:20935\td70b3c72ddf6d22311fe9a97a2fa8cb7:-:36296\td6c5250d88e7e71eecfee33ec0b8f2a3:-:36297\t" +
                "04047f19cc55dae88853913773854067:-:31242\t0f9ecde512ae2b799233f8acc0e7e43c:-:16635\t455dc7b6a868a2e6c14133330228d05c:+:70219\t" +
                "af469871ace61261d0c0ecb257129d53:+:20934\te6c8ed95a7bbd616cd7c73ad68b1885b:+:22105\tdc54493d608e10bbec3344cd34aa14fd:+:155693\t" +
                "0a38c98a1af292f078f3879a7d0d8179:+:70203\ta0db88f214f86d47c0fbecd8cbd16143:+:22163\t020e8351832417a8522fd2437e475ce1:-:8308\t" +
                "d84bea2a0e8bcf79fd1ba8781ad0818c:+:139499\t0bfa126c7d600678dd1dc9a66f80119b:-:35325\t8ad9073572cf3e16105f321a2727f6b7:-:35324\t" +
                "171bd10310c5affc6313856d7b017841:-:19667\tf304773b801a9249bf36bc4b6cafb9b4:-:34309\t1e77759e3c4f8437c2d6a10de91fbf9f:-:16641\t" +
                "233e06c196319b1b02c6210c0748dd62:+:5021"
        val parsed = mapReads.parseStringIntoMem(alignmentString)
        assertEquals("ST-E00317:129:HVMFTCCXX:7:1101:5944:1309", parsed.readName)
        assertEquals(1, parsed.readStart)
        assertEquals(150, parsed.readEnd)
        assertEquals(22, parsed.numHits)
        assertEquals(22, parsed.listMemHits.size)
        assertEquals("598e490a361fdaf13f1f7e4fd6aeeac4", parsed.listMemHits[0].contig)
        assertEquals("+", parsed.listMemHits[0].strand)
        assertEquals(22036, parsed.listMemHits[0].pos)
        assertEquals("2bb15f8e74f01c33164895f6aa0c1955", parsed.listMemHits[1].contig)
        assertEquals("+", parsed.listMemHits[1].strand)
        assertEquals(25669, parsed.listMemHits[1].pos)
        assertEquals("e0eec1384d9dfa4e62fcca27e4caadaa", parsed.listMemHits[2].contig)
        assertEquals("-", parsed.listMemHits[2].strand)
        assertEquals(35325, parsed.listMemHits[2].pos)

    }

    @Test
    fun testProcessMemsForRead() {
        val mapReads = MapReads()

        val hapIdToRefRangeMap = mapOf("hap1" to listOf(ReferenceRange("chr1",100,200)), "hap2" to listOf(ReferenceRange("chr1",100,200)), "hap3" to listOf(ReferenceRange("chr1",100,200)),
            "hap4" to listOf(ReferenceRange("chr1",100,200)), "hap5" to listOf(ReferenceRange("chr1",100,200)))

        //Make a simple 1 MEM hit
        val simpleMEMList = listOf(MEM("read1",0,150,3,listOf(MEMHit("hap1", "+", 100), MEMHit("hap2", "-", 200), MEMHit("hap3", "+", 300))))

        var simpleReadMapping = mutableMapOf<List<String>, Int>()
        mapReads.processMemsForRead(simpleMEMList,simpleReadMapping,5, hapIdToRefRangeMap)
        assertEquals(1, simpleReadMapping.size)
        assertEquals(listOf("hap1","hap2","hap3"), simpleReadMapping.keys.first())
        assertEquals(1, simpleReadMapping[listOf("hap1","hap2","hap3")])


        simpleReadMapping = mutableMapOf()
        mapReads.processMemsForRead(simpleMEMList,simpleReadMapping,2, hapIdToRefRangeMap)
        assertEquals(0, simpleReadMapping.size)


        //Create a few reads with the same hapIds hit
        val simpleMEMList2 = listOf(MEM("read2",0,150,3,listOf(MEMHit("hap1", "+", 100), MEMHit("hap2", "-", 200), MEMHit("hap3", "+", 300))))
        val simpleMEMList3 = listOf(
            MEM("read3",2,150,3,listOf(MEMHit("hap1", "+", 100), MEMHit("hap2", "-", 200), MEMHit("hap3", "+", 300))),
            MEM("read3",1,150,3,listOf(MEMHit("hap4", "+", 100), MEMHit("hap5", "-", 200)))
        )

        val simpleMEMList4 = listOf(
            MEM("read3",2,150,3,listOf(MEMHit("hap1", "+", 100), MEMHit("hap2", "-", 200), MEMHit("hap3", "+", 300))),
            MEM("read3",2,150,3,listOf(MEMHit("hap4", "+", 100), MEMHit("hap5", "-", 200)))
        )

        simpleReadMapping = mutableMapOf()
        mapReads.processMemsForRead(simpleMEMList,simpleReadMapping,6, hapIdToRefRangeMap)
        mapReads.processMemsForRead(simpleMEMList2,simpleReadMapping,6, hapIdToRefRangeMap)
        mapReads.processMemsForRead(simpleMEMList3,simpleReadMapping,6, hapIdToRefRangeMap)
        mapReads.processMemsForRead(simpleMEMList4,simpleReadMapping,6, hapIdToRefRangeMap)
        assertEquals(3, simpleReadMapping.size)
        assertTrue(simpleReadMapping.keys.contains(listOf("hap1","hap2","hap3")))
        assertEquals(2, simpleReadMapping[listOf("hap1","hap2","hap3")])
        assertTrue(simpleReadMapping.keys.contains(listOf("hap4","hap5")))
        assertEquals(1, simpleReadMapping[listOf("hap4","hap5")])
        assertTrue(simpleReadMapping.keys.contains(listOf("hap1","hap2","hap3","hap4","hap5")))
        assertEquals(1, simpleReadMapping[listOf("hap1","hap2","hap3","hap4","hap5")])


        //pass in empty list
        simpleReadMapping = mutableMapOf()
        assertThrows<NoSuchElementException> { mapReads.processMemsForRead(listOf(),simpleReadMapping,6, hapIdToRefRangeMap) }

    }

    @Test
    fun testCreateReadMappingsForFileReader() {
        val reader = bufferedReader("data/test/ropebwt/alignment.bed")

        val hapIdToRefRangeMap = mapOf("hap1" to listOf(ReferenceRange("chr1",100,200)),
            "hap2" to listOf(ReferenceRange("chr1",100,200)),
            "hap3" to listOf(ReferenceRange("chr1",100,200)),
            "hap4" to listOf(ReferenceRange("chr1",100,200)),
            "hap5" to listOf(ReferenceRange("chr1",100,200)))


        val mapReads = MapReads()
        val readMapping = mapReads.createReadMappingsForFileReader(reader, 5, hapIdToRefRangeMap)

        assertEquals(2, readMapping.size)
        assertTrue(readMapping.keys.contains(listOf("hap1","hap2","hap3")))
        assertEquals(2, readMapping[listOf("hap1","hap2","hap3")])
        assertTrue(readMapping.keys.contains(listOf("hap4","hap5")))
        assertEquals(1, readMapping[listOf("hap4","hap5")])

    }

    @Test
    fun testSetupMappingProcess() {
        val expected = bufferedReader("data/test/ropebwt/LineA_1_expected.bed").readLines()

        val mapReads = MapReads()
        val index = "data/test/ropebwt/testIndex.fmd"
        val readFile = "data/test/kmerReadMapping/simulatedReads/LineA_1.fq"
        //index: String, readFile: String, threads: Int, minMemLength: Int, maxNumHits: Int, condaEnvPrefix: String

        val observedAlignments = mapReads.setupMappingProcess(index, readFile, 5, 148, 5, "").readLines().filter { it.isNotEmpty() }

        assertEquals(expected.size, observedAlignments.size)
        expected.indices.forEach { assertEquals(expected[it], observedAlignments[it]) }
    }

    @Test
    fun testMapSingleReadFile() {
        resetDirs()
        val mapReads = MapReads()
        val index = "data/test/ropebwt/testIndex.fmd"
        val readFile = "data/test/kmerReadMapping/simulatedReads/LineA_1.fq"
        val outputDir = "${TestExtension.tempDir}ropebwtTest/"
        val outputFile = "$outputDir/LineA_1_readMapping.txt"


        val hvcfFiles = File(TestExtension.smallSeqInputDir).walkTopDown().filter { it.isFile }
            .filter { it.name.endsWith("h.vcf") || it.name.endsWith("h.vcf.gz") }.map { "${TestExtension.smallSeqInputDir}/${it.name}" }
            .toList()

        val graph = HaplotypeGraph(hvcfFiles)

        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()

        mapReads.mapSingleReadFile(index, "LineA" ,readFile, outputFile ,5, 148, 5, "", hapIdToRefRangeMap)

        val expectedLines = bufferedReader("data/test/ropebwt/LineA_1_readMapping.txt").readLines()
        val expectedReadMap = expectedLines.filter { !it.startsWith("#") }.filter { !it.startsWith("HapId") }.map { it.split("\t") }.associate { it[0] to it[1].toInt() }
        val observedLines = bufferedReader(outputFile).readLines()
        val observedReadMap = observedLines.filter { !it.startsWith("#") }.filter { !it.startsWith("HapId") }.map { it.split("\t") }.associate { it[0] to it[1].toInt() }

        assertEquals(expectedReadMap.size, observedReadMap.size)

        expectedReadMap.keys.forEach { assertEquals(expectedReadMap[it], observedReadMap[it]) }
    }

    @Test
    fun testMapAllReadFiles() {
        resetDirs()
        val mapReads = MapReads()
        //process all the reads
        val outputDir = "${TestExtension.tempDir}ropebwtTest/"
        val index = "data/test/ropebwt/testIndex.fmd"

        val hvcfFiles = File(TestExtension.smallSeqInputDir).walkTopDown().filter { it.isFile }
            .filter { it.name.endsWith("h.vcf") || it.name.endsWith("h.vcf.gz") }.map { "${TestExtension.smallSeqInputDir}/${it.name}" }
            .toList()

        val graph = HaplotypeGraph(hvcfFiles)

        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()

        val keyFileDataEntry = listOf(KeyFileData("LineA", "data/test/kmerReadMapping/simulatedReads/LineA_1.fq", "data/test/kmerReadMapping/simulatedReads/LineA_2.fq"))
        mapReads.mapAllReadFiles(index, keyFileDataEntry, outputDir, 5, 148, 5, "", hapIdToRefRangeMap )


        val expectedLines = bufferedReader("data/test/ropebwt/LineA_1_readMapping.txt").readLines()
        val expectedReadMap = expectedLines.filter { !it.startsWith("#") }.filter { !it.startsWith("HapId") }.map { it.split("\t") }.associate { it[0] to it[1].toInt() }
        val observedLines = bufferedReader("$outputDir/LineA_1_readMapping.txt").readLines()
        val observedReadMap = observedLines.filter { !it.startsWith("#") }.filter { !it.startsWith("HapId") }.map { it.split("\t") }.associate { it[0] to it[1].toInt() }

        assertEquals(expectedReadMap.size, observedReadMap.size)
        expectedReadMap.keys.forEach { assertEquals(expectedReadMap[it], observedReadMap[it]) }

        val expectedLines2 = bufferedReader("data/test/ropebwt/LineA_2_readMapping.txt").readLines()
        val expectedReadMap2 = expectedLines2.filter { !it.startsWith("#") }.filter { !it.startsWith("HapId") }.map { it.split("\t") }.associate { it[0] to it[1].toInt() }
        val observedLines2 = bufferedReader("$outputDir/LineA_2_readMapping.txt").readLines()
        val observedReadMap2 = observedLines2.filter { !it.startsWith("#") }.filter { !it.startsWith("HapId") }.map { it.split("\t") }.associate { it[0] to it[1].toInt() }

        assertEquals(expectedReadMap2.size, observedReadMap2.size)
        expectedReadMap2.keys.forEach { assertEquals(expectedReadMap2[it], observedReadMap2[it]) }
    }

    @Test
    fun testExportPathKeyFile() {
        val mapReads = MapReads()
        val outputDir = "${TestExtension.tempDir}ropebwtTest/"
        val readNameToFileMap = mutableMapOf<String,MutableList<String>>()
        readNameToFileMap["LineA"] = mutableListOf("$outputDir/LineA_1_readMapping.txt", "$outputDir/LineA_2_readMapping.txt")
        readNameToFileMap["LineB"] = mutableListOf("$outputDir/LineB_1_readMapping.txt")

        mapReads.exportPathKeyFile(outputDir, readNameToFileMap)

        val expectedLines = listOf("sampleName\tfilename", "LineA\t$outputDir/LineA_1_readMapping.txt", "LineA\t$outputDir/LineA_2_readMapping.txt", "LineB\t$outputDir/LineB_1_readMapping.txt")
        val observedLines = bufferedReader("$outputDir/pathKeyFile.txt").readLines()

        assertEquals(expectedLines.size, observedLines.size)
        expectedLines.indices.forEach { assertEquals(expectedLines[it], observedLines[it]) }
    }
}