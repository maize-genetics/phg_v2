package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import java.nio.file.Files
import kotlin.test.assertEquals
import kotlin.test.fail

class ReadMappingCountQcTest {

    companion object {
        val tempTestDir = TestExtension.testOutputDir
        val hvcfDir = "${tempTestDir}hvcfDir/"

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
            File(hvcfDir).deleteRecursively()

            File(TestExtension.tempDir).mkdirs()
            File(TestExtension.testOutputFastaDir).mkdirs()
            File(TestExtension.testOutputDir).mkdirs()
            File(tempTestDir).mkdirs()
            File(hvcfDir).mkdirs()
        }
    }

    @Test
    fun testCliktParams() {
        val readMappingCountQc = ReadMappingCountQc()

        val noHvcfDir = readMappingCountQc.test("--read-mapping-file dummyReadMapping.txt --target-sample-name dummySampleName --output-dir dummyOutputDir")
        assertEquals(1,noHvcfDir.statusCode)
        assertEquals(
            "Usage: read-mapping-count-qc [<options>]\n\n" +
                    "Error: missing option --hvcf-dir\n", noHvcfDir.stderr
        )

        val noReadMappingFile = readMappingCountQc.test("--hvcf-dir dummyHvcfDir --target-sample-name dummySampleName --output-dir dummyOutputDir")
        assertEquals(1,noReadMappingFile.statusCode)
        assertEquals(
            "Usage: read-mapping-count-qc [<options>]\n\n" +
                    "Error: missing option --read-mapping-file\n", noReadMappingFile.stderr
        )

        val noTargetSampleName = readMappingCountQc.test("--hvcf-dir dummyHvcfDir --read-mapping-file dummyReadMapping.txt --output-dir dummyOutputDir")
        assertEquals(1,noTargetSampleName.statusCode)
        assertEquals(
            "Usage: read-mapping-count-qc [<options>]\n\n" +
                    "Error: missing option --target-sample-name\n", noTargetSampleName.stderr
        )

        val noOutputDir = readMappingCountQc.test("--hvcf-dir dummyHvcfDir --read-mapping-file dummyReadMapping.txt --target-sample-name dummySampleName")
        assertEquals(1,noOutputDir.statusCode)
        assertEquals(
            "Usage: read-mapping-count-qc [<options>]\n\n" +
                    "Error: missing option --output-dir\n", noOutputDir.stderr
        )

    }

    //processReadMappingCounts(hvcfDir: String, readMappingFile: String, targetSampleName: String, outputDir: String)
    @Test
    fun testProcessReadMappingCounts() {
        val readMappingCountQc = ReadMappingCountQc()
        //Move the hvcfs into the temp directory

        Files.copy(File(TestExtension.smallseqLineAHvcfFile).toPath(), File("${hvcfDir}LineA.hvcf").toPath())
        Files.copy(File(TestExtension.smallseqLineBHvcfFile).toPath(), File("${hvcfDir}LineB.hvcf").toPath())
        Files.copy(File(TestExtension.smallseqRefHvcfFile).toPath(), File("${hvcfDir}Ref.hvcf").toPath())



        val readMappingFile = "data/test/ropebwt/LineA_1_readMapping.txt"
        val targetSampleName = "LineA"
        val outputDir = tempTestDir
        readMappingCountQc.processReadMappingCounts(hvcfDir, readMappingFile, targetSampleName, outputDir)

        val expectedFile = "data/test/kmerReadMapping/LineA_1_readMapping_LineA_counts.txt"
        val actualFile = "$outputDir/LineA_1_readMapping_LineA_counts.txt"
        val expected = bufferedReader(expectedFile).readLines().joinToString("\n")
        val actual = bufferedReader(actualFile).readLines().joinToString("\n")
        assertEquals(expected, actual)
    }

    @Test
    fun testConvertReadMappingToHapIdCounts() {
        val readMappingCountQc = ReadMappingCountQc()
        val readMappings = mapOf(listOf("hap1","hap2") to 1, listOf("hap2") to 2,
            listOf("hap1") to 3, listOf("hap3") to 4, listOf("hap4","hap5") to 6)
        val expected = mapOf("hap1" to 4, "hap2" to 3, "hap3" to 4, "hap4" to 6, "hap5" to 6)
        val actual = readMappingCountQc.convertReadMappingToHapIdCounts(readMappings)
        assertEquals(expected, actual)
    }

    @Test
    fun testBuildCountOutputFile() {
        //buildCountOutputFile(outputDir: String, readMappingFile: String, targetSampleName: String) : String
        val readMappingCountQc = ReadMappingCountQc()
        val outputDir = "dummyOutputDir"
        val readMappingFile = "dummyReadMapping.txt"
        val targetSampleName = "dummySampleName"
        val expected = "$outputDir/dummyReadMapping_dummySampleName_counts.txt"
        val actual = readMappingCountQc.buildCountOutputFile(outputDir, readMappingFile, targetSampleName)
        assertEquals(expected, actual)
    }

    /**
     * fun writeOutCounts(hapIdCounts: Map<String,Int>,
     *                        outputFile: String,
     *                        targetSampleName: String,
     *                        rangeToHapId: Map<ReferenceRange, List<String>>,
     *                        hapIdToSampleGamete: Map<String, List<SampleGamete>>) {
     */
    @Test
    fun testWriteOutCounts() {
        val readMappingCountQc = ReadMappingCountQc()
        val hapIdCounts = mapOf("hap1" to 4, "hap2" to 3, "hap3" to 3, "hap4" to 6, "hap5" to 6, "hap6" to 3)
        val outputFile = "${tempTestDir}/simpleOutput.txt"
        val targetSampleName = "sample1"
        val rangeToHapIds = mapOf(ReferenceRange("1",10,50) to listOf("hap1","hap2"),
            ReferenceRange("1",51,75) to listOf("hap3","hap4"),
            ReferenceRange("1",76,100) to listOf("hap5","hap6"))

        val hapIdToSampleGamete = mapOf(
            "hap1" to listOf(SampleGamete("sample1", 0), SampleGamete("sample2", 0)),
            "hap2" to listOf(SampleGamete("sample3", 0)),
            "hap3" to listOf(SampleGamete("sample1", 0), SampleGamete("sample3", 0)),
            "hap4" to listOf(SampleGamete("sample2", 0)),
            "hap5" to listOf(SampleGamete("sample1", 0)),
            "hap6" to listOf(SampleGamete("sample2", 0), SampleGamete("sample3", 0)))

        readMappingCountQc.writeOutCounts(hapIdCounts, outputFile, targetSampleName, rangeToHapIds, hapIdToSampleGamete)

        val expected = "refRange\tsample1_HapID\tsample1_HapCount\tHighestAltCount\tDifference\tOtherHapCounts\n" +
            "1:10-50\thap1\t4\t3\t1\t3_hap2\n" +
            "1:51-75\thap3\t3\t6\t-3\t6_hap4\n" +
            "1:76-100\thap5\t6\t3\t3\t3_hap6"

        val actual = bufferedReader(outputFile).readLines().joinToString("\n")
        assertEquals(expected, actual)


    }

    @Test
    fun testBuildOutputStringForHapIdsInRefRange() {
        val readMappingCountQc = ReadMappingCountQc()
        val hapIds = listOf("hap1","hap2","hap3")
        val hapIdCounts = mapOf("hap1" to 4, "hap2" to 3, "hap3" to 3, "hap4" to 6, "hap5" to 6)
        val hapIdToSampleGamete = mapOf("hap1" to listOf(SampleGamete("sample1", 0), SampleGamete("sample2", 0)),
            "hap2" to listOf(SampleGamete("sample3", 0)),
            "hap3" to listOf(SampleGamete("sample4", 0), SampleGamete("sample5", 0)))
        val rangeToHapIds = mapOf(ReferenceRange("1",10,50) to listOf("hap1","hap2"), ReferenceRange("1",51,75) to listOf("hap2","hap3"), ReferenceRange("1",76,100) to listOf("hap1","hap3"))
        val expected = "1:10-50\thap1\t4\t3\t1\t3_hap2\n"
        val actual = readMappingCountQc.buildOutputStringForHapIdsInRefRange(rangeToHapIds, ReferenceRange("1",10,50), hapIdToSampleGamete , hapIdCounts,"sample1")
        assertEquals(expected, actual)

        //test that the non Target is higher
        val expectedNonHigher = "1:10-50\thap2\t3\t4\t-1\t4_hap1\n"
        val actualNonHigher = readMappingCountQc.buildOutputStringForHapIdsInRefRange(rangeToHapIds, ReferenceRange("1",10,50), hapIdToSampleGamete , hapIdCounts,"sample3")
        assertEquals(expectedNonHigher, actualNonHigher)

        //Test what happens if the target sample is not in the current refRange
        val expectedNoTarget = "1:10-50\t\t0\t4\t-4\t4_hap1, 3_hap2\n"
        val actualNoTarget = readMappingCountQc.buildOutputStringForHapIdsInRefRange(rangeToHapIds, ReferenceRange("1",10,50), hapIdToSampleGamete , hapIdCounts,"sample6")
        assertEquals(expectedNoTarget, actualNoTarget)
    }

    @Test
    fun testComputeCountsForNonTargetHapids() {
        val readMappingCountQc = ReadMappingCountQc()
        //fun computeCountsForHapids(
        //        hapIds: List<String>,
        //        targetHapId: String,
        //        hapIdCounts: Map<String, Int>
        //    ) :List<Pair<Int,String>> {
        val hapIds = listOf("hap1","hap2","hap3")
        val targetHapId = "hap2"
        val hapIdCounts = mapOf("hap1" to 4, "hap2" to 3, "hap3" to 3, "hap4" to 6, "hap5" to 6)
        val expected = listOf(Pair(4,"hap1"), Pair(3,"hap3"))
        val actual = readMappingCountQc.computeCountsForNonTargetHapids(hapIds, targetHapId, hapIdCounts)
        assertEquals(expected, actual)

        val hapIdsMixed = listOf("hap3", "hap2", "hap1")
        val actualMixed = readMappingCountQc.computeCountsForNonTargetHapids(hapIdsMixed, targetHapId, hapIdCounts)
        assertEquals(expected, actualMixed)

        val hapIdsDup = listOf("hap1","hap2","hap3","hap1","hap2","hap3")
        val actualDup = readMappingCountQc.computeCountsForNonTargetHapids(hapIdsDup, targetHapId, hapIdCounts)
        assertEquals(expected, actualDup)

        val hapIdsNoTarget = listOf("hap1","hap3")
        val actualNoTarget = readMappingCountQc.computeCountsForNonTargetHapids(hapIdsNoTarget, targetHapId, hapIdCounts)
        assertEquals(listOf(Pair(4,"hap1"), Pair(3,"hap3")), actualNoTarget)
    }

    @Test
    fun testFindTargetHapIdForSample() {
        val readMappingCountQc = ReadMappingCountQc()

        val hapIds = listOf("hap1","hap2","hap3")
        val hapIdToSampleGamete = mapOf("hap1" to listOf(SampleGamete("sample1", 0), SampleGamete("sample2", 0)),
            "hap2" to listOf(SampleGamete("sample3", 0)),
            "hap3" to listOf(SampleGamete("sample4", 0), SampleGamete("sample5", 0)))
        val targetSampleName = "sample3"
        val expected = "hap2"
        val actual = readMappingCountQc.findTargetHapIdForSample(hapIds, hapIdToSampleGamete, targetSampleName)
        assertEquals(expected, actual)

        val missingTargetName = "sample6"
        val missingTargetExpected = ""
        val missingTargetActual = readMappingCountQc.findTargetHapIdForSample(hapIds, hapIdToSampleGamete, missingTargetName)
        assertEquals(missingTargetExpected, missingTargetActual)

        val overlapHapIdToSampleGameteMap = mapOf("hap1" to listOf(SampleGamete("sample1", 0), SampleGamete("sample2", 0)),
            "hap2" to listOf(SampleGamete("sample3", 0), SampleGamete("sample4", 0)),
            "hap3" to listOf(SampleGamete("sample3", 0), SampleGamete("sample5", 0)))

        val overlapExpected = "hap2"
        val overlapActual = readMappingCountQc.findTargetHapIdForSample(hapIds, overlapHapIdToSampleGameteMap, targetSampleName)
        assertEquals(overlapExpected, overlapActual)
    }
}