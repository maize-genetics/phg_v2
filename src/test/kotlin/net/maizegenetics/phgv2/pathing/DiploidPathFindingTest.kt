package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.testing.test
import htsjdk.variant.vcf.VCFFileReader
import junit.framework.TestCase.assertEquals
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.getBufferedWriter
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import java.nio.file.Paths
import kotlin.random.Random

@ExtendWith(TestExtension::class)
class DiploidPathFindingTest {
    companion object {
        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.testOutputDir).mkdirs()

            //erase any hvcf files in testVCFDir (which is created by TestExtension, if it does not already exist)
            val vcfDir = File(TestExtension.testVCFDir)
            val vcfFiles = vcfDir.listFiles { file -> file.name.endsWith(".h.vcf") }
            for (vcfFile in vcfFiles) {
                vcfFile.delete()
            }

            //copy hvcf files to vcfDir
            listOf(
                TestExtension.smallseqLineAHvcfFile,
                TestExtension.smallseqLineBHvcfFile,
                TestExtension.smallseqRefHvcfFile
            )
                .forEach { filepath ->
                    val origFile = File(filepath)
                    origFile.copyTo(vcfDir.resolve(origFile.name))
                }

        }

        @JvmStatic
        @AfterAll
        fun tearDown() {
            //comment out the following line to inspect the test results after the tests have been run
//            File(TestExtension.testOutputDir).deleteRecursively()
        }
    }

    @Test
    fun testCliktParameters() {
        //copied from HaploidPathFindingTest, then added validation check for inbreeding-coefficient
        //test that the parameters are being parsed correctly
        //test default parameters in HaploidPathFinding
        var result = HaploidPathFinding().test("")
        //expected errors:
        //Error: missing option --path-keyfile
        //Error: missing option --hvcf-dir
        //Error: missing option --reference-genome
        //Error: missing option --output-dir
        var errOut = result.stderr
        assert(errOut.contains("Error: missing option --path-keyfile"))
        assert(errOut.contains("Error: missing option --hvcf-dir"))
        assert(errOut.contains("Error: missing option --reference-genome"))
        assert(errOut.contains("Error: missing option --output-dir"))

        //test required parameter validation
        var testArgs = "--path-keyfile notafile --hvcf-dir notafile --reference-genome notafile --output-dir notafile"
        result = DiploidPathFinding().test(testArgs)
        errOut = result.stderr
        //expected errors:
        // Error: invalid value for --path-keyfile: notafile is not a valid file
        //Error: invalid value for --hvcf-dir: notafile is not a valid directory.
        //Error: invalid value for --reference-genome: notafile is not a valid file
        //Error: invalid value for --output-dir: notafile is not a valid directory.
        assert(errOut.contains("Error: invalid value for --path-keyfile: notafile is not a valid file"))
        assert(errOut.contains("Error: invalid value for --hvcf-dir: notafile is not a valid directory."))
        assert(errOut.contains("Error: invalid value for --reference-genome: notafile is not a valid file"))
        assert(errOut.contains("Error: invalid value for --output-dir: notafile is not a valid directory."))

        //test validation of optional parameters
        testArgs = "--path-keyfile ${TestExtension.testKeyFile} --hvcf-dir ${TestExtension.smallSeqInputDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.smallSeqInputDir} --prob-correct 0.4 " +
                "--min-gametes -1 --min-reads -1 --max-reads-per-kb -1 --min-coverage 0.4 --inbreeding-coefficient -1"

        result = DiploidPathFinding().test(testArgs)

        //Expected errors:
//        Error: invalid value for --prob-correct: prob-correct must be between 0.5 and 1.0
//        Error: invalid value for --min-gametes: min-gametes must be a positive integer
//        Error: invalid value for --min-reads: min-reads must be a positive integer.
//        Error: invalid value for --inbreeding-coefficient: inbreeding-coefficient must be between 0.0 and 1.0
//        Error: invalid value for --max-reads-per-kb: max-reads-per-kb must be a positive integer.
//        Error: invalid value for --min-coverage: min-coverage must be between 0.5 and 1.0

        errOut = result.stderr
        assert(errOut.contains("Error: invalid value for --prob-correct: prob-correct must be between 0.5 and 1.0"))
        assert(errOut.contains("Error: invalid value for --min-gametes: min-gametes must be a positive integer"))
        assert(errOut.contains("Error: invalid value for --min-reads: min-reads must be a positive integer."))
        assert(errOut.contains("Error: invalid value for --max-reads-per-kb: max-reads-per-kb must be a positive integer."))
        assert(errOut.contains("Error: invalid value for --min-coverage: min-coverage must be between 0.5 and 1.0"))
        assert(errOut.contains("Error: invalid value for --inbreeding-coefficient: inbreeding-coefficient must be between 0.0 and 1.0"))

    }

    @Test
    fun testDiploidPathFinding() {
        //plan for this test:
        //use a haplotype groph built from Ref, lineA, and lineB
        val vcfDir = File(TestExtension.testVCFDir)
        val listOfHvcfFilenames = vcfDir.listFiles().map { it.path }.filter { it.endsWith(".h.vcf") }
        val myGraph = HaplotypeGraph(listOfHvcfFilenames)

        //create a read mapping file
        val readMappingFile = TestExtension.testOutputDir + "testReadMapping.txt"
        createReadMappings(myGraph, readMappingFile)

        //create a keyfile
        val keyFilename = TestExtension.testOutputDir + "keyfileForPathTest.txt"
        getBufferedWriter(keyFilename).use { myWriter ->
            myWriter.write("SampleName\tReadMappingFiles\n")
            myWriter.write("TestLine\t$readMappingFile\n")
        }

        var pathFindingTestArgs = "--path-keyfile $keyFilename --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir} " +
                "--prob-same-gamete 0.8"

        val pathFindingResult = DiploidPathFinding().test(pathFindingTestArgs)

        //examine the resulting hvcf file
        val testVcf = "${TestExtension.testOutputDir}TestLine.h.vcf"
        val lineAgamete = SampleGamete("LineA")
        val lineBgamete = SampleGamete("LineB")
        val lineAHapids = myGraph.ranges().map { range -> myGraph.sampleToHapId(range, lineAgamete) ?: "." }.toSet()
        val lineBHapids = myGraph.ranges().map { range -> myGraph.sampleToHapId(range, lineBgamete) ?: "." }.toSet()

        VCFFileReader(Paths.get(testVcf), false).use { vcf ->
            //Chr 1: one haplotype should be all A for chr1, the other A before start=25000 and B after
            for (context in vcf) {
                val hapids = context.alleles.map { it.displayString.substringBefore(">").substringAfter("<") }
                when (context.contig) {
                    "1" -> when {
                        context.start < 25000 -> assert(hapids.all { lineAHapids.contains(it) })
                        else -> assert(hapids.any { lineAHapids.contains(it) } && hapids.any { lineBHapids.contains(it) })
                    }
                    else -> when {
                        context.start > 25000 -> assert(hapids.all { lineAHapids.contains(it) })
                        else -> assert(hapids.any { lineAHapids.contains(it) } && hapids.any { lineBHapids.contains(it) })
                    }
                }
            }
        }
    }

    @Test
    fun testUnorderedHaplotypePair() {

        val test1 = DiploidEmissionProbability.UnorderedHaplotypePair(Pair("aaa", "bbb"))
        val test2 = DiploidEmissionProbability.UnorderedHaplotypePair(Pair("bbb", "aaa"))
        val test3 = DiploidEmissionProbability.UnorderedHaplotypePair(Pair("aaa", "ccc"))

        assert(test1 == test1)
        assert(test1 == test2)
        assert(test1 != test3)

        val mySet = setOf(test1,test1,test2,test2,test3,test3)
        assertEquals(2, mySet.size)

    }

    /**
     * Create a set of read mappings for testing, write them to a file, and return the read counts for LineA and
     * total read counts. The read mappings should be (lineA,lineA) for the first 9 reference ranges of chromosome 1 then
     * (lineA,lineB). Chromosome 2 should be (lineA, lineB) for the first 9 ranges then (lineA,lineA).
     */
    private fun createReadMappings(graph: HaplotypeGraph, mappingFile: String): IntArray {
        val probMapToTarget = 0.99
        val probMapToOther = 0.1
        var readsWithA = 0

        //create a set of read mappings that are all lineA
        //then create set that is 5 ranges A, the rest B for chr 1 and 5 ranges B, the rest A for chr 2
        //then merge the two sets
        val readMap1 = mutableMapOf<List<String>, Int>()
        for (range in graph.ranges()) {
            val hapids = graph.hapIdToSampleGametes(range)
            //generate some mappings
            repeat(3) {
                val hapidList = mutableListOf<String>()
                for ((hapid, samples) in hapids.entries) {
                    val isLineA = samples.any { it.name == "LineA" }

                    if (isLineA && Random.nextDouble() < probMapToTarget) {
                        hapidList.add(hapid)
                        readsWithA += 2  //because the code adds a 2 count for each hapid list
                    }
                    if (!isLineA && Random.nextDouble() < probMapToOther) hapidList.add(hapid)
                }
                if (hapidList.size > 0) {
                    val reads = readMap1[hapidList] ?: 0
                    readMap1[hapidList] = reads + 2
                }
            }
        }

        val readMap2 = mutableMapOf<List<String>, Int>()
        for (range in graph.ranges()) {
            val hapids = graph.hapIdToSampleGametes(range)
            val target = when {
                range.contig == "1" && range.start <= 25000 -> "LineA"
                range.contig == "1" && range.start > 25000 -> "LineB"
                range.contig == "2" && range.start <= 25000 -> "LineB"
                else -> "LineA"
            }

            //generate some mappings
            repeat(3) {
                val hapidList = mutableListOf<String>()
                for ((hapid, samples) in hapids.entries) {
                    val isTarget = samples.any { it.name == target }
                    val isLineA = samples.any { it.name == "LineA" }

                    if (isTarget && Random.nextDouble() < probMapToTarget) {
                        hapidList.add(hapid)
                        if (isLineA) readsWithA += 2
                    }
                    if (!isTarget && Random.nextDouble() < probMapToOther) {
                        hapidList.add(hapid)
                        if (isLineA) readsWithA += 2
                    }
                }
                if (hapidList.size > 0) {
                    val reads = readMap2[hapidList] ?: 0
                    readMap2[hapidList] = reads + 2
                }
            }
        }

        val readMap = mergeReadMappings(listOf(readMap1, readMap2))
        exportReadMapping(mappingFile, readMap, "TestLine", Pair("file1", "file2"))

        return intArrayOf(readsWithA, readMap.values.sum())
    }


}