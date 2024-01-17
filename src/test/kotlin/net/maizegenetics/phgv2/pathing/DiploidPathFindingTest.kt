package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
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
            val vcfFiles = vcfDir.listFiles { file -> file.name.endsWith(".h.vcf")}
            for (vcfFile in vcfFiles) {
                vcfFile.delete()
            }

            //copy hvcf files to vcfDir
            listOf(TestExtension.smallseqLineAHvcfFile, TestExtension.smallseqLineBHvcfFile, TestExtension.smallseqRefHvcfFile)
                .forEach { filepath ->
                    val origFile = File(filepath)
                    origFile.copyTo(vcfDir.resolve(origFile.name))
                }

        }

        @JvmStatic
        @AfterAll
        fun tearDown() {
            //comment out the following line to inspect the test results after the tests have been run
            File(TestExtension.testOutputDir).deleteRecursively()
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
        result = HaploidPathFinding().test(testArgs)
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

        result = HaploidPathFinding().test(testArgs)
        //Expected errors:
        //Error: invalid value for --prob-correct: prob-correct must be between 0.5 and 1.0
        //Error: invalid value for --min-gametes: min-gametes must be a positive integer
        //Error: invalid value for --min-reads: min-reads must be a positive integer.
        //Error: invalid value for --max-reads-per-kb: max-reads-per-kb must be a positive integer.
        //Error: invalid value for --min-coverage: min-coverage must be between 0.5 and 1.0
        errOut = result.stderr
        assert(errOut.contains("Error: invalid value for --prob-correct: prob-correct must be between 0.5 and 1.0"))
        assert(errOut.contains("Error: invalid value for --min-gametes: min-gametes must be a positive integer"))
        assert(errOut.contains("Error: invalid value for --min-reads: min-reads must be a positive integer."))
        assert(errOut.contains("Error: invalid value for --max-reads-per-kb: max-reads-per-kb must be a positive integer."))
        assert(errOut.contains("Error: invalid value for --min-coverage: min-coverage must be between 0.5 and 1.0"))
        assert(errOut.contains("Error: invalid value for --inbreeding-coefficient: inbreeding-coefficient must be a positive double"))

    }

    @Test
    fun testDiploidPathFinding() {

    }

    /**
     * Create a set of read mappings for testing, write them to a file, and return the read counts for LineA and
     * total read counts. The read mappings should be (lineA,lineA) for the first 9 reference ranges of chromosome 1 then
     * (lineA,lineB). Chromosome 2 should be (lineA, lineB) for the first 9 ranges then (lineA,lineA).
     */
    private fun createReadMappings(graph: HaplotypeGraph, mappingFile: String): IntArray {
        val probMapToTarget = 0.99
        val probMapToOther = 0.2
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
                    val isLineA = samples.any {it.name == "LineA"}

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