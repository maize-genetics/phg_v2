package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.testing.test
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.getBufferedReader
import net.maizegenetics.phgv2.utils.getBufferedWriter
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.random.Random
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class HaploidPathFindingTest {

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
                "--min-gametes -1 --min-reads -1 --max-reads-per-kb -1 --min-coverage 0.4"

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

    }

    @Test
    fun testPathFinding() {
        //plan for this test:
        //build a haplotype groph from Ref, lineA, and lineB
        //for that need an hvcfdir with the files in it
        //make sure the test output directory exists
        val vcfDir = File(TestExtension.testVCFDir)

        //create a read mapping file
        val listOfHvcfFilenames = vcfDir.listFiles().map { it.path }.filter { it.endsWith(".h.vcf") }
        val myGraph = HaplotypeGraph(listOfHvcfFilenames)
        val readMappingFile = TestExtension.testOutputDir + "testReadMapping.txt"
        createReadMappings(myGraph, readMappingFile)

        //create a keyfile
        val keyFilename = TestExtension.testOutputDir + "keyfileForPathTest.txt"
        getBufferedWriter(keyFilename).use { myWriter ->
            myWriter.write("SampleName\tReadMappingFiles\n")
            myWriter.write("TestLine\t$readMappingFile\n")
        }

        var pathFindingTestArgs = "--path-keyfile $keyFilename --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir}"

        val pathFindingResult = HaploidPathFinding().test(pathFindingTestArgs)
        assertEquals(0, pathFindingResult.statusCode, "pathFinding status code was ${pathFindingResult.statusCode}")


        //are all the haplotypes in the path from LineA?
        var resultHvcfName = "${TestExtension.testOutputDir}TestLine.h.vcf"
        var vcfReader = VCFFileReader(File(resultHvcfName), false)

        val sampleA = SampleGamete("LineA")
        val haplotypesA = myGraph.ranges().map { myGraph.sampleToHapId(it, sampleA) }
        for (record in vcfReader) {
            val testHaplotype = record.genotypes["TestLine"].alleles[0].displayString
                .substringAfter("<").substringBefore(">")
            assert(haplotypesA.contains(testHaplotype)) {"$testHaplotype not a LineA haplotype"}
        }
        vcfReader.close()

        //test minGametes, minReads, maxReadsPerKb
        //set minReads = 10. This should generate an error message.
        //change the sample name
        getBufferedWriter(keyFilename).use { myWriter ->
            myWriter.write("SampleName\tReadMappingFiles\n")
            myWriter.write("TestLineMinReads\t$readMappingFile\n")
        }

        pathFindingTestArgs = "--path-keyfile $keyFilename --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir} " +
                "--min-reads 10"
        val minReadTestResult = HaploidPathFinding().test(pathFindingTestArgs)
        assertEquals(0, minReadTestResult.statusCode)
        assert(!File("${TestExtension.testOutputDir}TestLineMinReads.h.vcf").exists())

        //set minReads = 1, minGametes = 4. This should generate an error message.
        //change the sample name
        getBufferedWriter(keyFilename).use { myWriter ->
            myWriter.write("SampleName\tReadMappingFiles\n")
            myWriter.write("TestLineMinGametes\t$readMappingFile\n")
        }
        pathFindingTestArgs = "--path-keyfile $keyFilename --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir} " +
                "--min-reads 1 --min-gametes 4"
        val minGametesTestResult = HaploidPathFinding().test(pathFindingTestArgs)
        assertEquals(0, minGametesTestResult.statusCode)
        assert(!File("${TestExtension.testOutputDir}TestLineMinGametes.h.vcf").exists())

        //set minReads = 1. maxReadsPerKb = 4. This should filter out about half the reference ranges.
        //change the sample name
        getBufferedWriter(keyFilename).use { myWriter ->
            myWriter.write("SampleName\tReadMappingFiles\n")
            myWriter.write("TestLineMaxReads\t$readMappingFile\n")
        }
        pathFindingTestArgs = "--path-keyfile $keyFilename --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir} " +
                "--min-reads 1 --max-reads-per-kb 4"
        val maxReadsTestResult = HaploidPathFinding().test(pathFindingTestArgs)
        assertEquals(0, minReadTestResult.statusCode)
        resultHvcfName = "${TestExtension.testOutputDir}TestLineMaxReads.h.vcf"
        vcfReader = VCFFileReader(File(resultHvcfName), false)
        val numberOfRecords = vcfReader.count()
        assertEquals(16, numberOfRecords)
        vcfReader.close()


        //run same test with recombination (path switching)
        //create a read mapping file and a new keyfile
        //create a keyfile
        val switchKeyFilename = TestExtension.testOutputDir + "keyfileForPathTest.txt"
        getBufferedWriter(keyFilename).use { myWriter ->
            myWriter.write("SampleName\tReadMappingFiles\n")
            myWriter.write("TestLine2\t$readMappingFile\n")
        }

        createReadMappingsWithPathSwitches(myGraph, readMappingFile)
        val switchResult = HaploidPathFinding().test(pathFindingTestArgs)
        assertEquals(0, switchResult.statusCode, "pathFinding status code was ${pathFindingResult.statusCode}")

        //are the haplotypes from the expected target line
        //if range.start < 25000, LineA else LineB
        resultHvcfName = "${TestExtension.testOutputDir}TestLine2.h.vcf"
        vcfReader = VCFFileReader(File(resultHvcfName), false)

        val sampleB = SampleGamete("LineB")
        val haplotypesB = myGraph.ranges().map { myGraph.sampleToHapId(it, sampleB) }
        for (record in vcfReader) {
            val testHaplotype = record.genotypes["TestLine2"].alleles[0].displayString
                .substringAfter("<").substringBefore(">")
            if (record.start < 25000) assert(haplotypesA.contains(testHaplotype)) {"$testHaplotype for ${record.contig}:${record.start} not a LineA haplotype"}
            else assert(haplotypesB.contains(testHaplotype)) {"$testHaplotype for ${record.contig}:${record.start} not a LineB haplotype"}
        }
        vcfReader.close()


    }

    @Test
    fun testPathFindingWithLikelyParents() {
        //plan for this test:
        //build a haplotype groph from Ref, lineA, and lineB
        //for that need an hvcfdir with the files in it
        //make sure the test output directory exists
        val vcfDir = File(TestExtension.testVCFDir)

        //create a read mapping file
        val listOfHvcfFilenames = vcfDir.listFiles().map { it.path }.filter { it.endsWith(".h.vcf") }
        val myGraph = HaplotypeGraph(listOfHvcfFilenames)
        val readMappingFile = TestExtension.testOutputDir + "testReadMapping.txt"
        val readCounts = createReadMappings(myGraph, readMappingFile)
        val countOfA = readCounts[0]
        val expectedCoverage = readCounts[0].toDouble() / readCounts[1].toDouble()

        //create a keyfile
        val keyFilename = TestExtension.testOutputDir + "keyfileForPathTest.txt"
        getBufferedWriter(keyFilename).use { myWriter ->
            myWriter.write("SampleName\tReadMappingFiles\n")
            myWriter.write("TestLine\t$readMappingFile\n")
        }

        val likelyAncestorFile = TestExtension.testOutputDir + "testAncestors.txt"
        File(likelyAncestorFile).delete()
        val pathFindingTestArgs = "--path-keyfile $keyFilename --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir} " +
                "--use-likely-ancestors true --max-ancestors 2 --min-coverage 0.9 --likely-ancestor-file $likelyAncestorFile"

        val pathFindingResult = HaploidPathFinding().test(pathFindingTestArgs)
        assertEquals(0, pathFindingResult.statusCode, "pathFinding status code was ${pathFindingResult.statusCode}")

        //check contents of likelyAncestorFile
        getBufferedReader(likelyAncestorFile).use {
            val header = it.readLine()
            assertEquals("sample\tancestor\tgameteId\treads\tcoverage", header)
            val parentInfo = it.readLine().split("\t")
            assertEquals("TestLine", parentInfo[0])
            assertEquals("LineA", parentInfo[1])
            assertEquals("0", parentInfo[2])
            assertEquals(countOfA, parentInfo[3].toInt())
            assertEquals(expectedCoverage, parentInfo[4].toDouble())
        }



        //are all the haplotypes in the path from LineA?
        val resultHvcfName = "${TestExtension.testOutputDir}TestLine.h.vcf"
        val vcfReader = VCFFileReader(File(resultHvcfName), false)

        val sampleA = SampleGamete("LineA")
        val haplotypesA = myGraph.ranges().map { myGraph.sampleToHapId(it, sampleA) }
        for (record in vcfReader) {
            val testHaplotype = record.genotypes["TestLine"].alleles[0].displayString
                .substringAfter("<").substringBefore(">")
            assert(haplotypesA.contains(testHaplotype)) {"$testHaplotype not a LineA haplotype"}
        }

    }


    /**
     * Create a set of read mappings for testing, write them to a file, and return the read counts for LineA and
     * total read counts
     */
    private fun createReadMappings(graph: HaplotypeGraph, mappingFile: String): IntArray {
        val probMapToA = 0.99
        val probMapToOther = 0.2
        val readMap = mutableMapOf<List<String>, Int>()
        var readsWithA = 0
        for (range in graph.ranges()) {
            val hapids = graph.hapIdToSampleGametes(range)
            //generate some mappings
            repeat(4) {
                val hapidList = mutableListOf<String>()
                for ((hapid, samples) in hapids.entries) {
                    val isLineA = samples.any { it.name == "LineA" }

                    if (isLineA && Random.nextDouble() < probMapToA) {
                        hapidList.add(hapid)
                        readsWithA += 2  //because the code adds a 2 count for each hapid list
                    }
                    if (!isLineA && Random.nextDouble() < probMapToOther) hapidList.add(hapid)
                }
                if (hapidList.size > 0) {
                    val reads = readMap[hapidList] ?: 0
                    readMap[hapidList] = reads + 2
                }
            }
        }

        exportReadMapping(mappingFile, readMap, "TestLine", Pair("file1", "file2"))

        val totalCount = readMap.entries.sumOf { (_, count) -> count }
        return intArrayOf(readsWithA, totalCount)
    }

    private fun createReadMappingsWithPathSwitches(graph: HaplotypeGraph, mappingFile: String) {
        //additional conditions to test: minReads, maxReads


        val probMapToTarget = 0.99
        val probMapToOther = 0.2
        val readMap = mutableMapOf<List<String>, Int>()
        var readsWithA = 0
        for (range in graph.ranges()) {
            val hapids = graph.hapIdToSampleGametes(range)
            val target = if (range.start > 25000) "LineB" else "LineA"

            //generate some mappings
            repeat(4) {
                val hapidList = mutableListOf<String>()
                for ((hapid, samples) in hapids.entries) {
                    val isTarget = samples.any { it.name == target }

                    if (isTarget && Random.nextDouble() < probMapToTarget) {
                        hapidList.add(hapid)
                    }
                    if (!isTarget && Random.nextDouble() < probMapToOther) hapidList.add(hapid)
                }
                if (hapidList.size > 0) {
                    val reads = readMap[hapidList] ?: 0
                    readMap[hapidList] = reads + 2
                }
            }
        }

        exportReadMapping(mappingFile, readMap, "TestLine", Pair("file1", "file2"))

    }
}