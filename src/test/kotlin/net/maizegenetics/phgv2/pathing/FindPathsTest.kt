package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.testing.test
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.getBufferedReader
import net.maizegenetics.phgv2.utils.getBufferedWriter
import org.junit.jupiter.api.*
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import java.nio.file.Paths
import kotlin.random.Random
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class FindPathsTest {
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
        var result = FindPaths().test("")
        //expected errors:
        //Error: missing option --path-keyfile
        //Error: missing option --hvcf-dir
        //Error: missing option --reference-genome
        //Error: missing option --output-dir
        var errOut = result.stderr
        assertTrue(errOut.contains("Error: missing option --path-keyfile"))
        assertTrue(errOut.contains("Error: missing option --hvcf-dir"))
        assertTrue(errOut.contains("Error: missing option --reference-genome"))
        assertTrue(errOut.contains("Error: missing option --output-dir"))
        assertTrue(errOut.contains("Error: missing option --path-type"))

        //test wrong value for path-type
        var testArgs = "--path-type triploid --path-keyfile notafile --hvcf-dir notafile --reference-genome notafile --output-dir notafile"
        result = FindPaths().test(testArgs)
        errOut = result.stderr
        assertTrue(errOut.contains("Error: invalid value for --path-type: invalid choice: triploid. (choose from haploid, diploid)"))

        //test for invalid file names
        testArgs = "--path-type haploid --path-keyfile notafile --hvcf-dir notafile --reference-genome notafile --output-dir notafile"
        result = FindPaths().test(testArgs)
        errOut = result.stderr

        //expected errors:
        // Error: invalid value for --path-keyfile: notafile is not a valid file
        //Error: invalid value for --hvcf-dir: notafile is not a valid directory.
        //Error: invalid value for --reference-genome: notafile is not a valid file
        //Error: invalid value for --output-dir: notafile is not a valid directory.
        assertTrue(errOut.contains("Error: invalid value for --path-keyfile: notafile is not a valid file"))
        assertTrue(errOut.contains("Error: invalid value for --hvcf-dir: notafile is not a valid directory."))
        assertTrue(errOut.contains("Error: invalid value for --reference-genome: notafile is not a valid file"))
        assertTrue(errOut.contains("Error: invalid value for --output-dir: notafile is not a valid directory."))

        //test validation of optional parameters
        testArgs = "--path-keyfile ${TestExtension.testKeyFile} --hvcf-dir ${TestExtension.smallSeqInputDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.smallSeqInputDir} " +
                "--path-type haploid --prob-correct 0.4 " +
                "--min-gametes -1 --min-reads -1 --max-reads-per-kb -1 --min-coverage 0.4 --inbreeding-coefficient -1"

        result = FindPaths().test(testArgs)

        //Expected errors:
//        Error: invalid value for --prob-correct: prob-correct must be between 0.5 and 1.0
//        Error: invalid value for --min-gametes: min-gametes must be a positive integer
//        Error: invalid value for --min-reads: min-reads must be a positive integer.
//        Error: invalid value for --inbreeding-coefficient: inbreeding-coefficient must be between 0.0 and 1.0
//        Error: invalid value for --max-reads-per-kb: max-reads-per-kb must be a positive integer.
//        Error: invalid value for --min-coverage: min-coverage must be between 0.5 and 1.0

        errOut = result.stderr
        println(errOut)
        assertTrue(errOut.contains("Error: invalid value for --prob-correct: prob-correct must be between 0.5 and 1.0"))
        assertTrue(errOut.contains("Error: invalid value for --min-gametes: min-gametes must be a positive integer"))
        assertTrue(errOut.contains("Error: invalid value for --min-reads: min-reads must be a positive integer."))
        assertTrue(errOut.contains("Error: invalid value for --max-reads-per-kb: max-reads-per-kb must be a positive integer."))
        assertTrue(errOut.contains("Error: invalid value for --min-coverage: min-coverage must be between 0.5 and 1.0"))
        assertTrue(errOut.contains("Error: invalid value for --inbreeding-coefficient: inbreeding-coefficient must be between 0.0 and 1.0"))

    }

    @Test
    fun testHaploidPathFinding() {
        //plan for this test:
        //build a haplotype groph from Ref, lineA, and lineB
        //for that need an hvcfdir with the files in it
        //make sure the test output directory exists
        val vcfDir = File(TestExtension.testVCFDir)

        //create a read mapping file
        val listOfHvcfFilenames = vcfDir.listFiles().map { it.path }.filter { it.endsWith(".h.vcf") }
        val myGraph = HaplotypeGraph(listOfHvcfFilenames)
        val readMappingFile = TestExtension.testOutputDir + "testReadMapping.txt"
        createHaploidReadMappings(myGraph, readMappingFile)

        //create a keyfile
        val keyFilename = TestExtension.testOutputDir + "keyfileForPathTest.txt"
        getBufferedWriter(keyFilename).use { myWriter ->
            myWriter.write("SampleName\tReadMappingFiles\n")
            myWriter.write("TestLine\t$readMappingFile\n")
        }

        var pathFindingTestArgs = "--path-type haploid --path-keyfile $keyFilename --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir}"

        val pathFindingResult = FindPaths().test(pathFindingTestArgs)
        assertEquals(0, pathFindingResult.statusCode, "pathFinding status code was ${pathFindingResult.statusCode}")


        //are all the haplotypes in the path from LineA?
        var resultHvcfName = "${TestExtension.testOutputDir}TestLine.h.vcf"
        var vcfReader = VCFFileReader(File(resultHvcfName), false)

        val sampleA = SampleGamete("LineA")
        val haplotypesA = myGraph.ranges().map { myGraph.sampleToHapId(it, sampleA) }
        for (record in vcfReader) {
            val testHaplotype = record.genotypes["TestLine"].alleles[0].displayString
                .substringAfter("<").substringBefore(">")
            //this next line deals with null haplotypes, which in haplotypes A is null but is "" from VCF reader
            if (testHaplotype.isNotBlank()) assert(haplotypesA.contains(testHaplotype)) {"$testHaplotype not a LineA haplotype"}
        }
        vcfReader.close()

        //test minGametes, minReads, maxReadsPerKb
        //set minReads = 10. No file expected.
        //change the sample name
        getBufferedWriter(keyFilename).use { myWriter ->
            myWriter.write("SampleName\tReadMappingFiles\n")
            myWriter.write("TestLineMinReads\t$readMappingFile\n")
        }

        pathFindingTestArgs = "--path-type haploid --path-keyfile $keyFilename --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir} " +
                "--min-reads 10"
        val minReadTestResult = FindPaths().test(pathFindingTestArgs)
        assertEquals(0, minReadTestResult.statusCode)
        assert(!File("${TestExtension.testOutputDir}TestLineMinReads.h.vcf").exists())

        //set minReads = 1, minGametes = 4. This should generate an error message.
        //change the sample name
        getBufferedWriter(keyFilename).use { myWriter ->
            myWriter.write("SampleName\tReadMappingFiles\n")
            myWriter.write("TestLineMinGametes\t$readMappingFile\n")
        }
        pathFindingTestArgs = "--path-type haploid --path-keyfile $keyFilename --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir} " +
                "--min-reads 1 --min-gametes 4"
        val minGametesTestResult = FindPaths().test(pathFindingTestArgs)
        assertEquals(0, minGametesTestResult.statusCode)
        assert(!File("${TestExtension.testOutputDir}TestLineMinGametes.h.vcf").exists())

        //set minReads = 1. maxReadsPerKb = 4. This should filter out about half the reference ranges.
        //change the sample name
        getBufferedWriter(keyFilename).use { myWriter ->
            myWriter.write("SampleName\tReadMappingFiles\n")
            myWriter.write("TestLineMaxReads\t$readMappingFile\n")
        }
        pathFindingTestArgs = "--path-type haploid --path-keyfile $keyFilename --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir} " +
                "--min-reads 1 --max-reads-per-kb 4"
        val maxReadTestResults = FindPaths().test(pathFindingTestArgs)

        assertEquals(0, maxReadTestResults.statusCode)
        resultHvcfName = "${TestExtension.testOutputDir}TestLineMaxReads.h.vcf"
        vcfReader = VCFFileReader(File(resultHvcfName), false)
        val numberOfRecords = vcfReader.count()
        assertTrue(numberOfRecords in 18..20)
        vcfReader.close()


        //run a test with recombination (path switching)
        //create a read mapping file and a new keyfile
        //create a keyfile
        val switchKeyFilename = TestExtension.testOutputDir + "keyfileForPathTest.txt"
        getBufferedWriter(switchKeyFilename).use { myWriter ->
            myWriter.write("SampleName\tReadMappingFiles\n")
            myWriter.write("TestLine2\t$readMappingFile\n")
        }

        createReadMappingsWithPathSwitches(myGraph, readMappingFile)
        val switchTestArgs = "--path-type haploid --path-keyfile $switchKeyFilename --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir} "

        val switchResult = FindPaths().test(switchTestArgs)
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
            if (testHaplotype.isNotBlank()) {
                if (record.start < 25000) assert(haplotypesA.contains(testHaplotype)) { "$testHaplotype for ${record.contig}:${record.start} not a LineA haplotype" }
                else assert(haplotypesB.contains(testHaplotype)) { "$testHaplotype for ${record.contig}:${record.start} not a LineB haplotype" }
            }
        }
        vcfReader.close()



    }

    @Test
    fun testDiploidPathFinding() {
        //use a haplotype groph built from Ref, lineA, and lineB
        val vcfDir = File(TestExtension.testVCFDir)
        val listOfHvcfFilenames = vcfDir.listFiles().map { it.path }.filter { it.endsWith(".h.vcf") }
        val myGraph = HaplotypeGraph(listOfHvcfFilenames)

        //create a read mapping file
        val readMappingFile = TestExtension.testOutputDir + "testReadMapping_diploid.txt"
        createDiploidReadMappings(myGraph, readMappingFile)

        //create a keyfile
        val lineName = "TestLineD"
        val keyFilename = TestExtension.testOutputDir + "keyfileForDiploidPathTest.txt"
        getBufferedWriter(keyFilename).use { myWriter ->
            myWriter.write("SampleName\tReadMappingFiles\n")
            myWriter.write("$lineName\t$readMappingFile\n")
        }

        val pathFindingTestArgs = "--path-keyfile $keyFilename --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir} " +
                "--path-type diploid --prob-same-gamete 0.8"

        val pathFindingResult = FindPaths().test(pathFindingTestArgs)
        assert(pathFindingResult.statusCode == 0)

        //examine the resulting hvcf file
        val testVcf = "${TestExtension.testOutputDir}${lineName}.h.vcf"
        val lineAgamete = SampleGamete("LineA")
        val lineBgamete = SampleGamete("LineB")
        val refGamete = SampleGamete("Ref")
        val lineAHapids = myGraph.ranges().mapNotNull { range -> myGraph.sampleToHapId(range, lineAgamete) }.toSet()
        val lineBHapids = myGraph.ranges().mapNotNull { range -> myGraph.sampleToHapId(range, lineBgamete) }.toSet()
        val refHapids = myGraph.ranges().mapNotNull { range -> myGraph.sampleToHapId(range, refGamete) }.toSet()

        VCFFileReader(Paths.get(testVcf), false).use { vcf ->
            //Chr 1: one haplotype should be all A for chr1, the other A before start=25000 and ref after
            //Chr 2: one haplotype should be all A for chr1, the other B before start=25000 and A after
            for (context in vcf) {
                //filter out the ref allele
                val hapids = context.alternateAlleles.map { it.displayString.substringBefore(">").substringAfter("<") }.toList()
                when (context.contig) {
                    "1" -> when {
                        context.start < 25000 -> assertTrue(
                            hapids.all { lineAHapids.contains(it) },
                            "${context.contig}:${context.start}-${context.end}, hapids = $hapids")
                        context.start in (25000..50000) -> assertTrue(
                            hapids.any { lineAHapids.contains(it) } && hapids.any { refHapids.contains(it) },
                            "${context.contig}:${context.start}-${context.end}, hapids = $hapids")
                        else -> assertTrue(
                            hapids.all { refHapids.contains(it) },
                            "${context.contig}:${context.start}-${context.end}, hapids = $hapids")
                    }
                    else -> when {
                        context.start > 25000 -> assertTrue(
                            hapids.all { lineAHapids.contains(it) },
                            "${context.contig}:${context.start}-${context.end}, hapids = $hapids")
                        else -> assertTrue(
                            hapids.any { lineAHapids.contains(it) } && hapids.any { lineBHapids.contains(it) },
                            "${context.contig}:${context.start}-${context.end}, hapids = $hapids")
                    }
                }
            }
        }

        //test for exception when --use-likely-ancestors true
        val pathFindingTestArgsLikelyAncestor = "--path-type diploid --path-keyfile $keyFilename --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir} " +
                "--prob-same-gamete 0.8 --use-likely-ancestors true"

        val exception = assertThrows<IllegalArgumentException> { FindPaths().test(pathFindingTestArgsLikelyAncestor) }
        assert(exception.message!!.startsWith("UseLikelyAncestors is true but likelyAncestors will not be checked"))
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
        val readCounts = createHaploidReadMappings(myGraph, readMappingFile)
        val expectedCoverage = readCounts[0].toDouble() / readCounts[1].toDouble()

        //create a keyfile
        val keyFilename = TestExtension.testOutputDir + "keyfileForPathTest.txt"
        getBufferedWriter(keyFilename).use { myWriter ->
            myWriter.write("SampleName\tReadMappingFiles\n")
            myWriter.write("TestLine\t$readMappingFile\n")
        }

        val likelyAncestorFile = TestExtension.testOutputDir + "testAncestors.txt"
        File(likelyAncestorFile).delete()
        val pathFindingTestArgs = "--path-type haploid --path-keyfile $keyFilename --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir} " +
                "--use-likely-ancestors true --max-ancestors 2 --min-coverage 0.9 --likely-ancestor-file $likelyAncestorFile"

        val pathFindingResult = FindPaths().test(pathFindingTestArgs)
        assertEquals(0, pathFindingResult.statusCode, "pathFinding status code was ${pathFindingResult.statusCode}")

        //check contents of likelyAncestorFile
        getBufferedReader(likelyAncestorFile).use {
            val header = it.readLine()
            assertEquals("sample\tancestor\tgameteId\treads\tcoverage", header)
            val parentInfo = it.readLine().split("\t")
            assertEquals("TestLine", parentInfo[0])
            assertEquals("LineA", parentInfo[1])
            assertEquals("0", parentInfo[2])
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
            if (testHaplotype.isNotBlank()) assert(haplotypesA.contains(testHaplotype)) {"$testHaplotype not a LineA haplotype"}
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
        Assertions.assertEquals(2, mySet.size)
    }

    /**
     * Create a set of read mappings for testing, write them to a file, and return the read counts for LineA and
     * total read counts. The read mappings should be (lineA,lineA) for the first 9 reference ranges of chromosome 1 then
     * (lineA,lineB). Chromosome 2 should be (lineA, lineB) for the first 9 ranges then (lineA,lineA).
     */
    private fun createDiploidReadMappings(graph: HaplotypeGraph, mappingFile: String): IntArray {
        val probMapToTarget = 0.99
        val probMapToOther = 0.1
        var readsWithA = 0
        val repetitions = 5

        //create a set of read mappings that are all lineA
        //then create set that is 5 ranges A, the rest B for chr 1 and 5 ranges B, the rest A for chr 2
        //then merge the two sets
        //do not generate any reads for ranges that have no A or B haplotype
        val readMap1 = mutableMapOf<List<String>, Int>()

        for (range in graph.ranges()) {
            val hapids = graph.hapIdToSampleGametes(range)
            //generate some mappings if the range has a LineA haplotype
            val hasLineA = hapids.values.any{it.contains(SampleGamete("LineA"))}

            if (hasLineA) repeat(repetitions) {
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
                range.contig == "1" && range.start > 25000 -> "Ref"
                range.contig == "2" && range.start <= 25000 -> "LineB"
                else -> "LineA"
            }
            val hasTarget = hapids.values.any{it.contains(SampleGamete(target))}

            //generate some mappings
            if (hasTarget) repeat(repetitions) {
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

    /**
     * Create a set of read mappings for testing, write them to a file, and return the read counts for LineA and
     * total read counts
     */
    private fun createHaploidReadMappings(graph: HaplotypeGraph, mappingFile: String): IntArray {
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