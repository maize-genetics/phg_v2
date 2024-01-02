package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.random.Random


@ExtendWith(TestExtension::class)
class HaploidPathFindingTest {

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
        val vcfDir = File(TestExtension.testVCFDir)
        vcfDir.mkdirs()

        //erase any files in testVCFDir
        val vcfFiles = vcfDir.listFiles(){name -> name.endsWith("vcf")}
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

    fun getReadMappings(hvcfDir: File): Map<List<String>, Int> {
        val probMapToA = 0.95
        val probMapToOther = 0.5
        val listOfHvcfFilenames = hvcfDir.listFiles().map { it.path }
        val graph = HaplotypeGraph(listOfHvcfFilenames)
        val readMap = mutableMapOf<List<String>, Int>()
        for (range in graph.ranges()) {
            val hapids = graph.hapIdToSampleGametes(range)
            //generate some mappings
            repeat(4) {
                val hapidList = mutableListOf<String>()
                for ((hapid, samples) in hapids.entries) {
                    val isLineA = samples.any { it.name == "LineA" }

                    if (isLineA && Random.nextDouble() < probMapToA) hapidList.add("LineA")
                    if (!isLineA && Random.nextDouble() < probMapToOther) hapidList.add(samples[0].name)
                }
                if (hapidList.size > 0) readMap.put(hapidList, 2)
            }
        }

        return readMap
    }
}