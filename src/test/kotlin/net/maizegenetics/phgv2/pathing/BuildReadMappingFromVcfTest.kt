package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.testing.test
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.Position
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertTrue

/**
 * Test fixtures under data/test/vcfReadMapping/:
 *
 * index.txt covers:
 *   1:10  A -> [h1]        T -> [h1,h2]
 *   1:20  C -> [h3]                          (no entry for "G" - tests an unmatched allele)
 *   1:30  A -> [h6]        G -> [h4,h5]
 *   1:99  A -> [h7]                          (never touched by the impute VCF)
 *
 * imputeVcf.vcf has samples SampleX, SampleY, SampleZ, SampleEmpty:
 *   1:10 REF=A ALT=T  SampleX 0/1:12,8  SampleY 1/1 (no AD)  SampleZ 0/0 (no AD)  SampleEmpty ./.
 *   1:20 REF=C ALT=G  SampleX 0/0       SampleY 0/0          SampleZ 0/1 (no AD) SampleEmpty ./.
 *   1:30 REF=A ALT=G  SampleX 0/1:9 (AD too short for the ALT index)
 *                     SampleY 1/1:6,0 (ALT AD is 0)          SampleZ 0/0 (no AD) SampleEmpty ./.
 *   1:40 REF=A ALT=T  (not in the index at all) - all samples called, all skipped
 *
 * Hand-derived expected results (see comments on each test) drive both the default-mode and
 * --use-allele-depth-mode assertions below.
 */
@ExtendWith(TestExtension::class)
class BuildReadMappingFromVcfTest {

    companion object {
        const val indexFile = "data/test/vcfReadMapping/index.txt"
        const val imputeVcf = "data/test/vcfReadMapping/imputeVcf.vcf"

        const val step1HvcfDir = "data/test/vcfHaplotypeIndex/hvcfs"
        const val step1PanelVcf = "data/test/vcfHaplotypeIndex/testPanel.vcf"

        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.testOutputDir).deleteRecursively()
            File(TestExtension.testOutputDir).mkdirs()
        }

        @JvmStatic
        @AfterAll
        fun tearDown() {
            //comment out the following line to inspect the test results after the tests have been run
            File(TestExtension.testOutputDir).deleteRecursively()
        }

        fun imputeRecords(): Map<Int, VariantContext> {
            return VCFFileReader(File(imputeVcf), false).use { reader -> reader.associateBy { it.start } }
        }

        fun index(): Map<Position, VcfHaplotypeIndexPosition> = VcfHaplotypeIndexUtils.readIndex(indexFile)

        fun readMappingAsMap(file: String): Map<List<String>, Int> = AlignmentUtils.importReadMapping(file)
    }

    //region CLI validation

    @Test
    fun testCliktParams() {
        val cmd = BuildReadMappingFromVcf()

        val noIndexFile = cmd.test("--impute-vcf $imputeVcf --output-dir ${TestExtension.testOutputDir}")
        assertEquals(1, noIndexFile.statusCode)
        assertEquals(
            "Usage: build-read-mapping-from-vcf [<options>]\n\nError: missing option --index-file\n",
            noIndexFile.stderr
        )

        val noImputeVcf = cmd.test("--index-file $indexFile --output-dir ${TestExtension.testOutputDir}")
        assertEquals(1, noImputeVcf.statusCode)
        assertEquals(
            "Usage: build-read-mapping-from-vcf [<options>]\n\nError: missing option --impute-vcf\n",
            noImputeVcf.stderr
        )

        val noOutputDir = cmd.test("--index-file $indexFile --impute-vcf $imputeVcf")
        assertEquals(1, noOutputDir.statusCode)
        assertEquals(
            "Usage: build-read-mapping-from-vcf [<options>]\n\nError: missing option --output-dir\n",
            noOutputDir.stderr
        )
    }

    @Test
    fun testCliktParamsInvalidPaths() {
        val cmd = BuildReadMappingFromVcf()

        val badIndexFile = cmd.test("--index-file notAFile.txt --impute-vcf $imputeVcf --output-dir ${TestExtension.testOutputDir}")
        assertEquals(1, badIndexFile.statusCode)
        assertTrue(badIndexFile.stderr.contains("is not a valid file"))

        val badImputeVcf = cmd.test("--index-file $indexFile --impute-vcf notAFile.vcf --output-dir ${TestExtension.testOutputDir}")
        assertEquals(1, badImputeVcf.statusCode)
        assertTrue(badImputeVcf.stderr.contains("is not a valid file"))
    }

    //endregion

    //region alleleWeight

    @Test
    fun testAlleleWeightFlatModeAlwaysOne() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(10)
        val genotypeWithAD = record.getGenotype("SampleX")
        val genotypeWithoutAD = record.getGenotype("SampleY")
        val stats = ReadMappingStats()

        assertEquals(1, cmd.alleleWeight(genotypeWithAD, record.getAlternateAllele(0), record, useAlleleDepth = false, stats = stats))
        assertEquals(1, cmd.alleleWeight(genotypeWithoutAD, record.getAlternateAllele(0), record, useAlleleDepth = false, stats = stats))
        assertEquals(0, stats.adFallbackCount)
        assertEquals(0, stats.zeroDepthAlleleCalls)
    }

    @Test
    fun testAlleleWeightUsesADWhenPresentAndInRange() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(10)
        val genotype = record.getGenotype("SampleX") //AD=12,8
        val stats = ReadMappingStats()

        val refWeight = cmd.alleleWeight(genotype, record.reference, record, useAlleleDepth = true, stats = stats)
        val altWeight = cmd.alleleWeight(genotype, record.getAlternateAllele(0), record, useAlleleDepth = true, stats = stats)

        assertEquals(12, refWeight)
        assertEquals(8, altWeight)
        assertEquals(0, stats.adFallbackCount)
    }

    @Test
    fun testAlleleWeightFallsBackWhenADMissing() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(10)
        val genotype = record.getGenotype("SampleY") //no AD
        val stats = ReadMappingStats()

        val weight = cmd.alleleWeight(genotype, record.getAlternateAllele(0), record, useAlleleDepth = true, stats = stats)

        assertEquals(1, weight)
        assertEquals(1, stats.adFallbackCount)
        assertEquals(0, stats.zeroDepthAlleleCalls)
    }

    @Test
    fun testAlleleWeightFallsBackWhenAlleleIndexOutOfBounds() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(30)
        val genotype = record.getGenotype("SampleX") //AD=9, length 1 - only covers the REF allele
        val stats = ReadMappingStats()

        val refWeight = cmd.alleleWeight(genotype, record.reference, record, useAlleleDepth = true, stats = stats)
        assertEquals(9, refWeight)
        assertEquals(0, stats.adFallbackCount)

        val altWeight = cmd.alleleWeight(genotype, record.getAlternateAllele(0), record, useAlleleDepth = true, stats = stats)
        assertEquals(1, altWeight)
        assertEquals(1, stats.adFallbackCount)
    }

    @Test
    fun testAlleleWeightZeroDepthReturnsZero() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(30)
        val genotype = record.getGenotype("SampleY") //AD=6,0 - ALT depth is 0
        val stats = ReadMappingStats()

        val weight = cmd.alleleWeight(genotype, record.getAlternateAllele(0), record, useAlleleDepth = true, stats = stats)

        assertEquals(0, weight)
        assertEquals(1, stats.zeroDepthAlleleCalls)
        assertEquals(0, stats.adFallbackCount)
    }

    //endregion

    //region processRecord

    @Test
    fun testProcessRecordPositionNotInIndex() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(40)
        val stats = ReadMappingStats()

        val deltas = cmd.processRecord(record, index(), useAlleleDepth = false, stats = stats)

        assertTrue(deltas.isEmpty())
        assertEquals(1, stats.positionsNotInIndex)
        assertEquals(0, stats.positionsMatched)
    }

    @Test
    fun testProcessRecordHetCallProducesTwoIndependentIncrements() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(10)
        val stats = ReadMappingStats()

        val deltas = cmd.processRecord(record, index(), useAlleleDepth = false, stats = stats)

        //SampleX is 0/1 (A,T): allele A -> [h1] +1, allele T -> [h1,h2] +1 - two independent
        //buckets in the SAME sample's delta map, proving the confirmed diploid-counting decision.
        val sampleXDelta = deltas.getValue("SampleX")
        assertEquals(mapOf(listOf("h1") to 1, listOf("h1", "h2") to 1), sampleXDelta)
    }

    @Test
    fun testProcessRecordUnmatchedAlleleIsSkipped() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(20)
        val stats = ReadMappingStats()

        val deltas = cmd.processRecord(record, index(), useAlleleDepth = false, stats = stats)

        //SampleZ is 0/1 (C,G): C resolves to [h3], G has no entry in the index at 1:20.
        assertEquals(mapOf(listOf("h3") to 1), deltas.getValue("SampleZ"))
        assertEquals(1, stats.allelesNotInIndex)
    }

    @Test
    fun testProcessRecordNoCallSkipped() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(10)
        val stats = ReadMappingStats()

        val deltas = cmd.processRecord(record, index(), useAlleleDepth = false, stats = stats)

        assertTrue(!deltas.containsKey("SampleEmpty"))
        assertEquals(2, stats.noCallAlleles)
    }

    @Test
    fun testProcessRecordSharedAmbiguousHapIdSet() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(30)
        val stats = ReadMappingStats()

        val deltas = cmd.processRecord(record, index(), useAlleleDepth = false, stats = stats)

        //SampleZ is 0/0 (A,A) at 1:30 -> allele A resolves to the single hapId [h6].
        assertEquals(mapOf(listOf("h6") to 2), deltas.getValue("SampleZ"))
    }

    @Test
    fun testProcessRecordUseAlleleDepthWeighting() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(10)
        val stats = ReadMappingStats()

        val deltas = cmd.processRecord(record, index(), useAlleleDepth = true, stats = stats)

        //SampleX AD=12,8: allele A -> [h1] += 12, allele T -> [h1,h2] += 8 (not the flat +1/+1).
        assertEquals(mapOf(listOf("h1") to 12, listOf("h1", "h2") to 8), deltas.getValue("SampleX"))
    }

    //endregion

    //region buildReadMappings - default mode (hand-derived expected maps; see class KDoc)

    @Test
    fun testBuildReadMappingsDefaultMode() {
        val cmd = BuildReadMappingFromVcf()
        val (perSample, stats) = cmd.buildReadMappings(indexFile, imputeVcf, useAlleleDepth = false)

        assertEquals(mapOf(listOf("h1") to 1, listOf("h1", "h2") to 1, listOf("h3") to 2, listOf("h6") to 1, listOf("h4", "h5") to 1), perSample.getValue("SampleX"))
        assertEquals(mapOf(listOf("h1", "h2") to 2, listOf("h3") to 2, listOf("h4", "h5") to 2), perSample.getValue("SampleY"))
        assertEquals(mapOf(listOf("h1") to 2, listOf("h3") to 1, listOf("h6") to 2), perSample.getValue("SampleZ"))
        //SampleEmpty has zero resolvable evidence but still gets an entry.
        assertEquals(emptyMap<List<String>, Int>(), perSample.getValue("SampleEmpty"))

        assertEquals(4, stats.totalImputeRecords)
        assertEquals(3, stats.positionsMatched)
        assertEquals(1, stats.positionsNotInIndex)
        assertEquals(6, stats.noCallAlleles)
        assertEquals(1, stats.allelesNotInIndex)
        assertEquals(17, stats.gameteHits)
        assertEquals(17L, stats.totalEvidenceUnits)
        assertEquals(0, stats.adFallbackCount)
        assertEquals(0, stats.zeroDepthAlleleCalls)
    }

    //endregion

    //region buildReadMappings - --use-allele-depth mode

    @Test
    fun testBuildReadMappingsUseAlleleDepthMode() {
        val cmd = BuildReadMappingFromVcf()
        val (perSample, stats) = cmd.buildReadMappings(indexFile, imputeVcf, useAlleleDepth = true)

        assertEquals(
            mapOf(listOf("h1") to 12, listOf("h1", "h2") to 8, listOf("h3") to 2, listOf("h6") to 9, listOf("h4", "h5") to 1),
            perSample.getValue("SampleX")
        )
        //SampleY's [h4,h5] contribution is gone entirely - both G calls at 1:30 had AD=0.
        assertEquals(mapOf(listOf("h1", "h2") to 2, listOf("h3") to 2), perSample.getValue("SampleY"))
        //SampleZ never has AD anywhere in this fixture, so AD mode falls back to flat counts
        //everywhere - identical to default mode.
        assertEquals(mapOf(listOf("h1") to 2, listOf("h3") to 1, listOf("h6") to 2), perSample.getValue("SampleZ"))
        assertEquals(emptyMap<List<String>, Int>(), perSample.getValue("SampleEmpty"))

        //Position/allele-matching stats are identical to default mode - only weighting differs.
        assertEquals(3, stats.positionsMatched)
        assertEquals(1, stats.positionsNotInIndex)
        assertEquals(6, stats.noCallAlleles)
        assertEquals(1, stats.allelesNotInIndex)

        assertEquals(15, stats.gameteHits)
        assertEquals(41L, stats.totalEvidenceUnits)
        assertEquals(12, stats.adFallbackCount)
        assertEquals(2, stats.zeroDepthAlleleCalls)
    }

    //endregion

    //region End-to-end

    @Test
    fun testEndToEndDefaultMode() {
        val outputDir = "${TestExtension.testOutputDir}defaultMode/"
        val result = BuildReadMappingFromVcf().test("--index-file $indexFile --impute-vcf $imputeVcf --output-dir $outputDir")
        assertEquals(0, result.statusCode)

        listOf("SampleX", "SampleY", "SampleZ", "SampleEmpty").forEach {
            assertTrue(File("$outputDir${it}_readMapping.txt").isFile, "missing readMapping for $it")
        }

        assertEquals(mapOf(listOf("h1") to 2, listOf("h3") to 1, listOf("h6") to 2), readMappingAsMap("${outputDir}SampleZ_readMapping.txt"))
        assertEquals(emptyMap<List<String>, Int>(), readMappingAsMap("${outputDir}SampleEmpty_readMapping.txt"))
    }

    @Test
    fun testEndToEndUseAlleleDepthMode() {
        val outputDir = "${TestExtension.testOutputDir}adMode/"
        val result = BuildReadMappingFromVcf().test("--index-file $indexFile --impute-vcf $imputeVcf --output-dir $outputDir --use-allele-depth")
        assertEquals(0, result.statusCode)

        assertEquals(
            mapOf(listOf("h1") to 12, listOf("h1", "h2") to 8, listOf("h3") to 2, listOf("h6") to 9, listOf("h4", "h5") to 1),
            readMappingAsMap("${outputDir}SampleX_readMapping.txt")
        )
        assertEquals(mapOf(listOf("h1", "h2") to 2, listOf("h3") to 2), readMappingAsMap("${outputDir}SampleY_readMapping.txt"))
    }

    @Test
    fun testEndToEndWritesPathKeyFile() {
        val outputDir = "${TestExtension.testOutputDir}pathKeyFileTest/"
        val result = BuildReadMappingFromVcf().test("--index-file $indexFile --impute-vcf $imputeVcf --output-dir $outputDir")
        assertEquals(0, result.statusCode)

        val lines = File("${outputDir}pathKeyFile.txt").readLines()
        assertEquals("sampleName\tfilename", lines[0])
        assertEquals(5, lines.size) //header + 4 samples

        val rows = lines.drop(1).map { it.split("\t") }.associate { it[0] to it[1] }
        assertEquals(setOf("SampleX", "SampleY", "SampleZ", "SampleEmpty"), rows.keys)
        rows.forEach { (sampleName, filename) ->
            assertEquals("${outputDir}${sampleName}_readMapping.txt", filename)
            assertTrue(File(filename).isFile)
        }
    }

    @Test
    fun testEndToEndCreatesOutputDirIfMissing() {
        val outputDir = "${TestExtension.testOutputDir}autoCreated/"
        assertTrue(!File(outputDir).exists())

        val result = BuildReadMappingFromVcf().test("--index-file $indexFile --impute-vcf $imputeVcf --output-dir $outputDir")

        assertEquals(0, result.statusCode)
        assertTrue(File(outputDir).isDirectory)
        assertTrue(File("${outputDir}SampleX_readMapping.txt").isFile)
    }

    //endregion

    //region Pipeline integration (steps 1 -> 2 compose)

    @Test
    fun testPipelineIntegrationWithStep1() {
        val builtIndexFile = "${TestExtension.testOutputDir}pipelineIndex.txt"
        val indexResult = BuildVcfHaplotypeIndex().test(
            "--hvcf-dir $step1HvcfDir --panel-vcf $step1PanelVcf --output-file $builtIndexFile"
        )
        assertEquals(0, indexResult.statusCode)

        val outputDir = "${TestExtension.testOutputDir}pipelineReadMappings/"
        val readMappingResult = BuildReadMappingFromVcf().test(
            "--index-file $builtIndexFile --impute-vcf $step1PanelVcf --output-dir $outputDir"
        )
        assertEquals(0, readMappingResult.statusCode)

        //testPanel.vcf's samples are TestA, TestB, Ghost - Ghost isn't in the step-1 hvcfs, so it
        //never contributed any rows to the index, but it's still a VCF sample here and must still
        //get a (likely empty) readMapping file.
        listOf("TestA", "TestB", "Ghost").forEach {
            assertTrue(File("$outputDir${it}_readMapping.txt").isFile, "missing readMapping for $it")
        }

        //Every count in TestA's readMapping should be a positive Int, and the total evidence
        //written should be internally consistent with a fresh recount straight from the files.
        val testAMap = readMappingAsMap("${outputDir}TestA_readMapping.txt")
        assertTrue(testAMap.isNotEmpty())
        testAMap.values.forEach { assertTrue(it > 0) }

        val (_, stats) = BuildReadMappingFromVcf().buildReadMappings(builtIndexFile, step1PanelVcf, useAlleleDepth = false)
        assertEquals(stats.gameteHits.toLong(), stats.totalEvidenceUnits)
    }

    //endregion
}
