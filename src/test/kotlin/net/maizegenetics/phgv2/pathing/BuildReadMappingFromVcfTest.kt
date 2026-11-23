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
import kotlin.test.assertFailsWith
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

        /**
         * Runs processRecord for a record known to be in the index and collects its evidence
         * into sampleName -> (hapIds -> count), the shape the old processRecord returned.
         */
        fun deltas(
            cmd: BuildReadMappingFromVcf,
            record: VariantContext,
            useAlleleDepth: Boolean,
            stats: ReadMappingStats
        ): Map<String, Map<List<String>, Int>> {
            val indexPosition = index().getValue(Position(record.contig, record.start))
            val result = mutableMapOf<String, MutableMap<List<String>, Int>>()
            cmd.processRecord(record, indexPosition, useAlleleDepth, stats) { sampleName, hapIds, weight ->
                val sampleMap = result.getOrPut(sampleName) { mutableMapOf() }
                sampleMap[hapIds] = (sampleMap[hapIds] ?: 0) + weight
            }
            return result
        }

        /**
         * Runs buildReadMappings into a fresh directory and reads every sample's file back.
         */
        fun buildAndRead(
            index: String,
            vcf: String,
            outDirName: String,
            useAlleleDepth: Boolean,
            maxBufferedEntries: Long = BuildReadMappingFromVcf.DEFAULT_MAX_BUFFERED_ENTRIES
        ): Pair<Map<String, Map<List<String>, Int>>, ReadMappingStats> {
            val outDir = "${TestExtension.testOutputDir}$outDirName"
            File(outDir).mkdirs()
            val (sampleNames, stats) = BuildReadMappingFromVcf().buildReadMappings(index, vcf, outDir, useAlleleDepth, maxBufferedEntries)
            return sampleNames.associateWith { readMappingAsMap("$outDir/${it}_readMapping.txt.gz") } to stats
        }
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
        assertEquals(0L, stats.adFallbackCount)
        assertEquals(0L, stats.zeroDepthAlleleCalls)
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
        assertEquals(0L, stats.adFallbackCount)
    }

    @Test
    fun testAlleleWeightFallsBackWhenADMissing() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(10)
        val genotype = record.getGenotype("SampleY") //no AD
        val stats = ReadMappingStats()

        val weight = cmd.alleleWeight(genotype, record.getAlternateAllele(0), record, useAlleleDepth = true, stats = stats)

        assertEquals(1, weight)
        assertEquals(1L, stats.adFallbackCount)
        assertEquals(0L, stats.zeroDepthAlleleCalls)
    }

    @Test
    fun testAlleleWeightFallsBackWhenAlleleIndexOutOfBounds() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(30)
        val genotype = record.getGenotype("SampleX") //AD=9, length 1 - only covers the REF allele
        val stats = ReadMappingStats()

        val refWeight = cmd.alleleWeight(genotype, record.reference, record, useAlleleDepth = true, stats = stats)
        assertEquals(9, refWeight)
        assertEquals(0L, stats.adFallbackCount)

        val altWeight = cmd.alleleWeight(genotype, record.getAlternateAllele(0), record, useAlleleDepth = true, stats = stats)
        assertEquals(1, altWeight)
        assertEquals(1L, stats.adFallbackCount)
    }

    @Test
    fun testAlleleWeightZeroDepthReturnsZero() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(30)
        val genotype = record.getGenotype("SampleY") //AD=6,0 - ALT depth is 0
        val stats = ReadMappingStats()

        val weight = cmd.alleleWeight(genotype, record.getAlternateAllele(0), record, useAlleleDepth = true, stats = stats)

        assertEquals(0, weight)
        assertEquals(1L, stats.zeroDepthAlleleCalls)
        assertEquals(0L, stats.adFallbackCount)
    }

    //endregion

    //region processRecord

    @Test
    fun testProcessRecordHetCallProducesTwoIndependentIncrements() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(10)
        val stats = ReadMappingStats()

        val deltas = deltas(cmd, record, useAlleleDepth = false, stats = stats)

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

        val deltas = deltas(cmd, record, useAlleleDepth = false, stats = stats)

        //SampleZ is 0/1 (C,G): C resolves to [h3], G has no entry in the index at 1:20.
        assertEquals(mapOf(listOf("h3") to 1), deltas.getValue("SampleZ"))
        assertEquals(1L, stats.allelesNotInIndex)
    }

    @Test
    fun testProcessRecordNoCallSkipped() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(10)
        val stats = ReadMappingStats()

        val deltas = deltas(cmd, record, useAlleleDepth = false, stats = stats)

        assertTrue(!deltas.containsKey("SampleEmpty"))
        assertEquals(2L, stats.noCallAlleles)
    }

    @Test
    fun testProcessRecordSharedAmbiguousHapIdSet() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(30)
        val stats = ReadMappingStats()

        val deltas = deltas(cmd, record, useAlleleDepth = false, stats = stats)

        //SampleZ is 0/0 (A,A) at 1:30 -> allele A resolves to the single hapId [h6].
        assertEquals(mapOf(listOf("h6") to 2), deltas.getValue("SampleZ"))
    }

    @Test
    fun testProcessRecordUseAlleleDepthWeighting() {
        val cmd = BuildReadMappingFromVcf()
        val record = imputeRecords().getValue(10)
        val stats = ReadMappingStats()

        val deltas = deltas(cmd, record, useAlleleDepth = true, stats = stats)

        //SampleX AD=12,8: allele A -> [h1] += 12, allele T -> [h1,h2] += 8 (not the flat +1/+1).
        assertEquals(mapOf(listOf("h1") to 12, listOf("h1", "h2") to 8), deltas.getValue("SampleX"))
    }

    //endregion

    //region buildReadMappings - default mode (hand-derived expected maps; see class KDoc)

    @Test
    fun testBuildReadMappingsDefaultMode() {
        val (perSample, stats) = buildAndRead(indexFile, imputeVcf, "buildDefault", useAlleleDepth = false)

        assertEquals(mapOf(listOf("h1") to 1, listOf("h1", "h2") to 1, listOf("h3") to 2, listOf("h6") to 1, listOf("h4", "h5") to 1), perSample.getValue("SampleX"))
        assertEquals(mapOf(listOf("h1", "h2") to 2, listOf("h3") to 2, listOf("h4", "h5") to 2), perSample.getValue("SampleY"))
        assertEquals(mapOf(listOf("h1") to 2, listOf("h3") to 1, listOf("h6") to 2), perSample.getValue("SampleZ"))
        //SampleEmpty has zero resolvable evidence but still gets an entry.
        assertEquals(emptyMap<List<String>, Int>(), perSample.getValue("SampleEmpty"))

        assertEquals(4L, stats.totalImputeRecords)
        assertEquals(3L, stats.positionsMatched)
        assertEquals(1L, stats.positionsNotInIndex)
        assertEquals(6L, stats.noCallAlleles)
        assertEquals(1L, stats.allelesNotInIndex)
        assertEquals(17L, stats.gameteHits)
        assertEquals(17L, stats.totalEvidenceUnits)
        assertEquals(0L, stats.adFallbackCount)
        assertEquals(0L, stats.zeroDepthAlleleCalls)
    }

    //endregion

    //region buildReadMappings - --use-allele-depth mode

    @Test
    fun testBuildReadMappingsUseAlleleDepthMode() {
        val (perSample, stats) = buildAndRead(indexFile, imputeVcf, "buildAD", useAlleleDepth = true)

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
        assertEquals(3L, stats.positionsMatched)
        assertEquals(1L, stats.positionsNotInIndex)
        assertEquals(6L, stats.noCallAlleles)
        assertEquals(1L, stats.allelesNotInIndex)

        assertEquals(15L, stats.gameteHits)
        assertEquals(41L, stats.totalEvidenceUnits)
        assertEquals(12L, stats.adFallbackCount)
        assertEquals(2L, stats.zeroDepthAlleleCalls)
    }

    //endregion

    //region End-to-end

    @Test
    fun testEndToEndDefaultMode() {
        val outputDir = "${TestExtension.testOutputDir}defaultMode/"
        val result = BuildReadMappingFromVcf().test("--index-file $indexFile --impute-vcf $imputeVcf --output-dir $outputDir")
        assertEquals(0, result.statusCode)

        listOf("SampleX", "SampleY", "SampleZ", "SampleEmpty").forEach {
            assertTrue(File("$outputDir${it}_readMapping.txt.gz").isFile, "missing readMapping for $it")
        }

        assertEquals(mapOf(listOf("h1") to 2, listOf("h3") to 1, listOf("h6") to 2), readMappingAsMap("${outputDir}SampleZ_readMapping.txt.gz"))
        assertEquals(emptyMap<List<String>, Int>(), readMappingAsMap("${outputDir}SampleEmpty_readMapping.txt.gz"))
    }

    @Test
    fun testEndToEndUseAlleleDepthMode() {
        val outputDir = "${TestExtension.testOutputDir}adMode/"
        val result = BuildReadMappingFromVcf().test("--index-file $indexFile --impute-vcf $imputeVcf --output-dir $outputDir --use-allele-depth")
        assertEquals(0, result.statusCode)

        assertEquals(
            mapOf(listOf("h1") to 12, listOf("h1", "h2") to 8, listOf("h3") to 2, listOf("h6") to 9, listOf("h4", "h5") to 1),
            readMappingAsMap("${outputDir}SampleX_readMapping.txt.gz")
        )
        assertEquals(mapOf(listOf("h1", "h2") to 2, listOf("h3") to 2), readMappingAsMap("${outputDir}SampleY_readMapping.txt.gz"))
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
            assertEquals("${outputDir}${sampleName}_readMapping.txt.gz", filename)
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
        assertTrue(File("${outputDir}SampleX_readMapping.txt.gz").isFile)
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
            assertTrue(File("$outputDir${it}_readMapping.txt.gz").isFile, "missing readMapping for $it")
        }

        //Every count in TestA's readMapping should be a positive Int, and the total evidence
        //written should be internally consistent with a fresh recount straight from the files.
        val testAMap = readMappingAsMap("${outputDir}TestA_readMapping.txt.gz")
        assertTrue(testAMap.isNotEmpty())
        testAMap.values.forEach { assertTrue(it > 0) }

        val (_, stats) = buildAndRead(builtIndexFile, step1PanelVcf, "pipelineRecount", useAlleleDepth = false)
        assertEquals(stats.gameteHits, stats.totalEvidenceUnits)
    }

    //endregion

    //region Streaming merge: contig order checks, multi-contig/multi-range, buffered flushes

    private val multiContigIndex = """
        #PHG_VCF_HAPLOTYPE_INDEX
        #version=1.0
        #chrom	pos	allele	refRange	hapIds
        chr2	10	A	chr2:1-100	a1
        chr2	10	T	chr2:1-100	a2,a3
        chr2	150	C	chr2:101-200	b1,b2
        chr2	150	G	chr2:101-200	b3
        chr3	5	A	chr3:1-100	c1
        chr10	20	G	chr10:1-100	d1,d2
        chr10	20	T	chr10:1-100	d3
    """.trimIndent() + "\n"

    private val vcfHeader = """
        ##fileformat=VCFv4.2
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
    """.trimIndent() + "\n"

    private fun writeTemp(name: String, content: String): String {
        val dir = File("${TestExtension.testOutputDir}streaming/").apply { mkdirs() }
        return File(dir, name).apply { writeText(content) }.path
    }

    @Test
    fun testMultiContigMultiRangeMergeSkipsContigsMissingFromEitherFile() {
        val index = writeTemp("multiIndex.txt", multiContigIndex)
        //chr3 is only in the index; chrUn is only in the VCF; chr2:50 is not in the index
        val vcf = writeTemp("multi.vcf", vcfHeader + """
            chr2	10	.	A	T	.	.	.	GT	0/1	1/1
            chr2	50	.	A	T	.	.	.	GT	0/1	1/1
            chr2	150	.	C	G	.	.	.	GT	0/0	1/1
            chrUn	7	.	A	T	.	.	.	GT	0/1	1/1
            chr10	20	.	G	T	.	.	.	GT	0/1	0/0
        """.trimIndent() + "\n")

        val (perSample, stats) = buildAndRead(index, vcf, "multiDefault", useAlleleDepth = false)

        assertEquals(
            mapOf(listOf("a1") to 1, listOf("a2", "a3") to 1, listOf("b1", "b2") to 2, listOf("d1", "d2") to 1, listOf("d3") to 1),
            perSample.getValue("S1")
        )
        assertEquals(mapOf(listOf("a2", "a3") to 2, listOf("b3") to 2, listOf("d1", "d2") to 2), perSample.getValue("S2"))
        assertEquals(3L, stats.positionsMatched)
        assertEquals(2L, stats.positionsNotInIndex)

        //Flushing after every range must give identical files to one big flush at the end.
        val (flushedEveryRange, flushedStats) = buildAndRead(index, vcf, "multiFlushEveryRange", useAlleleDepth = false, maxBufferedEntries = 1)
        assertEquals(perSample, flushedEveryRange)
        assertTrue(flushedStats.flushes > 1, "expected multiple flushes, got ${flushedStats.flushes}")
    }

    @Test
    fun testContigOrderMismatchFailsBeforeCounting() {
        val index = writeTemp("orderIndex.txt", multiContigIndex) //chr2, chr3, chr10 (natural order)
        //lexicographic order, as some tools sort: chr10 before chr2
        val vcf = writeTemp("lexOrder.vcf", vcfHeader + """
            chr10	20	.	G	T	.	.	.	GT	0/1	0/0
            chr2	10	.	A	T	.	.	.	GT	0/1	1/1
        """.trimIndent() + "\n")

        val outDir = "${TestExtension.testOutputDir}orderMismatch"
        val error = assertFailsWith<IllegalArgumentException> {
            BuildReadMappingFromVcf().buildReadMappings(index, vcf, outDir, useAlleleDepth = false)
        }
        assertTrue(error.message!!.contains("Contig order differs"), error.message)
        assertTrue(error.message!!.contains("index has 'chr2'"), error.message)
        //nothing should have been written - the check runs before any output
        assertTrue(!File("$outDir/S1_readMapping.txt.gz").exists())
    }

    @Test
    fun testNonContiguousContigInVcfFails() {
        val index = writeTemp("contigIndex.txt", multiContigIndex)
        val vcf = writeTemp("nonContiguous.vcf", vcfHeader + """
            chr2	10	.	A	T	.	.	.	GT	0/1	1/1
            chr10	20	.	G	T	.	.	.	GT	0/1	0/0
            chr2	150	.	C	G	.	.	.	GT	0/0	1/1
        """.trimIndent() + "\n")

        val error = assertFailsWith<IllegalArgumentException> {
            BuildReadMappingFromVcf().buildReadMappings(index, vcf, "${TestExtension.testOutputDir}nonContiguous", useAlleleDepth = false)
        }
        assertTrue(error.message!!.contains("non-contiguous"), error.message)
    }

    @Test
    fun testUnsortedPositionsInVcfFail() {
        val index = writeTemp("posIndex.txt", multiContigIndex)
        val vcf = writeTemp("unsortedPos.vcf", vcfHeader + """
            chr2	150	.	C	G	.	.	.	GT	0/0	1/1
            chr2	10	.	A	T	.	.	.	GT	0/1	1/1
        """.trimIndent() + "\n")

        val error = assertFailsWith<IllegalArgumentException> {
            BuildReadMappingFromVcf().buildReadMappings(index, vcf, "${TestExtension.testOutputDir}unsortedPos", useAlleleDepth = false)
        }
        assertTrue(error.message!!.contains("not sorted by position"), error.message)
    }

    @Test
    fun testScanContigOrderPlainAndGzipMatch() {
        val vcfText = vcfHeader + """
            chr2	10	.	A	T	.	.	.	GT	0/1	1/1
            chr2	150	.	C	G	.	.	.	GT	0/0	1/1
            chr10	20	.	G	T	.	.	.	GT	0/1	0/0
        """.trimIndent() + "\n"
        val plain = writeTemp("scan.vcf", vcfText)
        val gz = "${TestExtension.testOutputDir}streaming/scan.vcf.gz"
        java.util.zip.GZIPOutputStream(File(gz).outputStream()).bufferedWriter().use { it.write(vcfText) }

        assertEquals(listOf("chr2", "chr10"), VcfHaplotypeIndexUtils.scanContigOrder(plain))
        assertEquals(listOf("chr2", "chr10"), VcfHaplotypeIndexUtils.scanContigOrder(gz))
        assertEquals(listOf("chr2", "chr3", "chr10"), VcfHaplotypeIndexUtils.scanContigOrder(writeTemp("scanIndex.txt", multiContigIndex)))
    }

    @Test
    fun testRequireMatchingContigOrderIgnoresUnsharedContigs() {
        //contigs that only one side has don't affect the check
        VcfHaplotypeIndexUtils.requireMatchingContigOrder(
            listOf("1", "2", "3", "10"), listOf("0", "1", "3", "scaffold_9", "10"), "index", "vcf"
        )
        assertFailsWith<IllegalArgumentException> {
            VcfHaplotypeIndexUtils.requireMatchingContigOrder(listOf("1", "2", "10"), listOf("1", "10", "2"), "index", "vcf")
        }
    }

    //endregion

    //region Accumulator, file writer, CLI buffering option, index-side and duplicate-record cases

    @Test
    fun testAccumulatorInternsEqualSetsAndCountsEntries() {
        val accumulator = ReadMappingAccumulator(numSamples = 2)
        val setA = listOf("h1", "h2")
        val setAOtherInstance = listOf("h1", "h2") //equal content, different object
        val setB = listOf("h3")

        accumulator.add(0, setA, 1)
        accumulator.add(0, setAOtherInstance, 2) //same set as setA -> same entry
        accumulator.add(0, setB, 1)
        accumulator.add(1, setB, 5)

        assertEquals(3L, accumulator.bufferedEntries) //(S0,A), (S0,B), (S1,B)
        assertEquals(listOf(true, true), accumulator.sampleHasEvidence.toList())

        val outDir = "${TestExtension.testOutputDir}accumulator"
        val writer = ReadMappingFileWriter(outDir, listOf("S0", "S1"), "source.vcf")
        accumulator.flush(writer)
        assertEquals(0L, accumulator.bufferedEntries)

        assertEquals(mapOf(setA to 3, setB to 1), readMappingAsMap("$outDir/S0_readMapping.txt.gz"))
        assertEquals(mapOf(setB to 5), readMappingAsMap("$outDir/S1_readMapping.txt.gz"))

        //a second flush after more evidence appends, rather than overwriting
        accumulator.add(1, listOf("h7"), 4)
        accumulator.flush(writer)
        assertEquals(mapOf(setB to 5, listOf("h7") to 4), readMappingAsMap("$outDir/S1_readMapping.txt.gz"))

        //file really is gzip (magic bytes 1f 8b) and now holds 3 gzip members (header + 2 flushes),
        //all of which were read back above
        val bytes = File("$outDir/S1_readMapping.txt.gz").readBytes()
        assertEquals(listOf(0x1f.toByte(), 0x8b.toByte()), bytes.take(2))
        val memberCount = (0 until bytes.size - 2).count { bytes[it] == 0x1f.toByte() && bytes[it + 1] == 0x8b.toByte() && bytes[it + 2] == 8.toByte() }
        assertTrue(memberCount >= 3, "expected at least 3 gzip members, found $memberCount")
    }

    @Test
    fun testFileWriterCreatesDirAndHeaderForEverySample() {
        val outDir = "${TestExtension.testOutputDir}fileWriter/nested"
        assertTrue(!File(outDir).exists())
        ReadMappingFileWriter(outDir, listOf("A", "B"), "my.vcf")

        listOf("A", "B").forEach { sample ->
            val lines = net.maizegenetics.phgv2.utils.getBufferedReader("$outDir/${sample}_readMapping.txt.gz").readLines()
            assertEquals(listOf("#sampleName=$sample", "#filename1=my.vcf", "HapIds\tcount"), lines)
            assertEquals(emptyMap<List<String>, Int>(), readMappingAsMap("$outDir/${sample}_readMapping.txt.gz"))
        }
    }

    @Test
    fun testCliMaxBufferedEntriesOption() {
        val badValue = BuildReadMappingFromVcf().test(
            "--index-file $indexFile --impute-vcf $imputeVcf --output-dir ${TestExtension.testOutputDir}badBuffer --max-buffered-entries 0"
        )
        assertEquals(1, badValue.statusCode)
        assertTrue(badValue.stderr.contains("--max-buffered-entries must be positive"), badValue.stderr)

        val outputDir = "${TestExtension.testOutputDir}cliBuffer/"
        val result = BuildReadMappingFromVcf().test(
            "--index-file $indexFile --impute-vcf $imputeVcf --output-dir $outputDir --max-buffered-entries 1"
        )
        assertEquals(0, result.statusCode)
        assertEquals(mapOf(listOf("h1") to 2, listOf("h3") to 1, listOf("h6") to 2), readMappingAsMap("${outputDir}SampleZ_readMapping.txt.gz"))
    }

    @Test
    fun testUnsortedIndexFailsDuringMerge() {
        val index = writeTemp("unsortedIndex.txt", """
            #chrom	pos	allele	refRange	hapIds
            chr2	150	C	chr2:101-200	b1
            chr2	10	A	chr2:1-100	a1
        """.trimIndent() + "\n")
        val vcf = writeTemp("forUnsortedIndex.vcf", vcfHeader + "chr2\t150\t.\tC\tG\t.\t.\t.\tGT\t0/0\t0/0\n")

        val error = assertFailsWith<IllegalArgumentException> {
            BuildReadMappingFromVcf().buildReadMappings(index, vcf, "${TestExtension.testOutputDir}unsortedIndex", useAlleleDepth = false)
        }
        assertTrue(error.message!!.contains("not sorted by position"), error.message)
    }

    @Test
    fun testDuplicateVcfRecordsAtSamePositionBothCount() {
        val index = writeTemp("dupIndex.txt", multiContigIndex)
        val vcf = writeTemp("dupRecords.vcf", vcfHeader + """
            chr2	10	.	A	T	.	.	.	GT	0/0	1/1
            chr2	10	.	A	T	.	.	.	GT	0/0	1/1
            chr2	150	.	C	G	.	.	.	GT	0/0	0/0
        """.trimIndent() + "\n")

        val (perSample, stats) = buildAndRead(index, vcf, "dupRecords", useAlleleDepth = false)
        assertEquals(mapOf(listOf("a1") to 4, listOf("b1", "b2") to 2), perSample.getValue("S1"))
        assertEquals(mapOf(listOf("a2", "a3") to 4, listOf("b1", "b2") to 2), perSample.getValue("S2"))
        assertEquals(3L, stats.positionsMatched)
    }

    @Test
    fun testImputeVcfWithNoSharedContigsWritesEmptyFiles() {
        val index = writeTemp("noSharedIndex.txt", multiContigIndex)
        val vcf = writeTemp("noShared.vcf", vcfHeader + "1\t10\t.\tA\tT\t.\t.\t.\tGT\t0/1\t1/1\n")

        val (perSample, stats) = buildAndRead(index, vcf, "noShared", useAlleleDepth = false)
        assertEquals(setOf("S1", "S2"), perSample.keys)
        perSample.values.forEach { assertEquals(emptyMap<List<String>, Int>(), it) }
        assertEquals(0L, stats.positionsMatched)
        assertEquals(1L, stats.positionsNotInIndex)
        assertEquals(listOf("S1", "S2"), stats.samplesWithNoEvidence)
    }

    //endregion
}
