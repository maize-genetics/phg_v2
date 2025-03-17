package net.maizegenetics.phgv2.cli

import biokotlin.genome.AssemblyVariantInfo
import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.testing.test
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.utils.*
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.*

@ExtendWith(TestExtension::class)
class CreateMafVCFTest {
    companion object {
        @JvmStatic
        @AfterAll
        fun tearDown() {
            File(TestExtension.testVCFDir).deleteRecursively()
            File(TestExtension.testTileDBURI).deleteRecursively()
        }
    }

    @Test
    fun testCliktParams() {
        val createMAFVCF = CreateMafVcf()

        val resultMissingBed =
            createMAFVCF.test("--db-path ${TestExtension.testTileDBURI} --maf-dir ${TestExtension.testMafDir} --reference-file ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingBed.statusCode, 1)
        assertEquals(
            "Usage: create-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --bed\n", resultMissingBed.output
        )
        val resultMissingRef =
            createMAFVCF.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} --maf-dir ${TestExtension.testMafDir} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingRef.statusCode, 1)
        assertEquals(
            "Usage: create-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --reference-file\n", resultMissingRef.output
        )

        val resultMissingOutput =
            createMAFVCF.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} --maf-dir ${TestExtension.testMafDir} --reference-file ${TestExtension.testRefFasta}")
        assertEquals(resultMissingOutput.statusCode, 1)
        assertEquals(
            "Usage: create-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --output-dir\n", resultMissingOutput.output
        )

        val resultMissingMafDir =
            createMAFVCF.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} --reference-file ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingMafDir.statusCode, 1)
        assertEquals(
            "Usage: create-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --maf-dir\n", resultMissingMafDir.output
        )

    }


    @Test
    fun testLoadingBEDFile() {
        val bedFile = "data/test/buildMAFVCF/B73_Test.bed"

        val ranges = loadRanges(bedFile)

        assertEquals(4, ranges.size)

        assertEquals(Position("chr1", 1), ranges[0].first)
        assertEquals(Position("chr1", 40), ranges[0].second)

        assertEquals(Position("chr7", 15), ranges[1].first)
        assertEquals(Position("chr7", 48), ranges[1].second)
        assertEquals(Position("chr7", 451), ranges[2].first)
        assertEquals(Position("chr7", 456), ranges[2].second)

        assertEquals(Position("chr10", 1), ranges[3].first)
        assertEquals(Position("chr10", 40), ranges[3].second)

    }

    @Test
    fun testBuildRefSeq() {
        val testRef = "data/test/buildMAFVCF/B73_Test.fa"

        val createMafVcf = CreateMafVcf()
        val refSeq = createMafVcf.buildRefGenomeSeq(testRef)

        assertEquals(3, refSeq.size)
        //check seqSizes
        assertEquals(461, refSeq["chr7"]!!.size())
        assertEquals(40, refSeq["chr1"]!!.size())
        assertEquals(40, refSeq["chr10"]!!.size())

        //check seqs
        val truthChr7Seq = "AAAAAAAAAAAAAAAGGGAATGTTAACCAAATGAATTGTCTCTTACGGTGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAGATGGGT"

        assertEquals(truthChr7Seq, refSeq["chr7"]!!.seq())
        assertEquals("GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT", refSeq["chr1"]!!.seq())
        assertEquals("GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT", refSeq["chr10"]!!.seq())
    }

    @Test
    fun testSimpleBuildMAFVCF() {
        //Need to create the agc record before we run this:
        val fastaCreateFileNamesFile = "data/test/buildMAFVCF/fastaCreateFileNames.txt"
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/buildMAFVCF/B73_Test.fa"

        //Create the tileDB datasets - these are verified in AgcCompress()
        Initdb().createDataSets(TestExtension.testTileDBURI,"")

        val agcCompress = AgcCompress()
        // Create the initial compressed file
        val agcResult = agcCompress.test("--fasta-list ${fastaCreateFileNamesFile} --db-path ${dbPath} --reference-file ${refFasta}")
        println(agcResult.output)

        val createMAFVCF = CreateMafVcf()
        val result = createMAFVCF.test("--db-path ${dbPath} --bed data/test/buildMAFVCF/B73_Test.bed --reference-file ${refFasta} --maf-dir data/test/buildMAFVCF/mafs/ -o ${TestExtension.testVCFDir}")
        println(result.output)
        //compare the contents of the output gVCF files to the expected output
        compareTwoGVCFFiles("data/test/buildMAFVCF/truthGVCFs/B97_truth.g.vcf", "${TestExtension.testVCFDir}/B97_ASM_Test.g.vcf.gz")

        //test that metrics file was created
        assertTrue(File("${TestExtension.testVCFDir}/VCFMetrics.tsv").exists())

        //Now we need to compare the hVCF's sequence with the sequence coming from the MAF files to make sure things match correctly as well as the boundaries
        val outputHVCF = "${TestExtension.testVCFDir}/B97_ASM_Test.h.vcf.gz"
        val outputHVCFReader = VCFFileReader(File(outputHVCF),false)

        val outputHeader = outputHVCFReader.header.metaDataInInputOrder.filter { it.key == "ALT" }

        //pull out the hashes from the output header.
        val outputHashes = outputHeader.map { it.toString().split("<")[1].split(",").first().split("=").last() }.toSet()

        //Manually get the sequences out of the reference file
        val seqsManuallyExtractedFromMAF = setOf("AAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT",
            "AGGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATG",
            "TAAGGA"
            )

        val mafHashes = seqsManuallyExtractedFromMAF.map { getChecksumForString(it, "Md5") }.toSet()
        //compare the hashsets
        assertEquals(mafHashes,outputHashes)

        //Check the coordinates of the GVCF and make sure these match the hvcfs.
        //These were manually found by comparing BED and GVCF so they are right
        //Ref coord chr1	0	40 cooresponds to chr6 98 .. 142
        //chr7	14	48 corresponds to chr4 4245 .. 4281
        //chr7	450	456 corresponds to chr4 5247 .. 5252
        //chr10	0	40 corresponds to chr6 1098 .. 1142
        val truthCoords = setOf(Pair(Position("chr6", 98), Position("chr6", 142)),
            Pair(Position("chr4", 4245), Position("chr4", 4281)),
            Pair(Position("chr4", 5247), Position("chr4", 5252)),
            )

        val seqCoords = outputHeader.map {
            val tokens = it.toString().split("<")[1]
                .split(",")
                .map { token -> token.split("=") }.associate { token -> Pair(token[0], token[1]) }

            val region = tokens["Regions"]!!
            val regionSplit = region.split(":")
            val chr = regionSplit[0]
            val bounds = regionSplit[1].split("-")
            val start = bounds[0].toInt()
            val end = bounds[1].toInt()
            Pair(Position(chr, start), Position(chr, end))
        }.toSet()

        assertEquals(truthCoords, seqCoords)

        //Check the refCoords of the variants
        val bedCoords = setOf(Pair(Position("chr1", 1), Position("chr1", 40)),
            Pair(Position("chr7", 15), Position("chr7", 48)),
            Pair(Position("chr7", 451), Position("chr7", 456)),
            Pair(Position("chr10", 1), Position("chr10", 40)))
        val variants = outputHVCFReader.iterator().toList()
        val refCoords = variants.map { Pair(Position(it.contig,it.start),Position(it.contig, it.getAttributeAsInt("END",0))) }.toSet()

        assertEquals(bedCoords, refCoords)
    }

    @Test
    fun verifyRefMatchesASM() {
        val createMAFVCF = CreateMafVcf()
        val testFasta = "data/test/buildMAFVCF/B97_ASM_Test.fa"
        val testSeqs = createMAFVCF.buildRefGenomeSeq(testFasta)
        //check first record
        //NucSeq is zerobased inclusive inclusive, VCF is 1-based inclusive inclusive.
        val firstSeq = testSeqs["chr6"]!![97..141].seq()
        assertEquals("AAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT",firstSeq)


        //Looking for gvcfCoords [1098,1142] which should be [1097,1141] in the NucSeq
        val secondSeq = testSeqs["chr6"]!![1097..1141].seq()
        assertEquals("AAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT",secondSeq)

        //Looking for gvcfCoords [4243,4283]
        val firstPartOfChr7_asm4 = testSeqs["chr4"]!![4242 .. 4282].seq()
        assertEquals("AAAGGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG", firstPartOfChr7_asm4) //Need to add an A at start here because MAF block begins with deletion so we need to start at position 12

        //Looking for gvcfCoords[5247,5257]
        val secondPartOfChr7_asm4 = testSeqs["chr4"]!![5246 .. 5256].seq()
        assertEquals("TAAGGATCCCT",secondPartOfChr7_asm4)

    }

    @Test
    fun testBedRegionContainedInVariant() {
        val createMafVCF = CreateMafVcf()
        val variant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10),"+")

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(createMafVCF.bedRegionContainedInVariant(fullyContainedBed, variant))
        assertTrue(createMafVCF.bedRegionContainedInVariant(containedBed, variant))
        assertFalse(createMafVCF.bedRegionContainedInVariant(partiallyContainedStart, variant))
        assertFalse(createMafVCF.bedRegionContainedInVariant(partiallyContainedEnd, variant))
        assertFalse(createMafVCF.bedRegionContainedInVariant(notContained, variant))
        assertFalse(createMafVCF.bedRegionContainedInVariant(notContained2, variant))
    }

    @Test
    fun testBedRegionContainedInVariantInfo() {
        val createMafVCF = CreateMafVcf()
        val variantInfos = AssemblyVariantInfo("chr1",5,10,"T","A","T",true)

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(createMafVCF.bedRegionContainedInVariantInfo(fullyContainedBed, variantInfos))
        assertTrue(createMafVCF.bedRegionContainedInVariantInfo(containedBed, variantInfos))
        assertFalse(createMafVCF.bedRegionContainedInVariantInfo(partiallyContainedStart, variantInfos))
        assertFalse(createMafVCF.bedRegionContainedInVariantInfo(partiallyContainedEnd, variantInfos))
        assertFalse(createMafVCF.bedRegionContainedInVariantInfo(notContained, variantInfos))
        assertFalse(createMafVCF.bedRegionContainedInVariantInfo(notContained2, variantInfos))
    }

    @Test
    fun testVariantFullyContained() {
        val createMafVCF = CreateMafVcf()
        val variant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10),"+")

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertTrue(createMafVCF.variantFullyContained(fullyContainedBed, variant))
        assertFalse(createMafVCF.variantFullyContained(containedBed, variant))
        assertFalse(createMafVCF.variantFullyContained(partiallyContainedStart, variant))
        assertFalse(createMafVCF.variantFullyContained(partiallyContainedEnd, variant))
        assertFalse(createMafVCF.variantFullyContained(notContained, variant))
        assertFalse(createMafVCF.variantFullyContained(notContained2, variant))

    }

    @Test
    fun testVariantInfoFullyContained() {
        val createMafVCF = CreateMafVcf()
        val variantInfo = AssemblyVariantInfo("chr1",5,10,"T","A","T",true)

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertTrue(createMafVCF.variantInfoFullyContained(fullyContainedBed, variantInfo))
        assertFalse(createMafVCF.variantInfoFullyContained(containedBed, variantInfo))
        assertFalse(createMafVCF.variantInfoFullyContained(partiallyContainedStart, variantInfo))
        assertFalse(createMafVCF.variantInfoFullyContained(partiallyContainedEnd, variantInfo))
        assertFalse(createMafVCF.variantInfoFullyContained(notContained, variantInfo))
        assertFalse(createMafVCF.variantInfoFullyContained(notContained2, variantInfo))
    }

    @Test
    fun testVariantPartiallyContainedStart() {
        val createMafVCF = CreateMafVcf()
        val variant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10),"+")

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(createMafVCF.variantPartiallyContainedStart(fullyContainedBed, variant))
        assertFalse(createMafVCF.variantPartiallyContainedStart(containedBed, variant))
        assertTrue(createMafVCF.variantPartiallyContainedStart(partiallyContainedStart, variant))
        assertFalse(createMafVCF.variantPartiallyContainedStart(partiallyContainedEnd, variant))
        assertFalse(createMafVCF.variantPartiallyContainedStart(notContained, variant))
        assertFalse(createMafVCF.variantPartiallyContainedStart(notContained2, variant))
    }

    @Test
    fun testVariantInfoPartiallyContainedStart() {
        val createMafVCF = CreateMafVcf()
        val variantInfo = AssemblyVariantInfo("chr1",5,10,"T","A","T",true)

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(createMafVCF.variantInfoPartiallyContainedStart(fullyContainedBed, variantInfo))
        assertFalse(createMafVCF.variantInfoPartiallyContainedStart(containedBed, variantInfo))
        assertTrue(createMafVCF.variantInfoPartiallyContainedStart(partiallyContainedStart, variantInfo))
        assertFalse(createMafVCF.variantInfoPartiallyContainedStart(partiallyContainedEnd, variantInfo))
        assertFalse(createMafVCF.variantInfoPartiallyContainedStart(notContained, variantInfo))
        assertFalse(createMafVCF.variantInfoPartiallyContainedStart(notContained2, variantInfo))

    }

    @Test
    fun testVariantPartiallyContainedEnd() {
        val createMafVCF = CreateMafVcf()
        val variant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10),"+")

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(createMafVCF.variantPartiallyContainedEnd(fullyContainedBed, variant))
        assertFalse(createMafVCF.variantPartiallyContainedEnd(containedBed, variant))
        assertFalse(createMafVCF.variantPartiallyContainedEnd(partiallyContainedStart, variant))
        assertTrue(createMafVCF.variantPartiallyContainedEnd(partiallyContainedEnd, variant))
        assertFalse(createMafVCF.variantPartiallyContainedEnd(notContained, variant))
        assertFalse(createMafVCF.variantPartiallyContainedEnd(notContained2, variant))
    }

    @Test
    fun testVariantInfoPartiallyContainedEnd() {
        val createMafVcf = CreateMafVcf()
        val variantInfo = AssemblyVariantInfo("chr1",5,10,"T","A","T",true)

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(createMafVcf.variantInfoPartiallyContainedEnd(fullyContainedBed, variantInfo))
        assertFalse(createMafVcf.variantInfoPartiallyContainedEnd(containedBed, variantInfo))
        assertFalse(createMafVcf.variantInfoPartiallyContainedEnd(partiallyContainedStart, variantInfo))
        assertTrue(createMafVcf.variantInfoPartiallyContainedEnd(partiallyContainedEnd, variantInfo))
        assertFalse(createMafVcf.variantInfoPartiallyContainedEnd(notContained, variantInfo))
        assertFalse(createMafVcf.variantInfoPartiallyContainedEnd(notContained2, variantInfo))

    }

    @Test
    fun testVariantAfterRegion() {
        val createMafVCF = CreateMafVcf()
        val variant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10),"+")

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(createMafVCF.variantAfterRegion(fullyContainedBed, variant))
        assertFalse(createMafVCF.variantAfterRegion(containedBed, variant))
        assertFalse(createMafVCF.variantAfterRegion(partiallyContainedStart, variant))
        assertFalse(createMafVCF.variantAfterRegion(partiallyContainedEnd, variant))
        assertFalse(createMafVCF.variantAfterRegion(notContained2, variant))
        assertTrue(createMafVCF.variantAfterRegion(notContained, variant))
    }

    @Test
    fun testVariantInfoAfterRegion() {
        val createMafVcf = CreateMafVcf()
        val variantInfo = AssemblyVariantInfo("chr1",5,10,"T","A","T",true)

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(createMafVcf.variantInfoAfterRegion(fullyContainedBed, variantInfo))
        assertFalse(createMafVcf.variantInfoAfterRegion(containedBed, variantInfo))
        assertFalse(createMafVcf.variantInfoAfterRegion(partiallyContainedStart, variantInfo))
        assertFalse(createMafVcf.variantInfoAfterRegion(partiallyContainedEnd, variantInfo))
        assertFalse(createMafVcf.variantInfoAfterRegion(notContained2, variantInfo))
        assertTrue(createMafVcf.variantInfoAfterRegion(notContained, variantInfo))

    }

    @Test
    fun testIsVariantResizable() {
        val createMafVCF = CreateMafVcf()
        val refBlockVariant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10),"+")

        val multiAllelicSNPVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",7), Pair("AAA", "TTT"), Position("chr1",10), Position("chr1",12),"+")

        val standardSNPVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",5), Pair("A", "T"), Position("chr1",10), Position("chr1",10),"+")

        val insertionVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",5), Pair("A", "TTT"), Position("chr1",10), Position("chr1",12),"+")
        val deletionVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",7), Pair("AAA", "T"), Position("chr1",10), Position("chr1",10),"+")

        assertTrue(createMafVCF.isVariantResizable(refBlockVariant))
        assertTrue(createMafVCF.isVariantResizable(multiAllelicSNPVariant))
        assertTrue(createMafVCF.isVariantResizable(standardSNPVariant)) //This is technically true as single bp variants just return the bp when resized
        assertFalse(createMafVCF.isVariantResizable(insertionVariant))
        assertFalse(createMafVCF.isVariantResizable(deletionVariant))
    }

    @Test
    fun testIsVariantInfoResizable() {
        val createMafVCF = CreateMafVcf()
        val refBlockVariant = AssemblyVariantInfo("chr1",5,10,"T","A","T",true)

        val multiAllelicSNPVariant = AssemblyVariantInfo("chr1",5,7,"AAA","TTT","AAA",true)

        val standardSNPVariant = AssemblyVariantInfo("chr1",5,5,"T","A","T",true)

        val insertionVariant = AssemblyVariantInfo("chr1",5,5,"TTT","A","TTT",true)
        val deletionVariant = AssemblyVariantInfo("chr1",5,7,"T","AAA","T",true)

        assertTrue(createMafVCF.isVariantInfoResizable(refBlockVariant))
        assertTrue(createMafVCF.isVariantInfoResizable(multiAllelicSNPVariant))
        assertTrue(createMafVCF.isVariantInfoResizable(standardSNPVariant)) //This is technically true as single bp variants just return the bp when resized
        assertFalse(createMafVCF.isVariantInfoResizable(insertionVariant))
        assertFalse(createMafVCF.isVariantInfoResizable(deletionVariant))
    }

    /**
     * This unit test checks different cases of VariantContexts and the returned values depending on the requested ref
     * position to resize to.
     */
    @Ignore
    @Test
    fun testResizeVariantContext() {
        val createMafVCF = CreateMafVcf()
        //Testing refBlocks
        val refBlockVariant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",15), Position("chr1",20),"+")

        //Need to have reveresed asm coords because that is how it looks with GVCFs coming from Biokotlin
        val refBlockVariantNegativeStrand = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",20), Position("chr1",15),"+")

        assertEquals(15, createMafVCF.resizeVariantContext(refBlockVariant, 5, "+"))
        assertEquals(20, createMafVCF.resizeVariantContext(refBlockVariantNegativeStrand, 5, "-"))


        assertEquals(17, createMafVCF.resizeVariantContext(refBlockVariant, 7, "+"))
        assertEquals(18, createMafVCF.resizeVariantContext(refBlockVariantNegativeStrand, 7, "-"))

        assertEquals(15, createMafVCF.resizeVariantContext(refBlockVariant, 2, "+"))
        assertEquals(20, createMafVCF.resizeVariantContext(refBlockVariant, 100, "+"))

        //Negative strand flips so we resize to the correct start/end on ASM
        assertEquals(15, createMafVCF.resizeVariantContext(refBlockVariantNegativeStrand, 100, "-"))
        assertEquals(20, createMafVCF.resizeVariantContext(refBlockVariantNegativeStrand, 2, "-"))

        assertEquals(-1, createMafVCF.resizeVariantContext(refBlockVariant, 7, "NOT_A_STRAND"))

        //Testing MultiAllelicPolymorphisms
        val multiAllelicSNPVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",10),
                                                Pair("AAAAA", "TTTTT"), Position("chr1",10), Position("chr1",15),"+")
        val multiAllelicSNPVariantNegativeStrand = createSNPVC("B97", Position("chr1",5),Position("chr1",10),
                                                Pair("AAAAA", "TTTTT"), Position("chr1",15), Position("chr1",10),"+")
        //Check that we can resize the variant to the correct start/end on the ASM
        assertEquals(10, createMafVCF.resizeVariantContext(multiAllelicSNPVariant, 5, "+"))
        assertEquals(15, createMafVCF.resizeVariantContext(multiAllelicSNPVariantNegativeStrand, 5, "-"))
        assertEquals(12, createMafVCF.resizeVariantContext(multiAllelicSNPVariant, 7, "+"))
        assertEquals(13, createMafVCF.resizeVariantContext(multiAllelicSNPVariantNegativeStrand, 7, "-"))

        //Check for positions out of the record
        assertEquals(15, createMafVCF.resizeVariantContext(multiAllelicSNPVariant, 100, "+"))
        assertEquals(10, createMafVCF.resizeVariantContext(multiAllelicSNPVariantNegativeStrand, 100, "-"))
        assertEquals(10, createMafVCF.resizeVariantContext(multiAllelicSNPVariant, 2, "+"))
        assertEquals(15, createMafVCF.resizeVariantContext(multiAllelicSNPVariantNegativeStrand, 2, "-"))

        val standardSNPVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",5), Pair("A", "T"), Position("chr1",10), Position("chr1",10),"+")
        //no matter what we request it should return 10 as it is only one bp of size
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 5, "+"))
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 5, "-"))
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 7, "+"))
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 7, "-"))
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 100, "+"))
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 100, "-"))
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 2, "+"))
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 2, "-"))

        val insertionVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",5), Pair("A", "TTT"), Position("chr1",10), Position("chr1",12),"+")
        //Everything should return -1 as its not resizable
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 5, "+"))
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 5, "-"))
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 7, "+"))
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 7, "-"))
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 100, "+"))
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 100, "-"))
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 2, "+"))
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 2, "-"))

        val deletionVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",7), Pair("AAA", "T"), Position("chr1",10), Position("chr1",10),"+")
        //Everything should return -1 as its not resizeable
        assertEquals(-1, createMafVCF.resizeVariantContext(deletionVariant, 5, "+"))
        assertEquals(-1, createMafVCF.resizeVariantContext(deletionVariant, 5, "-"))
        assertEquals(-1, createMafVCF.resizeVariantContext(deletionVariant, 7, "+"))
        assertEquals(-1, createMafVCF.resizeVariantContext(deletionVariant, 7, "-"))
        assertEquals(-1, createMafVCF.resizeVariantContext(deletionVariant, 100, "+"))
        assertEquals(-1, createMafVCF.resizeVariantContext(deletionVariant, 100, "-"))
        assertEquals(-1, createMafVCF.resizeVariantContext(deletionVariant, 2, "+"))
        assertEquals(-1, createMafVCF.resizeVariantContext(deletionVariant, 2, "-"))
    }

    @Test
    fun testResizeVariantInfo() {
        val createMafVCF = CreateMafVcf()
        //Testing refBlocks
        val refBlockVariant = AssemblyVariantInfo("chr1",5,10,"T","A","T",false, intArrayOf(),"chr1", 15, 20)

        //Need to have reversed asm coords because that is how it looks with GVCFs coming from Biokotlin
        val refBlockVariantNegativeStrand = AssemblyVariantInfo("chr1",5,10,"T","A","T",false, intArrayOf(),"chr1", 20, 15)

        assertEquals(15, createMafVCF.resizeVariantInfo(refBlockVariant, 5, "+"))
        assertEquals(20, createMafVCF.resizeVariantInfo(refBlockVariantNegativeStrand, 5, "-"))

        assertEquals(17, createMafVCF.resizeVariantInfo(refBlockVariant, 7, "+"))
        assertEquals(18, createMafVCF.resizeVariantInfo(refBlockVariantNegativeStrand, 7, "-"))

        assertEquals(15, createMafVCF.resizeVariantInfo(refBlockVariant, 2, "+"))
        assertEquals(20, createMafVCF.resizeVariantInfo(refBlockVariant, 100, "+"))

        //Negative strand flips so we resize to the correct start/end on ASM
        assertEquals(15, createMafVCF.resizeVariantInfo(refBlockVariantNegativeStrand, 100, "-"))
        assertEquals(20, createMafVCF.resizeVariantInfo(refBlockVariantNegativeStrand, 2, "-"))

        assertEquals(-1, createMafVCF.resizeVariantInfo(refBlockVariant, 7, "NOT_A_STRAND"))


        //Testing MultiAllelicPolymorphisms
        val multiAllelicSNPVariant = AssemblyVariantInfo("chr1",5,10,"AAAAA","TTTTT","AAAAA",true, intArrayOf(),"chr1", 10, 15)
        val multiAllelicSNPVariantNegativeStrand = AssemblyVariantInfo("chr1",5,10,"AAAAA","TTTTT","AAAAA",true, intArrayOf(),"chr1", 15, 10)
        //Check that we can resize the variant to the correct start/end on the ASM
        assertEquals(10, createMafVCF.resizeVariantInfo(multiAllelicSNPVariant, 5, "+"))
        assertEquals(15, createMafVCF.resizeVariantInfo(multiAllelicSNPVariantNegativeStrand, 5, "-"))

        assertEquals(12, createMafVCF.resizeVariantInfo(multiAllelicSNPVariant, 7, "+"))
        assertEquals(13, createMafVCF.resizeVariantInfo(multiAllelicSNPVariantNegativeStrand, 7, "-"))

        //Check for positions out of the record
        assertEquals(15, createMafVCF.resizeVariantInfo(multiAllelicSNPVariant, 100, "+"))
        assertEquals(10, createMafVCF.resizeVariantInfo(multiAllelicSNPVariantNegativeStrand, 100, "-"))
        assertEquals(10, createMafVCF.resizeVariantInfo(multiAllelicSNPVariant, 2, "+"))
        assertEquals(15, createMafVCF.resizeVariantInfo(multiAllelicSNPVariantNegativeStrand, 2, "-"))

        val standardSNPVariant = AssemblyVariantInfo("chr1",5,5,"T","A","T",true, intArrayOf(),"chr1", 10, 10)
        //no matter what we request it should return 10 as it is only one bp of size
        assertEquals(10, createMafVCF.resizeVariantInfo(standardSNPVariant, 5, "+"))
        assertEquals(10, createMafVCF.resizeVariantInfo(standardSNPVariant, 5, "-"))
        assertEquals(10, createMafVCF.resizeVariantInfo(standardSNPVariant, 7, "+"))
        assertEquals(10, createMafVCF.resizeVariantInfo(standardSNPVariant, 7, "-"))
        assertEquals(10, createMafVCF.resizeVariantInfo(standardSNPVariant, 100, "+"))
        assertEquals(10, createMafVCF.resizeVariantInfo(standardSNPVariant, 100, "-"))
        assertEquals(10, createMafVCF.resizeVariantInfo(standardSNPVariant, 2, "+"))
        assertEquals(10, createMafVCF.resizeVariantInfo(standardSNPVariant, 2, "-"))

        val insertionVariant = AssemblyVariantInfo("chr1",5,5,"TTT","A","TTT",true, intArrayOf(),"chr1", 10, 12)
        //Everything should return -1 as its not resizable
        assertEquals(-1, createMafVCF.resizeVariantInfo(insertionVariant, 5, "+"))
        assertEquals(-1, createMafVCF.resizeVariantInfo(insertionVariant, 5, "-"))
        assertEquals(-1, createMafVCF.resizeVariantInfo(insertionVariant, 7, "+"))
        assertEquals(-1, createMafVCF.resizeVariantInfo(insertionVariant, 7, "-"))
        assertEquals(-1, createMafVCF.resizeVariantInfo(insertionVariant, 100, "+"))
        assertEquals(-1, createMafVCF.resizeVariantInfo(insertionVariant, 100, "-"))
        assertEquals(-1, createMafVCF.resizeVariantInfo(insertionVariant, 2, "+"))
        assertEquals(-1, createMafVCF.resizeVariantInfo(insertionVariant, 2, "-"))

        //val deletionVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",7), Pair("AAA", "T"), Position("chr1",10), Position("chr1",10),"+")
        val deletionVariant = AssemblyVariantInfo("chr1",5,7,"T","AAA","T",true, intArrayOf(),"chr1", 10, 10)
        //Everything should return -1 as its not resizeable
        assertEquals(-1, createMafVCF.resizeVariantInfo(deletionVariant, 5, "+"))
        assertEquals(-1, createMafVCF.resizeVariantInfo(deletionVariant, 5, "-"))
        assertEquals(-1, createMafVCF.resizeVariantInfo(deletionVariant, 7, "+"))
        assertEquals(-1, createMafVCF.resizeVariantInfo(deletionVariant, 7, "-"))
        assertEquals(-1, createMafVCF.resizeVariantInfo(deletionVariant, 100, "+"))
        assertEquals(-1, createMafVCF.resizeVariantInfo(deletionVariant, 100, "-"))
        assertEquals(-1, createMafVCF.resizeVariantInfo(deletionVariant, 2, "+"))
        assertEquals(-1, createMafVCF.resizeVariantInfo(deletionVariant, 2, "-"))


    }


    @Test
    fun testBuildRegionStrings() {
        val createMafVcf = CreateMafVcf()


        //Create a list of variantContext records with refBlocks and SNPs to resize and create boundaries
        val variantContexts = listOf(
            createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
                    Position("chr1",5), Position("chr1",10),
                    Position("chr1",5), Position("chr1",10),"+"),
            createSNPVC("B97", Position("chr1",11),Position("chr1",11), Pair("A", "T"),
                Position("chr1",11), Position("chr1",11),"+"),
            createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
                    Position("chr1",12), Position("chr1",15),
                    Position("chr1",12), Position("chr1",15),"+"),
            createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
                    Position("chr1",19), Position("chr1",20),
                    Position("chr1",19), Position("chr1",20),"+"),
            createSNPVC("B97", Position("chr1",21),Position("chr1",21), Pair("A", "T"),
                    Position("chr1",21), Position("chr1",21),"+"),
            createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
                    Position("chr1",22), Position("chr1",25),
                    Position("chr1",22), Position("chr1",25),"+"),
            )

        val regionStrings = createMafVcf.buildNewAssemblyRegions(7,23,variantContexts)
            .map { "${it.first.contig}:${it.first.position}-${it.second.position}" }

        //List<Pair<Position,Position>>
        assertEquals(2, regionStrings.size)
        assertEquals("chr1:7-15", regionStrings[0])
        assertEquals("chr1:19-23", regionStrings[1])

        // test condition where one VariantContext spans the entire reference range

        val oneVariantContext = listOf(createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",60), Position("chr1",300),
            Position("chr1",160), Position("chr1",400),"+"))

        val containedRegionStrings = createMafVcf.buildNewAssemblyRegions(180, 350, oneVariantContext)
            .map { "${it.first.contig}:${it.first.position}-${it.second.position}" }

        assertEquals(1, containedRegionStrings.size)
        assertEquals("chr1:180-350", containedRegionStrings[0])


    }

    @Test
    fun testMergeConsecutiveRegions() {
        val createMafVcf = CreateMafVcf()

        val regions = mutableListOf(Pair(Position("chr1",1),Position("chr1",10)),
            Pair(Position("chr1",11),Position("chr1",20)),
            Pair(Position("chr1",21),Position("chr1",30)),
            Pair(Position("chr1",31),Position("chr1",40)),
            Pair(Position("chr1",41),Position("chr1",50)),
            Pair(Position("chr1",51),Position("chr1",60)),
            Pair(Position("chr1",61),Position("chr1",70)),
            Pair(Position("chr1",71),Position("chr1",80)),
            Pair(Position("chr1",81),Position("chr1",90)),
            Pair(Position("chr1",91),Position("chr1",100)),
            Pair(Position("chr1",201),Position("chr1",210)),
            Pair(Position("chr1",211),Position("chr1",220)),
            Pair(Position("chr1",221),Position("chr1",230)),
            Pair(Position("chr1",231),Position("chr1",240)),
            Pair(Position("chr1",241),Position("chr1",250)),
            Pair(Position("chr1",251),Position("chr1",260)),
            Pair(Position("chr1",261),Position("chr1",270)),
            Pair(Position("chr1",271),Position("chr1",280)))

        val mergedRegions = createMafVcf.mergeConsecutiveRegions(regions)

        assertEquals(2, mergedRegions.size)
        assertEquals(Pair(Position("chr1",1),Position("chr1",100)), mergedRegions[0])
        assertEquals(Pair(Position("chr1",201),Position("chr1",280)), mergedRegions[1])


        //make inverted regions
        val invertedRegions = mutableListOf(Pair(Position("chr1",280),Position("chr1",271)),
            Pair(Position("chr1",270),Position("chr1",261)),
            Pair(Position("chr1",260),Position("chr1",251)),
            Pair(Position("chr1",250),Position("chr1",241)),
            Pair(Position("chr1",240),Position("chr1",231)),
            Pair(Position("chr1",230),Position("chr1",221)),
            Pair(Position("chr1",220),Position("chr1",211)),
            Pair(Position("chr1",210),Position("chr1",201)),
            Pair(Position("chr1",100),Position("chr1",91)),
            Pair(Position("chr1",90),Position("chr1",81)),
            Pair(Position("chr1",80),Position("chr1",71)),
            Pair(Position("chr1",70),Position("chr1",61)),
            Pair(Position("chr1",60),Position("chr1",51)),
            Pair(Position("chr1",50),Position("chr1",41)),
            Pair(Position("chr1",40),Position("chr1",31)),
            Pair(Position("chr1",30),Position("chr1",21)),
            Pair(Position("chr1",20),Position("chr1",11)),
            Pair(Position("chr1",10),Position("chr1",1)))

        val mergedInvertedRegions = createMafVcf.mergeConsecutiveRegions(invertedRegions)
        assertEquals(2, mergedInvertedRegions.size)
        assertEquals(Pair(Position("chr1",280),Position("chr1",201)), mergedInvertedRegions[0])
        assertEquals(Pair(Position("chr1",100),Position("chr1",1)), mergedInvertedRegions[1])


        //test mixed positive then inverted regions
        val mixedRegions = mutableListOf(Pair(Position("chr1",1),Position("chr1",10)),
            Pair(Position("chr1",11),Position("chr1",20)),
            Pair(Position("chr1",21),Position("chr1",30)),
            Pair(Position("chr1",31),Position("chr1",40)),
            Pair(Position("chr1",41),Position("chr1",50)),
            Pair(Position("chr1",51),Position("chr1",60)),
            Pair(Position("chr1",61),Position("chr1",70)),
            Pair(Position("chr1",71),Position("chr1",80)),
            Pair(Position("chr1",81),Position("chr1",90)),
            Pair(Position("chr1",91),Position("chr1",100)),
            Pair(Position("chr1",280),Position("chr1",271)),
            Pair(Position("chr1",270),Position("chr1",261)),
            Pair(Position("chr1",260),Position("chr1",251)),
            Pair(Position("chr1",250),Position("chr1",241)),
            Pair(Position("chr1",240),Position("chr1",231)),
            Pair(Position("chr1",230),Position("chr1",221)),
            Pair(Position("chr1",220),Position("chr1",211)),
            Pair(Position("chr1",210),Position("chr1",201)))

        val mergedMixedRegions = createMafVcf.mergeConsecutiveRegions(mixedRegions)
        assertEquals(2, mergedMixedRegions.size)
        assertEquals(Pair(Position("chr1",1),Position("chr1",100)), mergedMixedRegions[0])
        assertEquals(Pair(Position("chr1",280),Position("chr1",201)), mergedMixedRegions[1])


        //test mixed inverted then positive regions
        val mixedRegions2 = mutableListOf(Pair(Position("chr1",280),Position("chr1",271)),
            Pair(Position("chr1",270),Position("chr1",261)),
            Pair(Position("chr1",260),Position("chr1",251)),
            Pair(Position("chr1",250),Position("chr1",241)),
            Pair(Position("chr1",240),Position("chr1",231)),
            Pair(Position("chr1",230),Position("chr1",221)),
            Pair(Position("chr1",220),Position("chr1",211)),
            Pair(Position("chr1",210),Position("chr1",201)),
            Pair(Position("chr1",1),Position("chr1",10)),
            Pair(Position("chr1",11),Position("chr1",20)),
            Pair(Position("chr1",21),Position("chr1",30)),
            Pair(Position("chr1",31),Position("chr1",40)),
            Pair(Position("chr1",41),Position("chr1",50)),
            Pair(Position("chr1",51),Position("chr1",60)),
            Pair(Position("chr1",61),Position("chr1",70)),
            Pair(Position("chr1",71),Position("chr1",80)),
            Pair(Position("chr1",81),Position("chr1",90)),
            Pair(Position("chr1",91),Position("chr1",100)))

        val mergedMixedRegions2 = createMafVcf.mergeConsecutiveRegions(mixedRegions2)
        assertEquals(2, mergedMixedRegions2.size)
        assertEquals(Pair(Position("chr1",280),Position("chr1",201)), mergedMixedRegions2[0])
        assertEquals(Pair(Position("chr1",1),Position("chr1",100)), mergedMixedRegions2[1])


        //Test merging when they are consecutive but switch strands
        val mixedRegions3 = mutableListOf(Pair(Position("chr1",1),Position("chr1",10)),
            Pair(Position("chr1",11),Position("chr1",20)),
            Pair(Position("chr1",21),Position("chr1",30)),
            Pair(Position("chr1",31),Position("chr1",40)),
            Pair(Position("chr1",41),Position("chr1",50)),
            Pair(Position("chr1",51),Position("chr1",60)),
            Pair(Position("chr1",61),Position("chr1",70)),
            Pair(Position("chr1",71),Position("chr1",80)),
            Pair(Position("chr1",81),Position("chr1",90)),
            Pair(Position("chr1",91),Position("chr1",100)),
            Pair(Position("chr1",150),Position("chr1",141)),
            Pair(Position("chr1",140),Position("chr1",121)),
            Pair(Position("chr1",120),Position("chr1",101)))
        val mergedMixedRegions3 = createMafVcf.mergeConsecutiveRegions(mixedRegions3)
        assertEquals(2, mergedMixedRegions3.size)
        assertEquals(Pair(Position("chr1",1),Position("chr1",100)), mergedMixedRegions3[0])
        assertEquals(Pair(Position("chr1",150),Position("chr1",101)), mergedMixedRegions3[1])


        //have inverted regions that are single bp regions
        val singleRegions = mutableListOf(Pair(Position("chr1",1),Position("chr1",30)),
                Pair(Position("chr1",31),Position("chr1",31)),
                Pair(Position("chr1",32),Position("chr1",60)))

        val mergedSingleRegions = createMafVcf.mergeConsecutiveRegions(singleRegions)
        assertEquals(1, mergedSingleRegions.size)
        assertEquals(Pair(Position("chr1",1),Position("chr1",60)), mergedSingleRegions[0])

        val singleRegionsInverted = mutableListOf(Pair(Position("chr1",60),Position("chr1",32)),
                Pair(Position("chr1",31),Position("chr1",31)),
                Pair(Position("chr1",30),Position("chr1",1)))

        val mergedSingleRegionsInverted = createMafVcf.mergeConsecutiveRegions(singleRegionsInverted)
        assertEquals(1, mergedSingleRegionsInverted.size)
        assertEquals(Pair(Position("chr1",60),Position("chr1",1)), mergedSingleRegionsInverted[0])
    }

    @Test
    fun testConvertVariantContextToPositionRange() {
        val createMafVcf = CreateMafVcf()

        val refBlockVariant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",15), Position("chr1",20),"+")

        val range = createMafVcf.convertVariantContextToPositionRange(refBlockVariant)

        assertEquals(Position("chr1",15), range.first)
        assertEquals(Position("chr1",20), range.second)


        //test a reverse complimented variant
        val refBlockVariant2 = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",20), Position("chr1",15),"+")

        val range2 = createMafVcf.convertVariantContextToPositionRange(refBlockVariant2)

        assertEquals(Position("chr1",20), range2.first)
        assertEquals(Position("chr1",15), range2.second)

    }

    @Test
    fun testResizePositionRange() {
        val createMafVcf = CreateMafVcf()

        val range = Pair(Position("chr1",5), Position("chr1",20))
        val resizedRange = createMafVcf.resizePositionRange(range, 15, true)

        assertEquals(Position("chr1",15), resizedRange.first)
        assertEquals(Position("chr1",20), resizedRange.second)

        val resizedRange2 = createMafVcf.resizePositionRange(range, 15, false)

        assertEquals(Position("chr1",5), resizedRange2.first)
        assertEquals(Position("chr1",15), resizedRange2.second)

    }

    @Test
    fun testGVCFResizingIssue() {
        //Issue found with one of Matt's cassava lines
        //I think the issue is that we have a refBlock variant which is before the end of the bed region which then attempts to get resized.
        //The resulting coordiantes fall within the negative boundaries which is not correct.

        // |-------BED------------|
        //           |--GCVF--|

        val createMafVcf = CreateMafVcf()

        //create a bedRange
        //Chromosome04	35928654	35931063
        val bedRange = Pair(Position("Chromosome04",35928654), Position("Chromosome04",35931063))


        //Test forward strand - working initially

        val refBlockVariantForward = createRefRangeVC(mapOf("Chromosome04" to NucSeq("A".repeat(35929999))),"cassava_test",
            Position("Chromosome04",35920437), Position("Chromosome04",35929999),
            Position("chr8",4514), Position("chr8",14076),"+")

        val resizedVariantStartForward = createMafVcf.resizeVariantContext(refBlockVariantForward, bedRange.first.position, "+")
        assertEquals(12731, resizedVariantStartForward)

        val resizedVariantEndForward = createMafVcf.resizeVariantContext(refBlockVariantForward, bedRange.second.position, "+")
        assertEquals(14076, resizedVariantEndForward)

        //Make a variant context based on this:
//        Chromosome04	35920437	.	C	<NON_REF>	.	.	ASM_Chr=chr8;ASM_End=4514;ASM_Start=14076;ASM_Strand=-;END=35929999	GT:AD:DP:PL	0:30,0:30:0,90,90

        //Need to have reversed asm coords because that is how it looks with GVCFs coming from Biokotlin
        val refBlockVariant = createRefRangeVC(mapOf("Chromosome04" to NucSeq("A".repeat(35929999))),"cassava_test",
            Position("Chromosome04",35920437), Position("Chromosome04",35929999),
            Position("chr8",14076), Position("chr8",4514),"+")


        val resizedVariantStart = createMafVcf.resizeVariantContext(refBlockVariant, bedRange.first.position, "-")
        assertEquals(5859, resizedVariantStart)

        val resizedVariantEnd = createMafVcf.resizeVariantContext(refBlockVariant, bedRange.second.position, "-")
        assertEquals(4514, resizedVariantEnd)

    }

    /**
     * Function to compare the output gVCF file with the expected gVCF
     * It compares the alleles, depths, genotypes and the ASM metadata.
     */
    fun compareTwoGVCFFiles(truthGVCFFile : String, generatedFile: String) {
        //Load in the output GVCF  and the truth GVCF and verify that the output is correct
        val truthVariantIterator = VCFFileReader(File(truthGVCFFile),false).iterator()
        val truthVariants = mutableListOf<VariantContext>()
        while(truthVariantIterator.hasNext()) {
            truthVariants.add(truthVariantIterator.next())
        }
        val truthMap = truthVariants.associateBy { Position(it.contig, it.start) }

        val outputVariantIterator = VCFFileReader(File(generatedFile), false).iterator()
        val outputVariants = mutableListOf<VariantContext>()
        while(outputVariantIterator.hasNext()) {
            outputVariants.add(outputVariantIterator.next())
        }

        assertEquals(truthVariants.size, outputVariants.size,"Number of Variants does not match:")

        for(variant in outputVariants) {
            if(!truthMap.containsKey(Position(variant.contig, variant.start))) {
                fail("No matching variant found: ${variant.contig}:${variant.start}")
            }
            val matchingTruth = truthMap[Position(variant.contig, variant.start)]!!

            //Check END
            assertEquals(matchingTruth.end, variant.end,
                "End position does not match: outputFile: ${variant.contig}:${variant.start}-${variant.end}, " +
                        "Truth: ${matchingTruth.contig}:${matchingTruth.start}-${matchingTruth.end}")

            //Check alleles
            for(element in matchingTruth.alleles.zip(variant.alleles)) {
                assertEquals(element.first, element.second, "Alleles do not match: ${variant.contig}:${variant.start}-${variant.end}")
            }
            //Check GT
            assertEquals( matchingTruth.getGenotype(0).genotypeString,variant.getGenotype(0).genotypeString,
                "GT Fields do not match: ${variant.contig}:${variant.start}-${variant.end}")
            //Check AD
            for(element in matchingTruth.getGenotype(0).ad.zip(variant.getGenotype(0).ad)) {
                assertEquals(element.first, element.second, "AD fields do not match: ${variant.contig}:${variant.start}-${variant.end}")
            }
            //Check ASM Contig
            assertEquals(matchingTruth.getAttribute("ASM_Chr"), variant.getAttribute("ASM_Chr"),
                "ASM_Contig does not match: ${variant.contig}:${variant.start}-${variant.end}")
            //Check ASM Start
            assertEquals(matchingTruth.getAttribute("ASM_Start"), variant.getAttribute("ASM_Start"),
                "ASM_Start does not match: ${variant.contig}:${variant.start}-${variant.end}")
            //Check ASM END
            assertEquals(matchingTruth.getAttribute("ASM_End"), variant.getAttribute("ASM_End"),
                "ASM_End does not match: ${variant.contig}:${variant.start}-${variant.end}")
            //Check ASM Strand
            assertEquals(matchingTruth.getAttribute("ASM_Strand"), variant.getAttribute("ASM_Strand"),
                "ASM_Strand does not match: ${variant.contig}:${variant.start}-${variant.end}")

        }
    }

}