package net.maizegenetics.phgv2.cli

import biokotlin.seq.NucSeq
import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.testing.test
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.createRefRangeVC
import net.maizegenetics.phgv2.utils.createSNPVC
import net.maizegenetics.phgv2.utils.getChecksumForString
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
            createMAFVCF.test("--db-path ${TestExtension.testTileDBURI} --maf-dir ${TestExtension.testMafDir} --reference ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingBed.statusCode, 1)
        assertEquals(
            "Usage: create-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --bed: --bed must not be blank\n", resultMissingBed.output
        )
        val resultMissingRef =
            createMAFVCF.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} --maf-dir ${TestExtension.testMafDir} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingRef.statusCode, 1)
        assertEquals(
            "Usage: create-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --reference: --reference must not be blank\n", resultMissingRef.output
        )

        val resultMissingOutput =
            createMAFVCF.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} --maf-dir ${TestExtension.testMafDir} --reference ${TestExtension.testRefFasta}")
        assertEquals(resultMissingOutput.statusCode, 1)
        assertEquals(
            "Usage: create-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --output-dir: --output-dir/-o must not be blank\n", resultMissingOutput.output
        )

        val resultMissingMafDir =
            createMAFVCF.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} --reference ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingMafDir.statusCode, 1)
        assertEquals(
            "Usage: create-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --maf-dir: --maf-dir must not be blank\n", resultMissingMafDir.output
        )

        val resultMissingDbPath =
            createMAFVCF.test("--bed ${TestExtension.testBEDFile} --maf-dir ${TestExtension.testMafDir} --reference ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingDbPath.statusCode, 1)
        assertEquals(
            "Usage: create-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --db-path: --db-path must not be blank\n", resultMissingDbPath.output
        )

    }


    @Test
    fun testLoadingBEDFile() {
        val bedFile = "data/test/buildMAFVCF/B73_Test.bed"

        val createMAFVCF = CreateMafVcf()
        val ranges = createMAFVCF.loadRanges(bedFile)

        assertEquals(4, ranges.size)

        assertEquals(Position("chr1", 1), ranges[0].first)
        assertEquals(Position("chr1", 40), ranges[0].second)
        assertEquals(Position("chr10", 1), ranges[1].first)
        assertEquals(Position("chr10", 40), ranges[1].second)
        assertEquals(Position("chr7", 15), ranges[2].first)
        assertEquals(Position("chr7", 48), ranges[2].second)
        assertEquals(Position("chr7", 451), ranges[3].first)
        assertEquals(Position("chr7", 456), ranges[3].second)

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


        val agcCompress = AgcCompress()
        // Create the initial compressed file
        agcCompress.test("--fasta-list ${fastaCreateFileNamesFile} --db-path ${dbPath} --ref-fasta ${refFasta}")


        val createMAFVCF = CreateMafVcf()
//        createMAFVCF.test("--db-path data/test/buildMAFVCF/B97_ASM_Test.fa --bed data/test/buildMAFVCF/B73_Test.bed --reference data/test/buildMAFVCF/B73_Test.fa --maf-dir data/test/buildMAFVCF/mafs/ -o ${TestExtension.testVCFDir}")
        createMAFVCF.test("--db-path ${dbPath} --bed data/test/buildMAFVCF/B73_Test.bed --reference ${refFasta} --maf-dir data/test/buildMAFVCF/mafs/ -o ${TestExtension.testVCFDir}")

        //compare the contents of the output gVCF files to the expected output
        compareTwoGVCFFiles("data/test/buildMAFVCF/truthGVCFs/B97_truth.g.vcf", "${TestExtension.testVCFDir}/B97_ASM_Test.g.vcf")

        //Now we need to compare the hVCF's sequence with the sequence coming from the MAF files to make sure things match correctly as well as the boundaries
        val outputHVCF = "${TestExtension.testVCFDir}/B97_ASM_Test.h.vcf"
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
        //chr7	14	48 cooresponds to chr4 4245 .. 4281
        //chr7	450	456 cooresponds to chr4 5247 .. 5252
        //chr10	0	40 cooresponds to chr6 1098 .. 1142
        val truthCoords = setOf(Pair(Position("chr6", 98), Position("chr6", 142)),
            Pair(Position("chr4", 4245), Position("chr4", 4281)),
            Pair(Position("chr4", 5247), Position("chr4", 5252)),
            )

        val seqCoords = outputHeader.map {
            val tokens = it.toString().split("<")[1]
                .split(",")
                .filter { token -> token.startsWith("Asm") }
                .map { token -> token.split("=") }.associate { token -> Pair(token[0], token[1]) }

            val chr = tokens["Asm_Contig"]!!
            val start = tokens["Asm_Start"]!!.toInt()
            val end = tokens["Asm_End"]!!.toInt()
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
            Position("chr1",5), Position("chr1",10))

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
    fun testVariantFullyContained() {
        val createMafVCF = CreateMafVcf()
        val variant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10))

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
    fun testVariantPartiallyContainedStart() {
        val createMafVCF = CreateMafVcf()
        val variant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10))

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
    fun testVariantPartiallyContainedEnd() {
        val createMafVCF = CreateMafVcf()
        val variant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10))

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
    fun testVariantAfterRegion() {
        val createMafVCF = CreateMafVcf()
        val variant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10))

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
    fun testIsVariantResizable() {
        val createMafVCF = CreateMafVcf()
        val refBlockVariant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10))

        val multiAllelicSNPVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",7), Pair("AAA", "TTT"), Position("chr1",10), Position("chr1",12))

        val standardSNPVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",5), Pair("A", "T"), Position("chr1",10), Position("chr1",10))

        val insertionVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",5), Pair("A", "TTT"), Position("chr1",10), Position("chr1",12))
        val deletionVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",7), Pair("AAA", "T"), Position("chr1",10), Position("chr1",10))

        assertTrue(createMafVCF.isVariantResizable(refBlockVariant))
        assertTrue(createMafVCF.isVariantResizable(multiAllelicSNPVariant))
        assertTrue(createMafVCF.isVariantResizable(standardSNPVariant)) //This is technically true as single bp variants just return the bp when resized
        assertFalse(createMafVCF.isVariantResizable(insertionVariant))
        assertFalse(createMafVCF.isVariantResizable(deletionVariant))
    }

    /**
     * This unit test checks different cases of VariantContexts and the returned values depending on the requested ref
     * position to resize to.
     */
    @Test
    fun testResizeVariantContext() {
        val createMafVCF = CreateMafVcf()
        //Testing refBlocks
        val refBlockVariant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",15), Position("chr1",20))

        assertEquals(15, createMafVCF.resizeVariantContext(refBlockVariant, 5, "+"))
        assertEquals(20, createMafVCF.resizeVariantContext(refBlockVariant, 5, "-"))


        assertEquals(17, createMafVCF.resizeVariantContext(refBlockVariant, 7, "+"))
        assertEquals(18, createMafVCF.resizeVariantContext(refBlockVariant, 7, "-"))

        assertEquals(15, createMafVCF.resizeVariantContext(refBlockVariant, 2, "+"))
        assertEquals(20, createMafVCF.resizeVariantContext(refBlockVariant, 100, "+"))

        //Negative strand flips so we resize to the correct start/end on ASM
        assertEquals(15, createMafVCF.resizeVariantContext(refBlockVariant, 100, "-"))
        assertEquals(20, createMafVCF.resizeVariantContext(refBlockVariant, 2, "-"))

        assertEquals(-1, createMafVCF.resizeVariantContext(refBlockVariant, 7, "NOT_A_STRAND"))

        //Testing MultiAllelicPolymorphisms
        val multiAllelicSNPVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",10),
                                                Pair("AAAAA", "TTTTT"), Position("chr1",10), Position("chr1",15))
        //Check that we can resize the variant to the correct start/end on the ASM
        assertEquals(10, createMafVCF.resizeVariantContext(multiAllelicSNPVariant, 5, "+"))
        assertEquals(15, createMafVCF.resizeVariantContext(multiAllelicSNPVariant, 5, "-"))
        assertEquals(12, createMafVCF.resizeVariantContext(multiAllelicSNPVariant, 7, "+"))
        assertEquals(13, createMafVCF.resizeVariantContext(multiAllelicSNPVariant, 7, "-"))

        //Check for positions out of the record
        assertEquals(15, createMafVCF.resizeVariantContext(multiAllelicSNPVariant, 100, "+"))
        assertEquals(10, createMafVCF.resizeVariantContext(multiAllelicSNPVariant, 100, "-"))
        assertEquals(10, createMafVCF.resizeVariantContext(multiAllelicSNPVariant, 2, "+"))
        assertEquals(15, createMafVCF.resizeVariantContext(multiAllelicSNPVariant, 2, "-"))

        val standardSNPVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",5), Pair("A", "T"), Position("chr1",10), Position("chr1",10))
        //no matter what we request it should return 10 as it is only one bp of size
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 5, "+"))
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 5, "-"))
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 7, "+"))
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 7, "-"))
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 100, "+"))
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 100, "-"))
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 2, "+"))
        assertEquals(10, createMafVCF.resizeVariantContext(standardSNPVariant, 2, "-"))

        val insertionVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",5), Pair("A", "TTT"), Position("chr1",10), Position("chr1",12))
        //Everything should return -1 as its not resizable
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 5, "+"))
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 5, "-"))
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 7, "+"))
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 7, "-"))
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 100, "+"))
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 100, "-"))
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 2, "+"))
        assertEquals(-1, createMafVCF.resizeVariantContext(insertionVariant, 2, "-"))

        val deletionVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",7), Pair("AAA", "T"), Position("chr1",10), Position("chr1",10))
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