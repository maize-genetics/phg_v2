package net.maizegenetics.phgv2.utils

import biokotlin.seq.NucSeq
import com.google.common.io.Files
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFAltHeaderLine
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeaderVersion
import io.kotest.core.spec.style.AnnotationSpec
import org.apache.logging.log4j.core.util.FileUtils
import org.junit.jupiter.api.*
import org.junit.jupiter.api.Assertions.assertNotEquals
import java.io.File
import java.lang.IllegalArgumentException
import java.lang.IllegalStateException
import kotlin.test.assertEquals
import kotlin.test.Ignore
//import net.maizegenetics.phgv2.utils.VariantLoadingUtils

class VariantLoadingUtilsTest {
    companion object {
        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"

        @JvmStatic
        @BeforeAll
        fun setup() {
            File(tempDir).mkdirs()
        }

        // Comment out the tearDown() you need to look at the logs files or other files
        // created by the junit tests
        @JvmStatic
        @AfterAll
        fun teardown() {
            File(tempDir).deleteRecursively()
        }
    }

    @Ignore
    @Test
    fun testBgzipAndIndex() {
        // Copy the sample gvcf file to the test directory
        val origGvcfFile = "data/test/smallseq/sample.gvcf"
        val testGvcfFile = tempDir + "sample.gvcf"
        Files.copy(File(origGvcfFile), File(testGvcfFile))

        // call bgzipAndIndexGVCFfile to zip and index the file
        bgzipAndIndexGVCFfile(testGvcfFile)

        // check that the compressed file exists
        val gvcfFileZipped = File(testGvcfFile + ".gz")
        assertEquals(true, gvcfFileZipped.exists())

        // check that the index file exists
        val gvcfFileZippedIndex = File(testGvcfFile + ".gz.csi")
        assertEquals(true, gvcfFileZippedIndex.exists())
    }

    @Test
    fun testExportVariantContext() {
        //Create some simple VariantContext records into a list
        val refMap = mapOf("1" to NucSeq("AAAAAAAAAAA"), "2" to NucSeq("TTTTTTTTTT"))
        val variants = listOf(
            createRefRangeVC(refMap, "Line1", Position("1",1), Position("1",1), Position("1",1), Position("1",1)),
            createSNPVC("Line1", Position("1",2), Position("1",2), Pair("A","T"), Position("1",2), Position("1",2)),
            createRefRangeVC(refMap,"Line1", Position("1",3), Position("1",6), Position("1",3), Position("1",6)),
            createSNPVC("Line1", Position("1",7), Position("1",7), Pair("A","G"), Position("1",2), Position("1",2)),
            createRefRangeVC(refMap,"Line1", Position("1",8), Position("1",10), Position("1",8), Position("1",10))
        )

        //Export the variants to a file
        val testFile = "src/test/kotlin/net/maizegenetics/phgv2/testData/testExportVariantContext.vcf"
        exportVariantContext("Line1", variants, testFile, refMap, setOf())

        //Load the file back in and check that the variants are as expected
        val loadedVariants = VCFFileReader(File(testFile), false).iterator().asSequence().toList()

        assertEquals(5, loadedVariants.size, "Unexpected number of variants")
        for(i in 0 until variants.size) {
            compareVariants(variants[i], loadedVariants[i])
        }

    }

    fun compareVariants(expectedVariant: VariantContext, loadedVariant: VariantContext) {
        //compare the two variant contexts
        assertEquals(expectedVariant.contig, loadedVariant.contig, "Contig not as expected")
        assertEquals(expectedVariant.start, loadedVariant.start, "Start position not as expected")
        assertEquals(expectedVariant.end, loadedVariant.end, "End position not as expected")
        assertEquals(expectedVariant.genotypes[0].sampleName, loadedVariant.genotypes[0].sampleName, "Sample name not as expected")
        assertEquals(expectedVariant.alleles[0].baseString, loadedVariant.alleles[0].baseString, "Reference allele not as expected")
        assertEquals(expectedVariant.alternateAlleles[0].baseString, loadedVariant.alternateAlleles[0].baseString, "Alternate allele not as expected")
        assertEquals(expectedVariant.getAttributeAsString("ASM_Chr",""), loadedVariant.getAttributeAsString("ASM_Chr",""), "ASM_Chr not as expected")
        assertEquals(expectedVariant.getAttributeAsInt("ASM_Start",0), loadedVariant.getAttributeAsInt("ASM_Start",0), "ASM_Start not as expected")
        assertEquals(expectedVariant.getAttributeAsInt("ASM_End",0), loadedVariant.getAttributeAsInt("ASM_End",0), "ASM_End not as expected")

    }

    @Test
    fun testCreateGenericHeader() {
        // this also tests createGenericHeaderLineSet
        val taxa = listOf("LineA", "LineB", "LineC")
        val altHeaders = (1 .. 10).map{VCFAltHeaderLine("<ID=${it}, Description=\"Description for ${it}\">", VCFHeaderVersion.VCF4_2)}.toSet()

        //create the headers with the taxa and the alternative headers:
        val testHeader = createGenericHeader(taxa, altHeaders)

        // Check that the header lines are as expected
        val headerAltIdLines = testHeader.idHeaderLines.map { it.toString() }.filter { it.startsWith("ALT") }
        assertEquals(10, headerAltIdLines.size, "Unexpected number of ALT header lines")
        for (i in 1 .. 10) {
            val expectedHeader = "ALT=<ID=${i},Description=\"Description for ${i}\">"
            assertEquals(true, headerAltIdLines.contains(expectedHeader), "Expected header line not found:\n${expectedHeader}")
        }

        //Check the samples as well
        val headerSampleLines = testHeader.sampleNamesInOrder
        assertEquals(3, headerSampleLines.size, "Unexpected number of samples")
        for (i in 0 .. 2) {
            val expectedHeader = taxa[i]
            assertEquals(true, headerSampleLines.contains(expectedHeader), "Expected header line not found:\n${expectedHeader}")
        }
    }

    @Test
    fun testCreateGenericHeaderLineSet() {
        val testHeaderLineSet = createGenericHeaderLineSet()
        assertEquals(13, testHeaderLineSet.size, "Unexpected number of header lines")

        // Setup the expected header strings
        val expectedHeaderLines = listOf(
            "FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "FORMAT=<ID=AD,Number=3,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">",
            "FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">",
            "FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
            "FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">",
            "INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
            "INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">",
            "INFO=<ID=AF,Number=3,Type=Integer,Description=\"Allele Frequency\">",
            "INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">",
            "INFO=<ID=ASM_Chr,Number=1,Type=String,Description=\"Assembly chromosome\">",
            "INFO=<ID=ASM_Start,Number=1,Type=Integer,Description=\"Assembly start position\">",
            "INFO=<ID=ASM_End,Number=1,Type=Integer,Description=\"Assembly end position\">",
            "INFO=<ID=ASM_Strand,Number=1,Type=String,Description=\"Assembly strand\">"
        )
        val vcfHeaderStrings = testHeaderLineSet.map { it.toString() }
        // Check that the header lines are as expected
        for(header in expectedHeaderLines) {
            assertEquals(true, vcfHeaderStrings.contains(header), "Expected header line not found:\n${header}")
        }
    }

    @Test
    fun testAddSequenceDictionary() {
        //fun addSequenceDictionary(vcfheader : VCFHeader, refGenomeSequence: Map<String,NucSeq>) {
        //Create a generic vcfHeader
        val testHeader = createGenericHeader(listOf("LineA", "LineB", "LineC"), setOf())
        //create a dummy refGenomeSequence
        val refGenomeSequence = mapOf("1" to NucSeq("AAAAAAAAAAA"), "2" to NucSeq("TTTTTTTTTT"))
        addSequenceDictionary(testHeader, refGenomeSequence)

        //Compare TestHeader's seqDict to expected
        testHeader.sequenceDictionary.sequences.forEach { seq ->
            assertEquals(true, refGenomeSequence.containsKey(seq.sequenceName), "Sequence name not found in refGenomeSequence")
            assertEquals(refGenomeSequence[seq.sequenceName]!!.size(), seq.sequenceLength, "Sequence length not as expected")
        }
    }

    @Test
    fun testGetChecksumForString() {
        // this mostly verifies that different values are returned for different strings
        val testString1 = "AGCGGTTAAGGGGTTACACACACACATGTGTGTTTTTGGGGGGGGGGGGGGAAAAAAAAACACACACAC"
        val testString2 = "CCCCCCCTTTTTTTAAAAAAAGTGATCGATCGTACGTACGTACTACTACGTACGTACGTACTACGTACA"
        // Get checksums for each string
        val chrom = "1"
        val position = 1
        val testString1Checksum = getChecksumForString(testString1 )
        val testString2Checksum = getChecksumForString(testString2)
        // Verify that the checksums are different
        assertNotEquals(testString1Checksum, testString2Checksum)

        assertThrows<IllegalStateException> { getChecksumForString(testString1,"fake_protocol") }
    }

    @Test
    fun testCreateSNPVC() {
        //createSNPVC(assemblyTaxon: String, startPosition: Position, endPosition: Position, calls: Pair<String, String>,
        //                asmStart: Position, asmEnd: Position)
        val vc1 = createSNPVC("Line1", Position("1",2), Position("1",2), Pair("A","T"), Position("1",2), Position("1",2))

        assertEquals("1", vc1.contig)
        assertEquals(2, vc1.start)
        assertEquals(2, vc1.end)
        assertEquals("Line1", vc1.genotypes[0].sampleName)
        assertEquals("A", vc1.alleles[0].baseString)
        assertEquals("T", vc1.alternateAlleles[0].baseString)
        assertEquals("1", vc1.getAttributeAsString("ASM_Chr",""))
        assertEquals(2, vc1.getAttributeAsInt("ASM_Start",0))
        assertEquals(2, vc1.getAttributeAsInt("ASM_End",0))


        assertThrows<IllegalStateException> {
            //Check to make sure the reference Coords are not backwards
            createSNPVC("Line1", Position("1",4), Position("1",2), Pair("A","T"), Position("1",3), Position("1",2))
        }
        assertThrows<IllegalStateException> {
            //Check to make sure the assembly coords exist on same ASM chrom
            createSNPVC("Line1", Position("1",2), Position("1",2), Pair("A","T"), Position("1",2), Position("2",3))
        }

    }

    @Test
    fun testCreateRefRangeVC() {
        val refSeq = mapOf("1" to NucSeq("AAAAAAAAAAA"), "2" to NucSeq("TTTTTTTTTTT"))
        val vc1 = createRefRangeVC(refSeq, "Line1", Position("1",2), Position("1",5), Position("1",3), Position("1",6))

        assertEquals("1", vc1.contig)
        assertEquals(2, vc1.start)
        assertEquals(5, vc1.end)
        assertEquals("Line1", vc1.genotypes[0].sampleName)
        assertEquals("1", vc1.getAttributeAsString("ASM_Chr",""))
        assertEquals(3, vc1.getAttributeAsInt("ASM_Start",0))
        assertEquals(6, vc1.getAttributeAsInt("ASM_End",0))

        assertThrows<IllegalStateException> {
            //Check that it throws an error when the refCoordinates are backwards
            createRefRangeVC(refSeq, "Line1", Position("1",10), Position("1",5), Position("1",10), Position("1",5))
        }
        assertThrows<IllegalStateException> {
            //Check that it throws an error when the assembly coords span 2 chroms
            createRefRangeVC(refSeq, "Line1", Position("1",2), Position("1",5), Position("1",10), Position("2",5))
        }

    }

    @Test
    fun testVerifyIntervalRanges( ) {
        // NOTE: Ensure multiline strings are separated with an actual tab character!
        val testDataDir = "src/test/kotlin/net/maizegenetics/phgv2/testData/"
        val anchorFile = "${testDataDir}/testAnchorFile.txt"
        val incorrectAnchorFile = "${testDataDir}/testAnchorFileIncorrect.txt"

        File(anchorFile).bufferedWriter().use {
            // Lines 4 and 5 overlap line 3
            // Line 8 overlaps line 7 (First chr2 anchor)
            // Line 13 overlaps line 12 (but different chrom, so should not count)
            val anchorContents = """
                chr1	10575	13198
                chr1	20460	29234
                chr1	141765	145627
                chr1	143661	146302
                chr1	144363	148863
                chr1	175219	177603
                chr2	81176	84375
                chr2	82776	87024
                chr2	108577	113286
                chr3	116671	122941
                chr3	140393	145805
                chr3	159053	164568
                chr5	158053	164568
            """.trimIndent()
            it.write(anchorContents)
        }

        val overlaps1 = verifyIntervalRanges(anchorFile)
        Assertions.assertEquals(3, overlaps1.size)

        val containsEntry1 = overlaps1.contains("chr1\t143661\t146302");
        Assertions.assertEquals(containsEntry1, true);

        val containsEntry2 = overlaps1.contains("chr5\t158053\t164568");
        Assertions.assertEquals(false, containsEntry2)

        // Do again, but have the intervals out of order in the anchor file.
        File(anchorFile).bufferedWriter().use {
            // Line 3 overlaps line 1
            // Line 5 overlaps lines 3 and 4
            // Line 8 overlaps line 7
            // Line 13 overlaps line 12 (but different chrom, so should not count)
            val anchorContents = """
                chr1	141765	145627
                chr1	20460	29234
                chr1	143661	146302
                chr1	10575	13198
                chr1	144363	148863
                chr1	175219	177603
                chr2	82776	87024
                chr2	81176	84375
                chr2	108577	113286
                chr3	116671	122941
                chr3	140393	145805
                chr3	159053	164568
                chr5	158053	164568
            """.trimIndent()
            it.write(anchorContents)
        }

        val overlaps2 = verifyIntervalRanges(anchorFile)
        Assertions.assertEquals(3, overlaps2.size)

        val containsEntry3 = overlaps2.contains("chr1\t143661\t146302");
        Assertions.assertEquals(containsEntry3, true);

        val containsEntry4 = overlaps2.contains("chr5\t158053\t164568");
        Assertions.assertEquals(containsEntry4, false);

        // Write anchor file with NO overlaps
        // Do again, but have the intervals out of order in the anchor file.
        File(anchorFile).bufferedWriter().use {
            val anchorContents = """
                chr1	10575	13198
                chr1	20460	29234
                chr1	141765	145627
                chr1	146661	147302
                chr1	148363	150863
                chr1	175219	177603
                chr2	81176	81375
                chr2	82776	87024
                chr2	108577	113286
                chr3	116671	122941
                chr3	140393	145805
                chr3	159053	164568
                chr5	158053	164568
            """.trimIndent()
            it.write(anchorContents)
        }

        val overlaps3 = verifyIntervalRanges(anchorFile);
        Assertions.assertEquals(0, overlaps3.size)

        assertThrows<IllegalArgumentException> { verifyIntervalRanges(incorrectAnchorFile) }
    }




}