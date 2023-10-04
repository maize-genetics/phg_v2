package net.maizegenetics.phgv2.utils

import org.junit.jupiter.api.Assertions
import org.junit.jupiter.api.Assertions.assertNotEquals
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.assertEquals
//import net.maizegenetics.phgv2.utils.VariantLoadingUtils

class VariantLoadingUtilsTest {
    @Test
    fun testExportVariantContext() {
        // THis calls addSequenceDictionary, so that is tested as part of this
        // TBD
    }

    @Test
    fun testCreateGenericHeader() {
        // this also tests createGenericHeaderLineSet
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
    }

    @Test
    fun testCreateSNPVC() {

    }

    @Test
    fun testCreateRefRangeVC() {

    }

    @Test
    fun testVerifyIntervalRanges( ) {
        // NOTE: Ensure multiline strings are separated with an actual tab character!
        val testDataDir = "src/test/kotlin/net/maizegenetics/phgv2/testData/"
        val anchorFile = "${testDataDir}/testAnchorFile.txt"
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
    }




}