package net.maizegenetics.phgv2.utils

import org.junit.jupiter.api.Assertions.assertNotEquals
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
//import net.maizegenetics.phgv2.utils.VariantLoadingUtils

class VariantLoadingUtilsTest {
    @Test
    fun testGetCheckSumForString() {
        // this mostly verifies that different values are returned for different strings
        val testString1 = "AGCGGTTAAGGGGTTACACACACACATGTGTGTTTTTGGGGGGGGGGGGGGAAAAAAAAACACACACAC"
        val testString2 = "CCCCCCCTTTTTTTAAAAAAAGTGATCGATCGTACGTACGTACTACTACGTACGTACGTACTACGTACA"
        // Get checksums for each string
        val chrom = "1"
        val position = 1
        val posClass = Position(chrom,position)
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
    fun testVerifyIntervalRanges() {

    }

    @Test
    fun testCreateGenericHeader() {

    }

}