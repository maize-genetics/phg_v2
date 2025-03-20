package net.maizegenetics.phgv2.utils

import biokotlin.genome.AssemblyVariantInfo
import net.maizegenetics.phgv2.cli.CreateMafVcf
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class AssemblyVariantInfoUtilsTest {

    @Test
    fun testIsVariantInfoResizable() {
        val refBlockVariant = AssemblyVariantInfo("chr1",5,10,"T","A","T",true)

        val multiAllelicSNPVariant = AssemblyVariantInfo("chr1",5,7,"AAA","TTT","AAA",true)

        val standardSNPVariant = AssemblyVariantInfo("chr1",5,5,"T","A","T",true)

        val insertionVariant = AssemblyVariantInfo("chr1",5,5,"TTT","A","TTT",true)
        val deletionVariant = AssemblyVariantInfo("chr1",5,7,"T","AAA","T",true)

        assertTrue(AssemblyVariantInfoUtils.isVariantInfoResizable(refBlockVariant))
        assertTrue(AssemblyVariantInfoUtils.isVariantInfoResizable(multiAllelicSNPVariant))
        assertTrue(AssemblyVariantInfoUtils.isVariantInfoResizable(standardSNPVariant)) //This is technically true as single bp variants just return the bp when resized
        assertFalse(AssemblyVariantInfoUtils.isVariantInfoResizable(insertionVariant))
        assertFalse(AssemblyVariantInfoUtils.isVariantInfoResizable(deletionVariant))
    }



    @Test
    fun testResizeVariantInfo() {
        //Testing refBlocks
        val refBlockVariant = AssemblyVariantInfo("chr1",5,10,"T","A","T",false, intArrayOf(),"chr1", 15, 20)

        //Need to have reversed asm coords because that is how it looks with GVCFs coming from Biokotlin
        val refBlockVariantNegativeStrand = AssemblyVariantInfo("chr1",5,10,"T","A","T",false, intArrayOf(),"chr1", 20, 15)

        assertEquals(15, AssemblyVariantInfoUtils.resizeVariantInfo(refBlockVariant, 5, "+"))
        assertEquals(20, AssemblyVariantInfoUtils.resizeVariantInfo(refBlockVariantNegativeStrand, 5, "-"))

        assertEquals(17, AssemblyVariantInfoUtils.resizeVariantInfo(refBlockVariant, 7, "+"))
        assertEquals(18, AssemblyVariantInfoUtils.resizeVariantInfo(refBlockVariantNegativeStrand, 7, "-"))

        assertEquals(15, AssemblyVariantInfoUtils.resizeVariantInfo(refBlockVariant, 2, "+"))
        assertEquals(20, AssemblyVariantInfoUtils.resizeVariantInfo(refBlockVariant, 100, "+"))

        //Negative strand flips so we resize to the correct start/end on ASM
        assertEquals(15, AssemblyVariantInfoUtils.resizeVariantInfo(refBlockVariantNegativeStrand, 100, "-"))
        assertEquals(20, AssemblyVariantInfoUtils.resizeVariantInfo(refBlockVariantNegativeStrand, 2, "-"))

        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(refBlockVariant, 7, "NOT_A_STRAND"))


        //Testing MultiAllelicPolymorphisms
        val multiAllelicSNPVariant = AssemblyVariantInfo("chr1",5,10,"AAAAA","TTTTT","AAAAA",true, intArrayOf(),"chr1", 10, 15)
        val multiAllelicSNPVariantNegativeStrand = AssemblyVariantInfo("chr1",5,10,"AAAAA","TTTTT","AAAAA",true, intArrayOf(),"chr1", 15, 10)
        //Check that we can resize the variant to the correct start/end on the ASM
        assertEquals(10, AssemblyVariantInfoUtils.resizeVariantInfo(multiAllelicSNPVariant, 5, "+"))
        assertEquals(15, AssemblyVariantInfoUtils.resizeVariantInfo(multiAllelicSNPVariantNegativeStrand, 5, "-"))

        assertEquals(12, AssemblyVariantInfoUtils.resizeVariantInfo(multiAllelicSNPVariant, 7, "+"))
        assertEquals(13, AssemblyVariantInfoUtils.resizeVariantInfo(multiAllelicSNPVariantNegativeStrand, 7, "-"))

        //Check for positions out of the record
        assertEquals(15, AssemblyVariantInfoUtils.resizeVariantInfo(multiAllelicSNPVariant, 100, "+"))
        assertEquals(10, AssemblyVariantInfoUtils.resizeVariantInfo(multiAllelicSNPVariantNegativeStrand, 100, "-"))
        assertEquals(10, AssemblyVariantInfoUtils.resizeVariantInfo(multiAllelicSNPVariant, 2, "+"))
        assertEquals(15, AssemblyVariantInfoUtils.resizeVariantInfo(multiAllelicSNPVariantNegativeStrand, 2, "-"))

        val standardSNPVariant = AssemblyVariantInfo("chr1",5,5,"T","A","T",true, intArrayOf(),"chr1", 10, 10)
        //no matter what we request it should return 10 as it is only one bp of size
        assertEquals(10, AssemblyVariantInfoUtils.resizeVariantInfo(standardSNPVariant, 5, "+"))
        assertEquals(10, AssemblyVariantInfoUtils.resizeVariantInfo(standardSNPVariant, 5, "-"))
        assertEquals(10, AssemblyVariantInfoUtils.resizeVariantInfo(standardSNPVariant, 7, "+"))
        assertEquals(10, AssemblyVariantInfoUtils.resizeVariantInfo(standardSNPVariant, 7, "-"))
        assertEquals(10, AssemblyVariantInfoUtils.resizeVariantInfo(standardSNPVariant, 100, "+"))
        assertEquals(10, AssemblyVariantInfoUtils.resizeVariantInfo(standardSNPVariant, 100, "-"))
        assertEquals(10, AssemblyVariantInfoUtils.resizeVariantInfo(standardSNPVariant, 2, "+"))
        assertEquals(10, AssemblyVariantInfoUtils.resizeVariantInfo(standardSNPVariant, 2, "-"))

        val insertionVariant = AssemblyVariantInfo("chr1",5,5,"TTT","A","TTT",true, intArrayOf(),"chr1", 10, 12)
        //Everything should return -1 as its not resizable
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(insertionVariant, 5, "+"))
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(insertionVariant, 5, "-"))
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(insertionVariant, 7, "+"))
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(insertionVariant, 7, "-"))
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(insertionVariant, 100, "+"))
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(insertionVariant, 100, "-"))
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(insertionVariant, 2, "+"))
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(insertionVariant, 2, "-"))

        //val deletionVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",7), Pair("AAA", "T"), Position("chr1",10), Position("chr1",10),"+")
        val deletionVariant = AssemblyVariantInfo("chr1",5,7,"T","AAA","T",true, intArrayOf(),"chr1", 10, 10)
        //Everything should return -1 as its not resizeable
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(deletionVariant, 5, "+"))
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(deletionVariant, 5, "-"))
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(deletionVariant, 7, "+"))
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(deletionVariant, 7, "-"))
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(deletionVariant, 100, "+"))
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(deletionVariant, 100, "-"))
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(deletionVariant, 2, "+"))
        assertEquals(-1, AssemblyVariantInfoUtils.resizeVariantInfo(deletionVariant, 2, "-"))


    }

    @Test
    fun testBedRegionContainedInVariantInfo() {
        val variantInfos = AssemblyVariantInfo("chr1",5,10,"T","A","T",true)

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(AssemblyVariantInfoUtils.bedRegionContainedInVariantInfo(fullyContainedBed, variantInfos))
        assertTrue(AssemblyVariantInfoUtils.bedRegionContainedInVariantInfo(containedBed, variantInfos))
        assertFalse(AssemblyVariantInfoUtils.bedRegionContainedInVariantInfo(partiallyContainedStart, variantInfos))
        assertFalse(AssemblyVariantInfoUtils.bedRegionContainedInVariantInfo(partiallyContainedEnd, variantInfos))
        assertFalse(AssemblyVariantInfoUtils.bedRegionContainedInVariantInfo(notContained, variantInfos))
        assertFalse(AssemblyVariantInfoUtils.bedRegionContainedInVariantInfo(notContained2, variantInfos))
    }

    @Test
    fun testVariantInfoFullyContained() {
        val variantInfo = AssemblyVariantInfo("chr1",5,10,"T","A","T",true)

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertTrue(AssemblyVariantInfoUtils.variantInfoFullyContained(fullyContainedBed, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoFullyContained(containedBed, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoFullyContained(partiallyContainedStart, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoFullyContained(partiallyContainedEnd, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoFullyContained(notContained, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoFullyContained(notContained2, variantInfo))
    }

    @Test
    fun testVariantInfoPartiallyContainedStart() {
        val variantInfo = AssemblyVariantInfo("chr1",5,10,"T","A","T",true)

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(AssemblyVariantInfoUtils.variantInfoPartiallyContainedStart(fullyContainedBed, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoPartiallyContainedStart(containedBed, variantInfo))
        assertTrue(AssemblyVariantInfoUtils.variantInfoPartiallyContainedStart(partiallyContainedStart, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoPartiallyContainedStart(partiallyContainedEnd, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoPartiallyContainedStart(notContained, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoPartiallyContainedStart(notContained2, variantInfo))

    }

    @Test
    fun testVariantInfoPartiallyContainedEnd() {
        val variantInfo = AssemblyVariantInfo("chr1",5,10,"T","A","T",true)

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(AssemblyVariantInfoUtils.variantInfoPartiallyContainedEnd(fullyContainedBed, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoPartiallyContainedEnd(containedBed, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoPartiallyContainedEnd(partiallyContainedStart, variantInfo))
        assertTrue(AssemblyVariantInfoUtils.variantInfoPartiallyContainedEnd(partiallyContainedEnd, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoPartiallyContainedEnd(notContained, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoPartiallyContainedEnd(notContained2, variantInfo))

    }

    @Test
    fun testVariantInfoAfterRegion() {
        val variantInfo = AssemblyVariantInfo("chr1",5,10,"T","A","T",true)

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(AssemblyVariantInfoUtils.variantInfoAfterRegion(fullyContainedBed, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoAfterRegion(containedBed, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoAfterRegion(partiallyContainedStart, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoAfterRegion(partiallyContainedEnd, variantInfo))
        assertFalse(AssemblyVariantInfoUtils.variantInfoAfterRegion(notContained2, variantInfo))
        assertTrue(AssemblyVariantInfoUtils.variantInfoAfterRegion(notContained, variantInfo))

    }
}