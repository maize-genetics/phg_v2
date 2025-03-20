package net.maizegenetics.phgv2.utils

import biokotlin.seq.NucSeq
import net.maizegenetics.phgv2.cli.CreateMafVcf
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class VariantContextUtilsTest {

    @Test
    fun testIsVariantResizable() {
        val refBlockVariant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10),"+")

        val multiAllelicSNPVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",7), Pair("AAA", "TTT"), Position("chr1",10), Position("chr1",12),"+")

        val standardSNPVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",5), Pair("A", "T"), Position("chr1",10), Position("chr1",10),"+")

        val insertionVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",5), Pair("A", "TTT"), Position("chr1",10), Position("chr1",12),"+")
        val deletionVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",7), Pair("AAA", "T"), Position("chr1",10), Position("chr1",10),"+")

        assertTrue(VariantContextUtils.isVariantResizable(refBlockVariant))
        assertTrue(VariantContextUtils.isVariantResizable(multiAllelicSNPVariant))
        assertTrue(VariantContextUtils.isVariantResizable(standardSNPVariant)) //This is technically true as single bp variants just return the bp when resized
        assertFalse(VariantContextUtils.isVariantResizable(insertionVariant))
        assertFalse(VariantContextUtils.isVariantResizable(deletionVariant))
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
            Position("chr1",15), Position("chr1",20),"+")

        //Need to have reveresed asm coords because that is how it looks with GVCFs coming from Biokotlin
        val refBlockVariantNegativeStrand = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",20), Position("chr1",15),"+")

        assertEquals(15, VariantContextUtils.resizeVariantContext(refBlockVariant, 5, "+"))
        assertEquals(20, VariantContextUtils.resizeVariantContext(refBlockVariantNegativeStrand, 5, "-"))


        assertEquals(17, VariantContextUtils.resizeVariantContext(refBlockVariant, 7, "+"))
        assertEquals(18, VariantContextUtils.resizeVariantContext(refBlockVariantNegativeStrand, 7, "-"))

        assertEquals(15, VariantContextUtils.resizeVariantContext(refBlockVariant, 2, "+"))
        assertEquals(20, VariantContextUtils.resizeVariantContext(refBlockVariant, 100, "+"))

        //Negative strand flips so we resize to the correct start/end on ASM
        assertEquals(15, VariantContextUtils.resizeVariantContext(refBlockVariantNegativeStrand, 100, "-"))
        assertEquals(20, VariantContextUtils.resizeVariantContext(refBlockVariantNegativeStrand, 2, "-"))

        assertEquals(-1, VariantContextUtils.resizeVariantContext(refBlockVariant, 7, "NOT_A_STRAND"))

        //Testing MultiAllelicPolymorphisms
        val multiAllelicSNPVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",10),
            Pair("AAAAA", "TTTTT"), Position("chr1",10), Position("chr1",15),"+")
        val multiAllelicSNPVariantNegativeStrand = createSNPVC("B97", Position("chr1",5),Position("chr1",10),
            Pair("AAAAA", "TTTTT"), Position("chr1",15), Position("chr1",10),"+")
        //Check that we can resize the variant to the correct start/end on the ASM
        assertEquals(10, VariantContextUtils.resizeVariantContext(multiAllelicSNPVariant, 5, "+"))
        assertEquals(15, VariantContextUtils.resizeVariantContext(multiAllelicSNPVariantNegativeStrand, 5, "-"))
        assertEquals(12, VariantContextUtils.resizeVariantContext(multiAllelicSNPVariant, 7, "+"))
        assertEquals(13, VariantContextUtils.resizeVariantContext(multiAllelicSNPVariantNegativeStrand, 7, "-"))

        //Check for positions out of the record
        assertEquals(15, VariantContextUtils.resizeVariantContext(multiAllelicSNPVariant, 100, "+"))
        assertEquals(10, VariantContextUtils.resizeVariantContext(multiAllelicSNPVariantNegativeStrand, 100, "-"))
        assertEquals(10, VariantContextUtils.resizeVariantContext(multiAllelicSNPVariant, 2, "+"))
        assertEquals(15, VariantContextUtils.resizeVariantContext(multiAllelicSNPVariantNegativeStrand, 2, "-"))

        val standardSNPVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",5), Pair("A", "T"), Position("chr1",10), Position("chr1",10),"+")
        //no matter what we request it should return 10 as it is only one bp of size
        assertEquals(10, VariantContextUtils.resizeVariantContext(standardSNPVariant, 5, "+"))
        assertEquals(10, VariantContextUtils.resizeVariantContext(standardSNPVariant, 5, "-"))
        assertEquals(10, VariantContextUtils.resizeVariantContext(standardSNPVariant, 7, "+"))
        assertEquals(10, VariantContextUtils.resizeVariantContext(standardSNPVariant, 7, "-"))
        assertEquals(10, VariantContextUtils.resizeVariantContext(standardSNPVariant, 100, "+"))
        assertEquals(10, VariantContextUtils.resizeVariantContext(standardSNPVariant, 100, "-"))
        assertEquals(10, VariantContextUtils.resizeVariantContext(standardSNPVariant, 2, "+"))
        assertEquals(10, VariantContextUtils.resizeVariantContext(standardSNPVariant, 2, "-"))

        val insertionVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",5), Pair("A", "TTT"), Position("chr1",10), Position("chr1",12),"+")
        //Everything should return -1 as its not resizable
        assertEquals(-1, VariantContextUtils.resizeVariantContext(insertionVariant, 5, "+"))
        assertEquals(-1, VariantContextUtils.resizeVariantContext(insertionVariant, 5, "-"))
        assertEquals(-1, VariantContextUtils.resizeVariantContext(insertionVariant, 7, "+"))
        assertEquals(-1, VariantContextUtils.resizeVariantContext(insertionVariant, 7, "-"))
        assertEquals(-1, VariantContextUtils.resizeVariantContext(insertionVariant, 100, "+"))
        assertEquals(-1, VariantContextUtils.resizeVariantContext(insertionVariant, 100, "-"))
        assertEquals(-1, VariantContextUtils.resizeVariantContext(insertionVariant, 2, "+"))
        assertEquals(-1, VariantContextUtils.resizeVariantContext(insertionVariant, 2, "-"))

        val deletionVariant = createSNPVC("B97", Position("chr1",5),Position("chr1",7), Pair("AAA", "T"), Position("chr1",10), Position("chr1",10),"+")
        //Everything should return -1 as its not resizeable
        assertEquals(-1, VariantContextUtils.resizeVariantContext(deletionVariant, 5, "+"))
        assertEquals(-1, VariantContextUtils.resizeVariantContext(deletionVariant, 5, "-"))
        assertEquals(-1, VariantContextUtils.resizeVariantContext(deletionVariant, 7, "+"))
        assertEquals(-1, VariantContextUtils.resizeVariantContext(deletionVariant, 7, "-"))
        assertEquals(-1, VariantContextUtils.resizeVariantContext(deletionVariant, 100, "+"))
        assertEquals(-1, VariantContextUtils.resizeVariantContext(deletionVariant, 100, "-"))
        assertEquals(-1, VariantContextUtils.resizeVariantContext(deletionVariant, 2, "+"))
        assertEquals(-1, VariantContextUtils.resizeVariantContext(deletionVariant, 2, "-"))
    }

    @Test
    fun testBedRegionContainedInVariant() {
        val variant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10),"+")

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(VariantContextUtils.bedRegionContainedInVariant(fullyContainedBed, variant))
        assertTrue(VariantContextUtils.bedRegionContainedInVariant(containedBed, variant))
        assertFalse(VariantContextUtils.bedRegionContainedInVariant(partiallyContainedStart, variant))
        assertFalse(VariantContextUtils.bedRegionContainedInVariant(partiallyContainedEnd, variant))
        assertFalse(VariantContextUtils.bedRegionContainedInVariant(notContained, variant))
        assertFalse(VariantContextUtils.bedRegionContainedInVariant(notContained2, variant))
    }

    @Test
    fun testVariantFullyContained() {
        val variant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10),"+")

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertTrue(VariantContextUtils.variantFullyContained(fullyContainedBed, variant))
        assertFalse(VariantContextUtils.variantFullyContained(containedBed, variant))
        assertFalse(VariantContextUtils.variantFullyContained(partiallyContainedStart, variant))
        assertFalse(VariantContextUtils.variantFullyContained(partiallyContainedEnd, variant))
        assertFalse(VariantContextUtils.variantFullyContained(notContained, variant))
        assertFalse(VariantContextUtils.variantFullyContained(notContained2, variant))

    }

    @Test
    fun testVariantPartiallyContainedStart() {
        val variant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10),"+")

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(VariantContextUtils.variantPartiallyContainedStart(fullyContainedBed, variant))
        assertFalse(VariantContextUtils.variantPartiallyContainedStart(containedBed, variant))
        assertTrue(VariantContextUtils.variantPartiallyContainedStart(partiallyContainedStart, variant))
        assertFalse(VariantContextUtils.variantPartiallyContainedStart(partiallyContainedEnd, variant))
        assertFalse(VariantContextUtils.variantPartiallyContainedStart(notContained, variant))
        assertFalse(VariantContextUtils.variantPartiallyContainedStart(notContained2, variant))
    }


    @Test
    fun testVariantPartiallyContainedEnd() {
        val variant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10),"+")

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(VariantContextUtils.variantPartiallyContainedEnd(fullyContainedBed, variant))
        assertFalse(VariantContextUtils.variantPartiallyContainedEnd(containedBed, variant))
        assertFalse(VariantContextUtils.variantPartiallyContainedEnd(partiallyContainedStart, variant))
        assertTrue(VariantContextUtils.variantPartiallyContainedEnd(partiallyContainedEnd, variant))
        assertFalse(VariantContextUtils.variantPartiallyContainedEnd(notContained, variant))
        assertFalse(VariantContextUtils.variantPartiallyContainedEnd(notContained2, variant))
    }

    @Test
    fun testVariantAfterRegion() {
        val variant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",5), Position("chr1",10),"+")

        val fullyContainedBed = Pair(Position("chr1", 3), Position("chr1", 15))
        val containedBed = Pair(Position("chr1",6),Position("chr1",9))
        val partiallyContainedStart = Pair(Position("chr1",4),Position("chr1",9))
        val partiallyContainedEnd = Pair(Position("chr1",6),Position("chr1",11))
        val notContained = Pair(Position("chr1",1),Position("chr1",4))
        val notContained2 = Pair(Position("chr1",11),Position("chr1",15))

        assertFalse(VariantContextUtils.variantAfterRegion(fullyContainedBed, variant))
        assertFalse(VariantContextUtils.variantAfterRegion(containedBed, variant))
        assertFalse(VariantContextUtils.variantAfterRegion(partiallyContainedStart, variant))
        assertFalse(VariantContextUtils.variantAfterRegion(partiallyContainedEnd, variant))
        assertFalse(VariantContextUtils.variantAfterRegion(notContained2, variant))
        assertTrue(VariantContextUtils.variantAfterRegion(notContained, variant))
    }


    @Test
    fun testGVCFResizingIssue() {
        //Issue found with one of Matt's cassava lines
        //I think the issue is that we have a refBlock variant which is before the end of the bed region which then attempts to get resized.
        //The resulting coordiantes fall within the negative boundaries which is not correct.

        // |-------BED------------|
        //           |--GCVF--|

        //create a bedRange
        //Chromosome04	35928654	35931063
        val bedRange = Pair(Position("Chromosome04",35928654), Position("Chromosome04",35931063))


        //Test forward strand - working initially

        val refBlockVariantForward = createRefRangeVC(mapOf("Chromosome04" to NucSeq("A".repeat(35929999))),"cassava_test",
            Position("Chromosome04",35920437), Position("Chromosome04",35929999),
            Position("chr8",4514), Position("chr8",14076),"+")

        val resizedVariantStartForward = VariantContextUtils.resizeVariantContext(refBlockVariantForward, bedRange.first.position, "+")
        assertEquals(12731, resizedVariantStartForward)

        val resizedVariantEndForward = VariantContextUtils.resizeVariantContext(refBlockVariantForward, bedRange.second.position, "+")
        assertEquals(14076, resizedVariantEndForward)

        //Make a variant context based on this:
//        Chromosome04	35920437	.	C	<NON_REF>	.	.	ASM_Chr=chr8;ASM_End=4514;ASM_Start=14076;ASM_Strand=-;END=35929999	GT:AD:DP:PL	0:30,0:30:0,90,90

        //Need to have reversed asm coords because that is how it looks with GVCFs coming from Biokotlin
        val refBlockVariant = createRefRangeVC(mapOf("Chromosome04" to NucSeq("A".repeat(35929999))),"cassava_test",
            Position("Chromosome04",35920437), Position("Chromosome04",35929999),
            Position("chr8",14076), Position("chr8",4514),"+")


        val resizedVariantStart = VariantContextUtils.resizeVariantContext(refBlockVariant, bedRange.first.position, "-")
        assertEquals(5859, resizedVariantStart)

        val resizedVariantEnd = VariantContextUtils.resizeVariantContext(refBlockVariant, bedRange.second.position, "-")
        assertEquals(4514, resizedVariantEnd)

    }

}