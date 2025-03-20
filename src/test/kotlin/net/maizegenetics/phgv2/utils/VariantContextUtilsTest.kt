package net.maizegenetics.phgv2.utils

import biokotlin.seq.NucSeq
import org.junit.jupiter.api.Test
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class VariantContextUtilsTest {

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


}