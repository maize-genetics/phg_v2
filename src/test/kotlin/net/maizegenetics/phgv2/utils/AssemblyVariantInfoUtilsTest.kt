package net.maizegenetics.phgv2.utils

import biokotlin.genome.AssemblyVariantInfo
import org.junit.jupiter.api.Test
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class AssemblyVariantInfoUtilsTest {
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