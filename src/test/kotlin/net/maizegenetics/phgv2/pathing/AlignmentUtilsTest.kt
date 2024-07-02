package net.maizegenetics.phgv2.pathing

import org.junit.jupiter.api.Test
import org.junit.jupiter.api.Assertions.assertEquals

class AlignmentUtilsTest {

    @Test
    fun testFilterHapIdsByKmerCount() {
        val hapIds = listOf<String>()
        val result = AlignmentUtils.filterHapIdsByKmerCount(hapIds, 1.0)
        assertEquals(setOf<String>(), result, "Expected an empty set for an empty input list")
    }
    @Test
    fun testFilterHapIdsByKmerCountUniqueIds() {
        val hapIds = listOf("hap1", "hap2", "hap3")
        val result = AlignmentUtils.filterHapIdsByKmerCount(hapIds, 0.5)
        assertEquals(setOf("hap1", "hap2", "hap3"), result, "Expected all hapIds to be returned when all are unique")
    }
    @Test
    fun testFilterHapIdsByKmerCountBelowThreshold() {
        val hapIds = listOf("hap1", "hap1","hap1", "hap2")
        val result = AlignmentUtils.filterHapIdsByKmerCount(hapIds, 0.5)
        assertEquals(setOf("hap1"), result, "Expected only hapIds above the threshold to be returned")
    }

    @Test
    fun testFilterHapIdsByKmerCountMinPropCount() {
        val hapIds = listOf("hap1", "hap2", "hap2", "hap2", "hap3")
        val result = AlignmentUtils.filterHapIdsByKmerCount(hapIds, 0.66)
        assertEquals(setOf("hap2"), result, "Expected only hapIds that meet the min proportion of max count to be returned")
    }

    @Test
    fun testFilterHapIdsByKmerCountAllEqual() {
        val hapIds = listOf("hap1", "hap1", "hap1")
        val result = AlignmentUtils.filterHapIdsByKmerCount(hapIds, 1.0)
        assertEquals(setOf("hap1"), result, "Expected the single unique hapId to be returned when all hapIds are equal")
    }
}


