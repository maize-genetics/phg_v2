package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.api.SampleGamete

/**
 * This file contains utility functions for building or processing  kmer index code
 */

/**
 * This function is used to set a bit in a ULong for a sample in a sorted list of samples.
 * The bit is set at the position corresponding to the index of the sample in the list
 * and is used to store hits for the keepMap of the processGraphKmers function of BuildRanEfficientKmerIndex.kt
 */
fun setBitForSample(gamete: ULong, samples: List<String>, sample: String): ULong {
    val index = samples.indexOf(sample)
    if (index == -1) {
        throw IllegalArgumentException("Sample not found in the list.")
    }
    // Set the bit at the position corresponding to the index
    return gamete or (1uL shl index)
}

/**
 *  This function returns the indices of the set bits in a ULong number
 *  The iteration stops once the number becomes zero, it avoids unnecessary checks when bits are not set.
 *  It is used to get the samples that have hits in the keepMap of the processGraphKmers2 function
 */

fun getSetBitIndices(number: Long): List<Int> {
    val indices = mutableListOf<Int>()
    var value = number
    var index = 0

    while (value > 0L) {
        if ((value and 1L) != 0L) {
            indices.add(index)
        }
        value = value shr 1
        index++
    }

    return indices
}

/**
 * Given two ULong numbers representing a bit set of samples, and a map of samples to indexed positions,
 * this function returns a list of SampleGamete objects represented by the set bits in the ULong numbers.
 */

fun findSampleGametes(
    gameteLong0: Long,
    gameteLong1: Long,
    samplesMap: Map<String, Int>
): Set<SampleGamete>{
    val sampleGametes = mutableSetOf<SampleGamete>()

    // Combine gameteLong0 and gameteLong1 into a single map lookup
    val allSetBits = getSetBitIndices(gameteLong0) + getSetBitIndices(gameteLong1).map { it + 64 }

    // Map indices to SampleGametes using samplesMap
    for (bitIndex in allSetBits) {
        val sampleName = samplesMap.entries.find { it.value == bitIndex }?.key
        sampleName?.let { sampleGametes.add(parseSampleGamete(it)) }
    }

    return sampleGametes
}

/**
 * Given a string representation of a SampleGamete, this function returns a SampleGamete object.
 */
fun parseSampleGamete(sampleString: String): SampleGamete {
    val parts = sampleString.split(":")
    val name = parts[0]
    val gameteId = parts.getOrNull(1)?.toIntOrNull() ?: 0
    return SampleGamete(name, gameteId)
}