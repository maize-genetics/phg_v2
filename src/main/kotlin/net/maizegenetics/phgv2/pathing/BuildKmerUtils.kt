package net.maizegenetics.phgv2.pathing

/**
 * This file contains utility functions for building or processing  kmer index code
 */

// This function is used to set a bit in a ULong for a sample in a sorted list of samples.
// The bit is set at the position corresponding to the index of the sample in the list
// and is used to store hits for the keepMap of the processGraphKmers2 function
fun setBitForSample(gamete: ULong, samples: List<String>, sample: String): ULong {
    val index = samples.indexOf(sample)
    if (index == -1) {
        throw IllegalArgumentException("Sample not found in the list.")
    }
    // Set the bit at the position corresponding to the index
    return gamete or (1uL shl index)
}