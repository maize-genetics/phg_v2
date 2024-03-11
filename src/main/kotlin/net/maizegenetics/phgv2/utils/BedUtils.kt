package net.maizegenetics.phgv2.utils

import biokotlin.genome.SeqRangeSort
import biokotlin.util.bufferedReader

/**
 * Simple function to load a BED file in.  This will be replaced by a lightweight Biokotlin ranges class eventually.
 *
 * This will sort in alphabetical order first then will check if there are numbers in the chromosome name and will
 * sort those numerically. This means that chr10 will come after chr2.
 */
fun loadRanges(bedFileName: String): List<Pair<Position, Position>> {

    return bufferedReader(bedFileName).readLines().map { line ->
        val lineSplit = line.split("\t")
        val chrom = lineSplit[0]
        val start = lineSplit[1].toInt() + 1
        val end = lineSplit[2].toInt()
        Pair(Position(chrom, start), Position(chrom, end))
    }
        .sortedWith(compareBy(SeqRangeSort.alphaThenNumberSort) { positionRange: Pair<Position, Position> -> positionRange.first.contig })

}