package net.maizegenetics.phgv2.utils

import biokotlin.seq.NucSeq
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import com.google.common.collect.TreeRangeMap
import org.apache.logging.log4j.LogManager
import java.io.BufferedInputStream

private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.utils.GeneralUtilities")

/**
 * This function reads the output from a ProcessBuilder cammond.
 * The input is a BufferedInputStream.  This is read one line at a time
 * returning the results as a List<String>.   It is used to process output
 * from tiledbvcf list --uri <uri> and potentially other commands.
 */
fun inputStreamProcessing(tiledbList: BufferedInputStream): List<String> {
    val samples = mutableListOf<String>()  // this is the list of samples to be returned

    try {
        tiledbList.bufferedReader().use { br ->
            var line = br.readLine()
            while (line != null) {
                // skip blank lines.
                if (line != "") samples.add(line)
                line = br.readLine()
            }
        }
    } catch (exc: Exception) {
        myLogger.error("Error reading tiledb list output: ${exc.message}")
    } finally {
        tiledbList.close()
    }
    return samples
}

/**
 * This function takes a range map of Position, plus a new Range<Position> and gene string, and merges the
 * new range if it overlaps with an existing range.  If it does not overlap, it adds the new range to the map.
 * For merged ranges, the gene string is appended to the existing gene string.
 */
fun addRange(geneRange: RangeMap<Position, String>, range: Range<Position>, gene:String) {
    val overlaps: List<Map.Entry<Range<Position>, String>> =
        ArrayList<Map.Entry<Range<Position>, String>>(
            geneRange.subRangeMap(range).asMapOfRanges().entries
        )

    //if overlaps has length, merge ranges together
    if (overlaps.size != 0) {
        val overlappingEntry: Map.Entry<Range<Position>, String>? = geneRange.getEntry(
            overlaps[0].key.lowerEndpoint()
        )
        //then use the combined range and assign the call
        val newGene = "${overlappingEntry!!.value}-$gene"

        // Update overlappingEntry value with new merged gene value.
        // 2nd put is to ensure new entry is merged with the new value
        geneRange.put(overlappingEntry.key, newGene)
        geneRange.putCoalescing(range, newGene)
    } else {
        geneRange.put(range, gene)
    }
}

/**
 * This functions takes a RangeMap of Positions, plus a "pad" value, and returns a new RangeMap with the
 * ranges expanded by the pad value.  The reference sequence is used to determine chromosome sizes.
 * If there are not enough bases between ranges to add the specified padding, the existing base pairs
 * are split between the ranges, and padding is cut to what is available.  If the padding will exceed
 * the end of the chromsome, it is truncated as appropriate.
 */
fun createFlankingList(geneRange:RangeMap<Position,String>, numFlanking:Int, refSeq: Map<String, NucSeq>): RangeMap<Position,String> {

    // there are no embedded or overlapped entries in the geneRange map.
    // The goal is to add "numFlanking" flanking to start/end of each entry.
    // There may not be 2*numFlanking between each gene entry.  If not, split the difference
    // between:  half to geneA end, half to geneB start.  This is done in
    // findFlankingStartPos() and findFlankingEndPos()
    val flankingRange: RangeMap<Position, String> = TreeRangeMap.create()

    try {
        geneRange.asMapOfRanges().entries.forEach { range ->
            val data:Range<Position> = range.key
            val chr = data.lowerEndpoint().contig
            if (!(refSeq.keys.contains(chr))) {
                throw IllegalArgumentException("createFlankingList: chrom $chr not found in reference fasta.")
            }
            val chrLen = refSeq[chr]!!.size()

            // Find new start/end positions with specified number of flanking, add this position to the map
            var flankingStart = findFlankingStartPos(geneRange,data, numFlanking)
            // adjustment for 0-based occurred when creating GeneRanges
            flankingStart = Position(flankingStart.contig,flankingStart.position) // start was adjusted to be 0-based in geneRange
            val flankingEnd = findFlankingEndPos(geneRange,data, numFlanking, chrLen)
            flankingRange.put(Range.closed( flankingStart, flankingEnd), range.value)

        }
        if (geneRange.asMapOfRanges().size != flankingRange.asMapOfRanges().size) {
            println(
                "ERROR - original gene list size ${geneRange.asMapOfRanges().size} does NOT equal created flanking gene list size ${flankingRange.asMapOfRanges().size}"
            )
            throw IllegalArgumentException("ERROR - original gene list size ${geneRange.asMapOfRanges().size} does NOT equal created flanking gene list size ${flankingRange.asMapOfRanges().size}")
        }
        return flankingRange
    } catch (exc:Exception) {
        exc.printStackTrace()
        throw IllegalArgumentException("ERROR - createFlankingList faulted with message: ${exc.message}")
    }
}

/**
 * Using the geneRange map, the current range and a specified number of flanking bps,
 * this function returns the new start position for the range.  If there are not enough
 * bps between the current range and the previous range, the start position is adjusted
 * to split the difference between the two ranges.  If the start position is less than 0,
 * it is adjusted to 0.
 */
fun findFlankingStartPos(geneRange:RangeMap<Position,String>, data:Range<Position>, numFlanking:Int):Position{

    val flankCheck = numFlanking*2
    val chrom = data.lowerEndpoint().contig
    val curStart = data.lowerEndpoint().position
    val lowerFlank = data.lowerEndpoint().position - flankCheck // at this point, is ok if this is negative
    val flankStartTop = Position(chrom,lowerFlank)
    val flankUpToLower = Position(chrom,curStart-1)
    val startCheck = Range.closed(flankStartTop,flankUpToLower) // will be used to find overlaps

    // Want to add (user defined number) flanking on either side of each entry in geneRange map,
    // then add the new range positions to the flankingRange map to be returned.
    // Before adding the flanking bps, need to verify
    //  1.  there are necessary number bps between the start of the chrom and the start of this entry (or start = 1)
    //  2.  the distance between this entry and the current entry is >= numFlanking*2
    //  3.  the distance between this entry and the chrom len is at least numFlanking (otherwise chromlen = end)
    val overlapsStart:List<Map.Entry<Range<Position>,String>> =
        ArrayList<Map.Entry<Range<Position>, String>>(geneRange.subRangeMap(startCheck).asMapOfRanges().entries)
    var newLowerPos:Position? = null

    if (overlapsStart.size > 0) {
        // the last overlap should be the closest in position to the current range start
        val prevEnd = overlapsStart.get(overlapsStart.size-1).key.upperEndpoint().position
        val newFlankNum = curStart-prevEnd // flanking is adjusted to be half on either side of the previous range
        val newLowerInt = (curStart - newFlankNum/2)+1 // want start to be 1 past the end of the previous range
        newLowerPos = Position(chrom,newLowerInt)
    } else {
        val newLowerInt = if (curStart - numFlanking < 1)  1 else curStart - numFlanking
        newLowerPos = Position(chrom,newLowerInt)
    }
    return newLowerPos
}

/**
 * Using the geneRange map, the current range and a specified number of flanking bps,
 * this function returns the new end position for the range.  If there are not enough
 * bps between the current range and the next range, the end position is adjusted
 * to split the difference between the two ranges.  If the end position is greater than
 * the chromosome length, it is adjusted to the chromosome length.
 */
fun findFlankingEndPos(geneRange:RangeMap<Position,String>, data:Range<Position>, numFlanking:Int, chromLen:Int):Position{
    val chrom = data.lowerEndpoint().contig
    val flankCheck = numFlanking*2
    val curEnd = data.upperEndpoint().position
    val upperflank = data.upperEndpoint().position + flankCheck
    val curEndPlusOne = Position(chrom,curEnd+1)
    val flankEndTop = Position(chrom,upperflank)
    val endCheck = Range.closed(curEndPlusOne,flankEndTop)
    val overlapsEnd: List<Map.Entry<Range<Position>, String>> =
        ArrayList<Map.Entry<Range<Position>, String>>(geneRange.subRangeMap(endCheck).asMapOfRanges().entries)
    var newUpperPos:Position? = null

    if (overlapsEnd.size > 0) {
        // the first overlap should be the closest in position to the current range end
        val nextStart = overlapsEnd.get(0).key.lowerEndpoint().position
        val newFlankNum = nextStart - curEnd // adjust flanking to be half on either side of the next range
        val newUpperInt = (curEnd + newFlankNum / 2) // don't add 1 here, gets added at start above
        newUpperPos = Position(chrom, newUpperInt)
    } else {
        val newUpperInt = if (curEnd + numFlanking > chromLen) chromLen else curEnd + numFlanking
        newUpperPos = Position(chrom, newUpperInt)
    }
    return newUpperPos
}