package net.maizegenetics.phgv2.pathing.ropebwt

import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import com.google.common.collect.TreeRangeMap
import net.maizegenetics.phgv2.utils.Position
import kotlin.math.abs

//class LinearLookupFunction(val knotMap: RangeMap<Position, Pair<Position, Position>>) {
class LinearLookupFunction( knots: Map<String,Pair<DoubleArray,DoubleArray>>) {

    var knotMap : RangeMap<Position,Pair<Position,Position>>

    init {
        knotMap = convertKnotsToLinearSpline(knots)
    }

    fun convertKnotsToLinearSpline(knots: Map<String,Pair<DoubleArray,DoubleArray>>) : RangeMap<Position,Pair<Position,Position>> {
        //Need to make the splines into linear blocks

        //loop through each assembly chromosome:
        val knotMap = TreeRangeMap.create<Position,Pair<Position,Position>>()
        var counter = 0
        knots.forEach { (key, value) ->
            println("Processing $key: ${counter++}/${knots.size-1}")

            val asmPositions = value.first
            val refPositions = value.second

            if(asmPositions.size != refPositions.size) {
                throw IllegalStateException("ASM and REF positions are not the same size for $key")
            }

            //We need to create a range for each pair of asm and ref positions
            for (i in 0 until asmPositions.size - 1) {
                val asmStart = Position(key, asmPositions[i].toInt())
                val asmEnd = Position(key, asmPositions[i + 1].toInt())

                //Convert the ref positions to Position objects as they are encoded
                val refStartPos = PS4GUtils.decodePosition(refPositions[i].toInt())
                val refEndPos = PS4GUtils.decodePosition(refPositions[i + 1].toInt())

                if(refStartPos.contig != refEndPos.contig) {
                    continue // Skip if the contigs are not the same  We do not want a spline between them
                }

                knotMap.put(Range.closed(asmStart, asmEnd), Pair(refStartPos, refEndPos))
            }
        }
        return knotMap
    }

    fun value(position: Position): Position {
        //Find the range that contains the position

        val range = knotMap.getEntry(position)
        return if(range != null) {
            val asmRange = range.key
            val asmSt = asmRange.lowerEndpoint()
            val asmEnd = asmRange.upperEndpoint()
            resolveLinearInterpolation(position,Pair(asmSt,asmEnd), range.value)
        } else {
            Position("unknown", 0) // Return a default position if not found
        }
    }

    fun resolveLinearInterpolation(
        asmPosition: Position,
        asmRegion : Pair<Position, Position>,
        referenceRegion : Pair<Position, Position>
    ): Position {

        val offsetProp = (asmPosition.position - asmRegion.first.position).toDouble() /
                (asmRegion.second.position - asmRegion.first.position).toDouble()

        val offsetPos = abs(referenceRegion.second.position - referenceRegion.first.position) * offsetProp

        return if(referenceRegion.first.contig != referenceRegion.second.contig) {
            //If the contigs are not the same, we cannot interpolate
            Position("unknown", 0)
        }
        else if(referenceRegion.first.position == referenceRegion.second.position) {
            //If the positions are the same, all values are equal
            Position(referenceRegion.first.contig, referenceRegion.first.position)
        }
        else if(referenceRegion.first.position < referenceRegion.second.position) {
            //Positive strand
            val newPos = referenceRegion.first.position + offsetPos
            Position(referenceRegion.first.contig, newPos.toInt())
        } else {
            //Negative strand
            val newPos = referenceRegion.first.position - offsetPos
            Position(referenceRegion.first.contig, newPos.toInt())
        }
    }
}