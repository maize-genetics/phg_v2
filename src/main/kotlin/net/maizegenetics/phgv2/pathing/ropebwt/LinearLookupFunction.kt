package net.maizegenetics.phgv2.pathing.ropebwt

import com.google.common.collect.RangeMap
import net.maizegenetics.phgv2.utils.Position
import kotlin.math.abs

class LinearLookupFunction(val knotMap: RangeMap<Position, Pair<Position, Position>>) {
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
        asmRange : Pair<Position, Position>,
        referenceRange : Pair<Position, Position>
    ): Position {

        val offsetProp = (asmPosition.position - asmRange.first.position).toDouble() /
                (asmRange.second.position - asmRange.first.position).toDouble()

        val offsetPos = abs(referenceRange.second.position - referenceRange.first.position) * offsetProp

        if(referenceRange.first.contig != referenceRange.second.contig) {
            //If the contigs are not the same, we cannot interpolate
            Position("unknown", 0)
        }
        return if(referenceRange.first.position == referenceRange.second.position) {
            //If the positions are the same, all values are equal
            Position(referenceRange.first.contig, referenceRange.first.position)
        }
        else if(referenceRange.first.position < referenceRange.second.position) {
            //Positive strand
            val newPos = referenceRange.first.position + offsetPos
            Position(referenceRange.first.contig, newPos.toInt())
        } else {
            //Negative strand
            val newPos = referenceRange.first.position - offsetPos
            Position(referenceRange.first.contig, newPos.toInt())
        }
    }
}