package net.maizegenetics.phgv2.pathing.ropebwt

import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import com.google.common.collect.TreeRangeMap
import it.unimi.dsi.fastutil.longs.Long2LongAVLTreeMap
import it.unimi.dsi.fastutil.longs.Long2LongSortedMap
import net.maizegenetics.phgv2.utils.Position
import kotlin.math.abs

class LinearLookupFunction(knots: Map<String,List<Triple<Int,String,Int>>>, refChrIndexMap : Map<String,Int>) {

    val knotMap : Long2LongSortedMap
    val asmChrIndexMap : Map<String,Int>
    val refChrIndexMap : Map<String,Int>
    val indexRefChrMap : Map<Int,String> //Need this to do the reverse lookup to build a position to return

    init {
        //I tried to do this but it didn't work:
        //val (knotMap, asmChrIndexMap) = convertKnotsToLinearSplinePrimitiveAVL(knots, refChrIndexMap)
        val knotMapAndAsmChrMap = convertKnotsToLinearSplinePrimitiveAVL(knots, refChrIndexMap)
        knotMap = knotMapAndAsmChrMap.first
        asmChrIndexMap = knotMapAndAsmChrMap.second
        this.refChrIndexMap = refChrIndexMap
        this.indexRefChrMap = refChrIndexMap.entries.associate { (k, v) -> v to k }
    }


    /**
     * Function to setup the linear spline mapping using a Long2LongAVLTreeMap for efficiency.
     */
    fun convertKnotsToLinearSplinePrimitiveAVL(knots: Map<String,List<Triple<Int,String,Int>>>, refChrIndexMap: Map<String,Int>) : Pair<Long2LongAVLTreeMap,Map<String,Int>>{
        val treeMap = Long2LongAVLTreeMap()
        val asmChrIdLookup = mutableMapOf<String,Int>()
        knots.forEach { (key, value) ->
            for(triple in value) {
                if(!asmChrIdLookup.containsKey(key)) {
                    asmChrIdLookup[key] = asmChrIdLookup.size
                }

                val asmChromId = asmChrIdLookup[key]!!
                val refChromId = refChrIndexMap[triple.second]!!

                //Maybe update this to leave more bits for the position?
                val asmPosLong = (asmChromId.toLong() shl 32) or triple.first.toLong()
                val refPosLong = (refChromId.toLong() shl 32) or triple.third.toLong()

                treeMap[asmPosLong] = refPosLong
            }
        }

        return Pair(treeMap, asmChrIdLookup)
    }


    /**
     * Function to get the mapped position for a given input position.
     */
    fun value(position: Position): Position {
        if(!asmChrIndexMap.containsKey(position.contig)) {
            return Position("unknown", 0)
        }
        val encoded = (asmChrIndexMap[position.contig]!!.toLong() shl 32) or position.position.toLong()

        val floorKey = if (knotMap.containsKey(encoded)) {
            encoded
        } else {
            val headMap = knotMap.headMap(encoded)
            if(headMap.isEmpty() || headMap == null) {
                return Position("unknown", 0)
            }
            headMap.lastLongKey()
        }


        val ceilingKey = if( knotMap.containsKey(encoded)) {
            encoded
        }
        else {
            val tailMap = knotMap.tailMap(encoded)
            if(tailMap.isEmpty() || tailMap == null) {
                return Position("unknown", 0)
            }
            tailMap.firstLongKey()
        }


        //check that the chromosomes match
        val floorASMChromId = (floorKey shr 32).toInt()
        val ceilingASMChromId = (ceilingKey shr 32).toInt()
        if(floorASMChromId != ceilingASMChromId) {
            return Position("unknown", 0)
        }

        //KeyPos variables are the positions for the input assembly
        val floorKeyPos = floorKey and 0xFFFFFFFF
        val ceilingKeyPos = ceilingKey and 0xFFFFFFFF

        val floorEncodedValue = knotMap[floorKey]
        val ceilingEncodedValue = knotMap[ceilingKey]

        val floorChromId = (floorEncodedValue shr 32).toInt()
        val ceilingChromId = (ceilingEncodedValue shr 32).toInt()
        if(floorChromId != ceilingChromId) {
            return Position("unknown", 0)
        }
        //get the positions and remove the chrom idx bits
        val floorRefPos = floorEncodedValue and 0xFFFFFFFF
        val ceilingRefPos = ceilingEncodedValue and 0xFFFFFFFF

        val chromName = indexRefChrMap[floorChromId]

        if(floorRefPos == ceilingRefPos) {
            //If the positions are the same, all values are equal
            return Position(chromName!!,(floorRefPos and 0xFFFFFFFF).toInt())
        }


        val offsetProp = (position.position - floorKeyPos).toDouble() /
                (ceilingKeyPos - floorKeyPos).toDouble()
        val offsetPos = abs(ceilingRefPos - floorRefPos) * offsetProp


        val finalPosition = if(floorRefPos < ceilingRefPos) {
            //Positive strand
            (floorRefPos + offsetPos.toLong())
        } else {
            //Negative strand
            (floorRefPos - offsetPos.toLong())
        }


        return Position(chromName!!, (finalPosition and 0xFFFFFFFF).toInt())
    }

}