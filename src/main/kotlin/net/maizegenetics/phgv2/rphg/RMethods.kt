package net.maizegenetics.phgv2.rphg

import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.SampleGamete
import javax.xml.crypto.Data

class RMethods {
    fun getFullHapIdMatrix(g: HaplotypeGraph): StringMatrix {

        val allRanges = g.ranges()
        val sampleGametes = g.sampleGametesInGraph()
        val array2D = Array(g.samples().size) { Array(allRanges.size) { "" } }

        allRanges.forEach { range ->
            val hapIdRange = g.sampleGameteToHaplotypeId(range)
            hapIdRange.forEach {
                array2D[sampleGametes.indexOf(it.key)][allRanges.indexOf(range)] = it.value
            }
        }

        return StringMatrix (
            colNames = allRanges.map { "R$it" }.toTypedArray(),
            rowNames = sampleGametes.map { "${it.name}_G${it.gameteId + 1}" }.toTypedArray(),
            matrixData = array2D
        )
    }

    fun getRefRanges(g: HaplotypeGraph): DataFrame {
        val refRanges = g.ranges()
        return DataFrame(
            colNames = arrayOf("seqname", "start", "end"),
            rowNames = null,
            matrixData = arrayOf(
                refRanges.map { it.contig }.toTypedArray(),
                refRanges.map { it.start }.toTypedArray(),
                refRanges.map { it.end }.toTypedArray()
            )
        )
    }
    
    fun testStringMatrix(): StringMatrix {
        return MatrixWithNames(
            colNames = arrayOf("Col1", "Col2", "Col3"),
            rowNames = arrayOf("Row1", "Row2", "Row3", "Row4"),
            matrixData = arrayOf(
                arrayOf("a", "b", "c"),
                arrayOf("d", "e", "f"),
                arrayOf("g", "h", "i"),
                arrayOf("j", "k", "l")
            )
        )
    }

    fun testIntMatrix(): IntMatrix {
        return MatrixWithNames (
            arrayOf("i1", "i2", "i3", "i4"),
            arrayOf("r1", "r2"),
            arrayOf(
                intArrayOf(1, 2, 3, 4),
                intArrayOf(5, 6, 7, 8)
            )
        )
    }

    fun testDblMatrix(): DoubleMatrix {
        return MatrixWithNames(
            arrayOf("d1", "d2", "d3", "d4"),
            arrayOf("r1", "r2"),
            arrayOf(
                doubleArrayOf(1.3, 2.1, 3.1434, 4.2),
                doubleArrayOf(5.1, 6.234, 7.0, 8.67)
            )
        )
    }

    fun testDataFrame(): DataFrame {
        return MatrixWithNames(
            arrayOf("df1", "df2", "df3", "df4"),
            arrayOf("r1", "r2"),
            arrayOf(
                arrayOf(1, 2, 3),
                arrayOf("this", "is", "category")
            )
        )
    }
}

