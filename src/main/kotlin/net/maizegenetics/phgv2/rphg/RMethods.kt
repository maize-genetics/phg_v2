package net.maizegenetics.phgv2.rphg

import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.SampleGamete

class RMethods {
    fun hapIdMatrix(g: HaplotypeGraph): MutableList<Pair<SampleGamete, String>> {
        val allRanges = g.ranges()

        val returnState = mutableListOf<Pair<SampleGamete, String>>()
        allRanges.forEach { range ->
            val hapIdRange = g.sampleGameteToHaplotypeId(range)
            hapIdRange.forEach {
                returnState.add(Pair(it.key, it.value))
            }
        }

        return returnState
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

    val x = arrayOf(
        arrayOf(1, 2, 3),
        arrayOf("this", "is", "category")
    )
}

fun main() {
    val graph = HaplotypeGraph(
        listOf(
            "/Users/bm646-admin/Projects/phg_v2/data/test/smallseq/LineA.h.vcf",
            "/Users/bm646-admin/Projects/phg_v2/data/test/smallseq/LineB.h.vcf",
            "/Users/bm646-admin/Projects/phg_v2/data/test/smallseq/Ref.h.vcf"
        )
    )

    println(graph.numberOfRanges())
    println(graph.samples())

}