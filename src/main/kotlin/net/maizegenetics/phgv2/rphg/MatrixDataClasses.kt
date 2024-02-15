package net.maizegenetics.phgv2.rphg

typealias IntMatrix = MatrixWithNames<IntArray>
typealias DoubleMatrix = MatrixWithNames<DoubleArray>
typealias StringMatrix = MatrixWithNames<Array<String>>
typealias DataFrame = MatrixWithNames<Array<out Any>>

data class MatrixWithNames<T>(
    var colNames: Array<String>? = null,
    var rowNames: Array<String>? = null,
    var matrixData: Array<T>? = null
) {
    override fun toString(): String {
        if (matrixData?.isArrayOf<IntArray>() == true) {
            return "Int array of ${matrixData?.get(0)}"
        }
        return "Number of columns: ${matrixData?.get(0)}"
    }
}