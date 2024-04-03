package net.maizegenetics.phgv2.rphg

/**
 * General data class for storing simple data collections for
 * quick evaluation on the R side for rPHG2. Can be further modified
 * by specifying the type (`<T>`) (see type alias block).
 *
 * The three main R data collection objects are:
 *   * matrices - 2d array of one data type (e.g. character, numeric,
 *     integer, logical)
 *   * lists - a key-value map of name (key) and vector (value) -
 *     each key-value pair can be any basic type (see prior point)
 *   * data frames - modified list objects that are represented in
 *     tabular format. Each "column" can be of different data type
 *
 * Type aliases are used for simplifying the coding for type of data
 * returned for each possible R data object in the foreseeable future:
 *   * `RIntegerMatrix` - an R matrix consisting of type `integer`
 *      elements (`IntArray` in Kotlin)
 *   * `RNumericMatrix` - an R matrix consisting of type `numeric`
 *     elements (`DoubleArray` in Kotlin)
 *   * `RCharacterMatrix` - an R matrix consisting of type `character`
 *     elements (`Array<String>` in Kotlin)
 *   * `RList` - an R list consisting of a collection of variable
 *     data types (this can be easily converted into a dataframe)
 *
 * Additionally, arrays of column and row names (`colNames` and
 * `rowNames`) are added to the constructor to specify the names of
 * the dimensions on the R side (data frames generally only need
 * column names, but matrices usually use both column and row names)
 */

typealias RIntegerMatrix = MatrixWithNames<IntArray>
typealias RNumericMatrix = MatrixWithNames<DoubleArray>
typealias RCharacterMatrix = MatrixWithNames<Array<String>>
typealias RList = MatrixWithNames<Array<out Any>>

data class MatrixWithNames<T>(
    var colNames: Array<String>? = null,
    var rowNames: Array<String>? = null,
    var matrixData: Array<T>? = null
) {
    /**
     * Specify the type of data returned - useful for validity
     * checking on the R side:
     *
     * ``` r
     * exampleData <- interface$getAltHeadersFromGraph(g)
     * exampleData$toString() == "phgv2_r_list"
     * ```
     */
    override fun toString(): String {
        return when {
            matrixData?.isArrayOf<IntArray>() == true -> {
                "phgv2_int_matrix"
            }
            matrixData?.isArrayOf<Array<String>>() == true -> {
                "phgv2_string_matrix"
            }
            matrixData?.isArrayOf<DoubleArray>() == true -> {
                "phgv2_dbl_matrix"
            }
            else -> {
                "phgv2_r_list"
            }
        }
    }
}