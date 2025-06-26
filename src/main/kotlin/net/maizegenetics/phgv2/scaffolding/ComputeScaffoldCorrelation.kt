package net.maizegenetics.phgv2.scaffolding

import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.convert
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.double
import org.ejml.data.DMatrixSparseCSC
import org.ejml.data.DMatrixSparseTriplet
import org.ejml.ops.DConvertMatrixStruct
import java.io.File
import kotlin.math.abs
import kotlin.math.sqrt


enum class CORRELATION_TYPE {
    PEARSON,
    SPEARMAN
}

data class ScaffoldAlignmentMatrix(val gbsRowIdxMap: Map<String, Int>,
                                   val unitigColIdxMap: Map<String, Int>,
                                   val matrix: DMatrixSparseCSC)

class ComputeScaffoldCorrelation: CliktCommand(help = "Compute scaffold correlation") {

    val collectedCountsDir by option(help = "Input directory containing collected counts")
        .required()

    val outputMatrix by option(help = "Output file for scaffold correlation results")
        .required()

    val correlationThreshold by option(help = "Correlation threshold for filtering results")
        .double()
        .default(0.90)

    val correlationType by option(help = "Type of correlation to compute (PEARSON or SPEARMAN)")
        .convert { CORRELATION_TYPE.valueOf(it.uppercase()) }
        .default(CORRELATION_TYPE.PEARSON)

    override fun run() {
        // This command is a placeholder for computing scaffold correlation.
        // The actual implementation would involve reading scaffold data,
        // calculating correlations, and outputting the results.

        println("Computing scaffold correlation...")
        // Implement the logic to compute scaffold correlation here
        computeScaffoldCorrelation(collectedCountsDir, outputMatrix, correlationThreshold, correlationType)
    }

    fun computeScaffoldCorrelation(collectedCountsDir: String, outputMatrixFile: String, correlationThreshold: Double = 0.01, correlationType: CORRELATION_TYPE) {
        // Logic to read the collected counts from the directory,
        // compute scaffold correlations, and write the results to the output matrix file.

        println("Collected Counts Directory: $collectedCountsDir")
        println("Output Matrix File: $outputMatrixFile")

        // Implement the actual computation logic here
        val alignmentMatrix = buildAlignmentCountMatrix(collectedCountsDir)
        println("Computing correlations")
        val correlationMatrix = computeCorrelationIncludingZeros(alignmentMatrix.matrix, correlationThreshold, correlationType)

        println("Reversing index maps for output")
        val gbsIdxToName = reverseMatrixMap(alignmentMatrix.gbsRowIdxMap)
        val unitigIdxToName = reverseMatrixMap(alignmentMatrix.unitigColIdxMap)

        println("Writing out correlation matrix to $outputMatrixFile")
        writeOutCorrelationMatrix(correlationMatrix, gbsIdxToName, unitigIdxToName, outputMatrixFile)
    }

    fun buildAlignmentCountMatrix(collectedCountsDir: String): ScaffoldAlignmentMatrix {

        //Load all the files into data points : Triple<gbsSampleName, unitigId, count>
        println("Parsing Files in Directory: $collectedCountsDir")
        val triples = File(collectedCountsDir).listFiles().flatMap { file ->
            parseFileIntTriples(file)
        }

        println("Building dimension maps")
        //Build maps for each dimension so we can do an efficient lookup.
        val (gbsSampleToIdxMap, unitigToIdxMap) = buildDimensionMaps(triples)

        println("Building alignment matrix")
        val triplet = DMatrixSparseTriplet(gbsSampleToIdxMap.size, unitigToIdxMap.size, triples.size)
        for ((r, c, v) in triples) {
            val rIdx = gbsSampleToIdxMap[r] ?: error("GBS sample $r not found in index map")
            val cIdx = unitigToIdxMap[c] ?: error("Unitig $c not found in index map")
            triplet.addItem(rIdx, cIdx, v.toDouble())
        }
        return ScaffoldAlignmentMatrix(gbsSampleToIdxMap,
            unitigToIdxMap,
            DConvertMatrixStruct.convert(triplet, null as DMatrixSparseCSC?))
    }


    fun parseFileIntTriples(file: File): List<Triple<String,String,Int>> {
        val lines = bufferedReader("${file.path}").readLines().map { it.split("\t") }
        val header = lines.first()
        val columnIndicesMap = header.withIndex().associate { it.value to it.index }

        return lines.drop(1).map { line ->
            val gbsSampleName = line[columnIndicesMap["gbsSampleName"] ?: 0]
            val assemblySampleName = line[columnIndicesMap["assemblySampleName"] ?: 1]

            val unitigId = line[columnIndicesMap["unitigId"] ?: 2]
            val count = line[columnIndicesMap["count"] ?: 3].toIntOrNull() ?: 0
            Triple(gbsSampleName, "${unitigId}_${assemblySampleName}", count)
        }
    }

    fun buildDimensionMaps(triples: List<Triple<String,String,Int>>): Pair<Map<String,Int>, Map<String,Int>> {
        val gbsSampleToIdxMap = mutableMapOf<String, Int>()
        val unitigToIdxMap = mutableMapOf<String, Int>()

        var gbsSampleIndex = 0
        var unitigIndex = 0

        for ((gbsSampleName, unitigId, _) in triples) {
            if (!gbsSampleToIdxMap.containsKey(gbsSampleName)) {
                gbsSampleToIdxMap[gbsSampleName] = gbsSampleIndex++
            }
            if (!unitigToIdxMap.containsKey(unitigId)) {
                unitigToIdxMap[unitigId] = unitigIndex++
            }
        }
        return Pair(gbsSampleToIdxMap, unitigToIdxMap)
    }


    // Extract full column with zeros included
    fun extractFullColumn(matrix: DMatrixSparseCSC, col: Int): DoubleArray {
        val result = DoubleArray(matrix.numRows) { 0.0 }
        val start = matrix.col_idx[col]
        val end = matrix.col_idx[col + 1]
        for (i in start until end) {
            val row = matrix.nz_rows[i]
            result[row] = matrix.nz_values[i]
        }
        return result
    }

    // Pearson correlation between two full vectors
    fun pearsonCorr(x: DoubleArray, y: DoubleArray): Double {
        val meanX = x.average()
        val meanY = y.average()

        var num = 0.0
        var denX = 0.0
        var denY = 0.0

        for (i in x.indices) {
            val dx = x[i] - meanX
            val dy = y[i] - meanY
            num += dx * dy
            denX += dx * dx
            denY += dy * dy
        }

        val denom = sqrt(denX * denY)
        return if (denom != 0.0) num / denom else 0.0
    }

    fun rank(values: DoubleArray): DoubleArray {
        val indexed = values.mapIndexed { i, v -> i to v }
        val sorted = indexed.sortedBy { it.second }

        val ranks = DoubleArray(values.size)
        var i = 0
        while (i < sorted.size) {
            var j = i
            // Move j forward while there are ties
            while (j + 1 < sorted.size && sorted[j + 1].second == sorted[i].second) {
                j++
            }

            // Compute average rank for all tied entries
            val avgRank = (i..j).map { it + 1.0 }.average()
            for (k in i..j) {
                ranks[sorted[k].first] = avgRank
            }

            // Move to the next group
            i = j + 1
        }
        return ranks
    }

    // Spearman correlation between two vectors
    fun spearmanCorr(x: DoubleArray, y: DoubleArray): Double {
        val rx = rank(x)
        val ry = rank(y)

        val meanX = rx.average()
        val meanY = ry.average()

        var num = 0.0
        var denX = 0.0
        var denY = 0.0

        for (i in rx.indices) {
            val dx = rx[i] - meanX
            val dy = ry[i] - meanY
            num += dx * dy
            denX += dx * dx
            denY += dy * dy
        }

        val denom = sqrt(denX * denY)
        return if (denom != 0.0) num / denom else 0.0
    }

    // Compute full Pearson correlation matrix
    fun computeCorrelationIncludingZeros(matrix: DMatrixSparseCSC, correlationThreshold: Double = .01, type:CORRELATION_TYPE): DMatrixSparseCSC {
        val cols = matrix.numCols
        val resultTriples = mutableListOf<Triple<Int,Int,Double>>()

        val fullCols = Array(cols) { extractFullColumn(matrix, it) }

        for (i in 0 until cols) {
            for (j in i until cols) {
                val r = if(type == CORRELATION_TYPE.PEARSON) {
                    pearsonCorr(fullCols[i], fullCols[j])
                } else {
                    spearmanCorr(fullCols[i], fullCols[j])
                }
                if(abs(r) > correlationThreshold) {
                    resultTriples.add(Triple(i, j, r))
                }
            }
        }

        val resultTriplet = DMatrixSparseTriplet(cols, cols, resultTriples.size)
        for ((r, c, v) in resultTriples) {
            resultTriplet.addItem(r, c, v)
            if (r != c) {
                resultTriplet.addItem(c, r, v) // Ensure symmetry
            }
        }
        return DConvertMatrixStruct.convert(resultTriplet, null as DMatrixSparseCSC?)
    }

    fun writeOutCorrelationMatrix(correlationMatrix: DMatrixSparseCSC, gbsRowIdxMap: Map<Int,String>, unitigColIdxMap: Map<Int,String>, outputFile: String) {
        bufferedWriter(outputFile).use { writer ->

            writer.write("gbsSampleName\tunitigId\tcorrelation\n")
            val iterator = correlationMatrix.createCoordinateIterator()
            while (iterator.hasNext()) {
                val currentTriple = iterator.next()
                val rowIdx = currentTriple.row
                val colIdx = currentTriple.col
                val value = currentTriple.value

                val unitigId1 = unitigColIdxMap[rowIdx] ?: error("Unitig ID not found for index  $rowIdx")
                val unitigId2 = unitigColIdxMap[colIdx] ?: error("Unitig ID not found for index $colIdx")
                if (unitigId1.isNotEmpty() && unitigId2.isNotEmpty()) {
                    writer.write("$unitigId1\t$unitigId2\t$value\n")
                }
            }
        }
    }

    fun reverseMatrixMap( map: Map<String, Int>): Map<Int, String> {
        return map.entries.associate { (key, value) -> value to key }
    }

}