package net.maizegenetics.phgv2.pathing

import java.io.File
import java.nio.file.Files

/**
 * Parses a key file into lines.
 * @param keyFile   a keyFile
 * @param expectedColumnNames   the names of the columns expected in this file.
 * This is used to check that the expected columns are present.
 * @param expectedUniqueColumns a list of the columns that are expected to be jointly unique.
 *
 * For example, if a [expectedUniqueColumns] is a list of two names then the string formed by combining
 * those columns is expected to be unique. When the class is instantiated, the column headers are checked
 * to make sure the expected column names are present. The file is also checked for duplicates. If either check fails an
 * [IllegalArgumentException] is thrown.
 *
 * @property keyFileLines   a list of the keyFile lines, parsed into individual columns, without the header
 * @property columnNames    a list of the keyFile column names
 * @property columnNameMap  a map of column name -> index into each keyFile line
 *
 * For example, a sample name can be retrieved from keyFileLines using keyFileLines[lineIndex][columnNameMap["SampleName"]].
 *
 */
class KeyFile(val keyFile: File, val expectedColumnNames: List<String>, val expectedUniqueColumns: List<String> = listOf<String>()) {
    val keyFileLines: List<List<String>>
    val columnNames: List<String>
    val columnNameMap: Map<String,Int>

    init {
        val lines = Files.newBufferedReader(keyFile.toPath()).readLines()
        columnNames = lines[0].split("\t")
        keyFileLines = lines.drop(1).map { it.split("\t") }
        columnNameMap = columnNames.mapIndexed { index, name -> Pair(name, index) }.toMap()
        checkRequiredColumns()
        checkUniqueColumns()
    }

    private fun checkRequiredColumns() {
        for (expectedName in expectedColumnNames) {
            if (!columnNames.contains(expectedName)) {
                val expected = expectedColumnNames.joinToString(",")
                throw IllegalArgumentException("Expected columns named $expected but not all were present")
            }
        }
    }
    
    private fun checkUniqueColumns() {
        if (expectedUniqueColumns.isEmpty()) return
        val uniqueIndices = expectedUniqueColumns.map { columnNameMap[it] }.filter { it != null }
        if (uniqueIndices.size < expectedUniqueColumns.size) {
            throw IllegalArgumentException("list of expectedUniqueColumns was " +
                    "${expectedUniqueColumns.joinToString(",")} but not all of them were key file columns.")
        }

        val valuesThatShouldBeUnique = keyFileLines.map {line ->
            //none of the indices = null since null values were filtered out
            uniqueIndices.map{ index -> line[index!!] ?: "NULL" }.joinToString(",")
        }

        val nonUniqueValues = valuesThatShouldBeUnique.groupingBy { it }.eachCount().filter { (_, count) -> count > 1 }
        if (nonUniqueValues.isNotEmpty()) {
            val uniqueHeader = expectedUniqueColumns.joinToString(",")
            val errorList = nonUniqueValues.entries.map { (value, count) -> "$value ($count)"}
            val errorMessage = "Duplicates found in key file for $uniqueHeader (count)\n${errorList.joinToString("\n")}"
            throw IllegalArgumentException(errorMessage)
        }
    }

}