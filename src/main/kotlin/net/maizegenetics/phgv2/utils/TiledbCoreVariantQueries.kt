package net.maizegenetics.phgv2.utils

import io.tiledb.java.api.*
import io.tiledb.java.api.Array
import org.jetbrains.kotlinx.dataframe.DataFrame
import org.jetbrains.kotlinx.dataframe.api.*
import io.tiledb.java.api.Query

/**
 * This class contains queries to be run against a TileDB core arrays
 */

/**
 * This query is only for the variants array.  THis array should be populated with
 * hvcf data from both aligned asemblies and the imputation pipeline.
 * Queries to this array occur when the user has not specified RefChecksum or Regions
 * as attributes they would like returned.
 *
 * Initially returned a Map<String, List<Map<String, String>>> where the key is the RefRange
 * and the value is a list of maps where each map contains the fields SampleName, ID1, and ID2.
 * Converted this to a Kotlin DataFrame for easier processing and writing to a file.
 */
fun queryVariantsArray(
    arrayName: String,
    refRangeList: List<String>?,
    sampleList: List<String>?,
):DataFrame<*> {

    // TODO - determine correct output.  This currently returns
    // a map of a refRange to a list of maps, where each map contains
    // the fields SampleName, ID1, and ID2.
    val context = Context()
    //val result = mutableListOf<Map<String, String>>()
    val result = mutableMapOf<String, List<Map<String, String>>>()

    // Open the array in read mode
    val array = Array(context, arrayName, QueryType.TILEDB_READ)

    // Prepare buffers for RefRange - this will limit to the ranges in the user provided list
    val refRangeBuffers = if (refRangeList != null ) prepareVariableBuffer(context, refRangeList) else {
        Pair(io.tiledb.java.api.NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII),
            io.tiledb.java.api.NativeArray(context, 512, Datatype.TILEDB_UINT64))
    }

    // Prepare buffers for SampleName  - this will limit to the samples in the user provided list
    val sampleBuffers = if (sampleList != null) prepareVariableBuffer(context, sampleList) else {
        Pair(io.tiledb.java.api.NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII),
            io.tiledb.java.api.NativeArray(context, 512, Datatype.TILEDB_UINT64))
    }

    val id1Buffer = io.tiledb.java.api.NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII)
    val id1Offsets = io.tiledb.java.api.NativeArray(context, 512, Datatype.TILEDB_UINT64)

    val id2Buffer = io.tiledb.java.api.NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII)
    val id2Offsets = io.tiledb.java.api.NativeArray(context, 512, Datatype.TILEDB_UINT64)


    // create the query
    val query = Query(array, QueryType.TILEDB_READ).apply {
        setLayout(Layout.TILEDB_UNORDERED)

        // Set buffers for RefRange (dimension filter)
        setDataBuffer("RefRange", refRangeBuffers.first)
        setOffsetsBuffer("RefRange", refRangeBuffers.second)

        // Set buffers for SampleName (dimension retrieval)
        setDataBuffer("SampleName", sampleBuffers.first)
        setOffsetsBuffer("SampleName", sampleBuffers.second)

        // Set buffers for ID1/ID2 (attribute retrieval).  Offsets are not set
        // as these are fixed size buffers.
        setDataBuffer("ID1", id1Buffer)
        setDataBuffer("ID2", id2Buffer)

    }

    // Submit the query
    query.submit()

    // Use resultBufferElements to determine the number of valid entries
    val numSampleNames = query.resultBufferElements()["SampleName"]?.first?.toInt() ?: 0
    val numIDs = query.resultBufferElements()["ID1"]?.first?.toInt() ?: 0
    val numRefRanges = query.resultBufferElements()["RefRange"]?.first?.toInt() ?: 0
    println("LCJ - Num Results: SampleNames=$numSampleNames, IDs=$numIDs, RefRanges=$numRefRanges")

    // Extract SampleName data
    val sampleNameOffsetsArray = sampleBuffers.second.toJavaArray() as LongArray
    val sampleNameRawData = String(sampleBuffers.first.toJavaArray() as ByteArray)
    val sampleNames = sampleNameOffsetsArray.take(numSampleNames).mapIndexedNotNull { index, offset ->
        val end =
            if (index < sampleNameOffsetsArray.size - 1) {
                sampleNameOffsetsArray[index + 1].toInt()
            } else {
                sampleNameRawData.length  // last string ends at end of raw data
            }
        if (offset.toInt() >= 0 && offset.toInt() < end) {
            sampleNameRawData.substring(offset.toInt(), end).trimEnd('\u0000')
        } else {
            null // skip invalid offsets
        }

    }

    // Extract ID1 data
    val id1OffsetsArray = id1Offsets.toJavaArray() as LongArray
    val id1RawData = String(id1Buffer.toJavaArray() as ByteArray)

    val id1s = id1OffsetsArray.take(numIDs).mapIndexedNotNull { index, offset ->
        val end = if (index < numIDs - 1) {
            id1OffsetsArray[index + 1].toInt() // Next offset as end
        } else {
            id1RawData.length // Last string ends at the end of raw data
        }

        if (offset.toInt() >= 0 && offset.toInt() < end && end <= id1RawData.length) {
            id1RawData.substring(offset.toInt(), end).trimEnd('\u0000')
        } else {
            null // Skip invalid offsets
        }
    }

    // Extract ID2 data
    val id2OffsetsArray = id2Offsets.toJavaArray() as LongArray
    val id2RawData = String(id2Buffer.toJavaArray() as ByteArray)

    val id2s = id2OffsetsArray.take(numIDs).mapIndexedNotNull { index, offset ->
        val end = if (index < numIDs - 1) {
            id2OffsetsArray[index + 1].toInt() // Next offset as end
        } else {
            id2RawData.length // Last string ends at the end of raw data
        }

        if (offset.toInt() >= 0 && offset.toInt() < end && end <= id2RawData.length) {
            id2RawData.substring(offset.toInt(), end).trimEnd('\u0000')
        } else {
            null // Skip invalid offsets
        }
    }

    // Filter the results based on the refRanges.

    // Combine results grouped by RefRange
    val refRangeOffsetsArray = refRangeBuffers.second.toJavaArray() as LongArray
    val refRangeRawData = String(refRangeBuffers.first.toJavaArray() as ByteArray)
    val refRanges = refRangeOffsetsArray.take(numRefRanges).mapIndexedNotNull { index, offset ->
        val end = if (index < refRangeOffsetsArray.size - 1) refRangeOffsetsArray[index + 1].toInt() else refRangeRawData.length
        if (offset.toInt() < end) refRangeRawData.substring(offset.toInt(), end).trimEnd('\u0000') else null
    }

    refRanges.forEachIndexed { index, refRange ->
        val sampleName = sampleNames.getOrNull(index) ?: ""
        val id1 = id1s.getOrNull(index) ?: ""
        val id2 = id2s.getOrNull(index) ?: ""
        result[refRange] = result.getOrDefault(refRange, emptyList()) + mapOf(
            "SampleName" to sampleName,
            "ID1" to id1,
            "ID2" to id2
        )
    }

    // Close resources
    query.close()
    array.close()
    context.close()
     // print the result
    println("LCJ before converting to DataFrame- Result: $result")

    // convert to a dataFrame
    val dataFrame = convertQueryResultToDataFrame(result)
    //return result
    return dataFrame

}

/**
 * Convert the tiledb data to a Kotlin DataFrame for easier processing
 * and writing to a file.
 */
fun convertQueryResultToDataFrame(data: Map<String, List<Map<String, String>>>): DataFrame<*> {
    // Flatten the map into a list of rows
    val rows = data.flatMap { (refRange, records) ->
        records.map { record ->
            mapOf(
                "RefRange" to refRange,
                "SampleName" to record["SampleName"],
                "ID1" to record["ID1"],
                "ID2" to record["ID2"]
            )
        }
    }

    // Build the DataFrame
    return rows.toDataFrame()
}

/**
 * query distinct sample names from the filedb array.
 */
fun queryDistinctSampleNames(arrayName: String): Set<String> {
    val context = Context()

    // Open the TileDB array in read mode
    val array = Array(context, arrayName, QueryType.TILEDB_READ)

    // Prepare buffers for SampleName (variable-length strings)
    val sampleNameBuffer = NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII) // Adjust size as needed
    val sampleNameOffsets = NativeArray(context, 512, Datatype.TILEDB_UINT64) // Adjust size as needed

    // Set up the query
    val query = Query(array, QueryType.TILEDB_READ).apply {
        setLayout(Layout.TILEDB_GLOBAL_ORDER) // Read in global order
        setDataBuffer("SampleName", sampleNameBuffer)
        setOffsetsBuffer("SampleName", sampleNameOffsets)
    }

    // Submit the query
    query.submit()

    // Extract the number of results
    val numSampleNames = query.resultBufferElements()["SampleName"]?.first?.toInt() ?: 0

    // Extract SampleName data
    val offsetsArray = sampleNameOffsets.toJavaArray() as LongArray
    val rawData = String(sampleNameBuffer.toJavaArray() as ByteArray)

    val sampleNames = offsetsArray.take(numSampleNames).mapIndexedNotNull { index, offset ->
        val end = if (index < offsetsArray.size - 1) {
            offsetsArray[index + 1].toInt()
        } else {
            rawData.length
        }

        if (offset.toInt() >= 0 && offset.toInt() < end) {
            rawData.substring(offset.toInt(), end).trimEnd('\u0000')
        } else {
            null // Skip invalid offsets
        }
    }.toSet() // Use a set to ensure uniqueness

    // Close resources
    query.close()
    array.close()
    context.close()

    return sampleNames
}

fun queryDistinctRefRanges(arrayName: String): Set<String>{
    val context = Context()

    // Open the TileDB array in read mode
    val array = Array(context, arrayName, QueryType.TILEDB_READ)

    // Prepare buffers for SampleName (variable-length strings)
    val refRangeBuffer = NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII) // Adjust size as needed
    val refRangeOffsets = NativeArray(context, 512, Datatype.TILEDB_UINT64) // Adjust size as needed

    // Set up the query
    val query = Query(array, QueryType.TILEDB_READ).apply {
        setLayout(Layout.TILEDB_GLOBAL_ORDER) // Read in global order
        setDataBuffer("RefRange", refRangeBuffer)
        setOffsetsBuffer("RefRange", refRangeOffsets)
    }

    // Submit the query
    query.submit()

    // Extract the number of results
    val numRefRanges = query.resultBufferElements()["RefRange"]?.first?.toInt() ?: 0

    // Extract SampleName data
    val offsetsArray = refRangeOffsets.toJavaArray() as LongArray
    val rawData = String(refRangeBuffer.toJavaArray() as ByteArray)

    val refRanges = offsetsArray.take(numRefRanges).mapIndexedNotNull { index, offset ->
        val end = if (index < offsetsArray.size - 1) {
            offsetsArray[index + 1].toInt()
        } else {
            rawData.length
        }

        if (offset.toInt() >= 0 && offset.toInt() < end) {
            rawData.substring(offset.toInt(), end).trimEnd('\u0000')
        } else {
            null // Skip invalid offsets
        }
    }.toSet() // Use a set to ensure uniqueness

    // Close resources
    query.close()
    array.close()
    context.close()

    return refRanges

}
// THis returns the SampleName and ID1 for a list of RefRanges
fun queryByRefRanges(
    arrayName: String,
    refRangeList: List<String>
): List<Map<String, String>> {
    val context = Context()

    // Open the TileDB array in read mode
    val array = Array(context, arrayName, QueryType.TILEDB_READ)

    println("queryByRefRanges: refRangeList size = ${refRangeList.size}")
    // print the ranges in the list
    refRangeList.forEach { println("queryByRefRanges: refRange = $it") }

    // Prepare buffers for RefRange using prepareVariableBuffer
    val refRangeBuffers = prepareVariableBuffer(context, refRangeList)

    // Prepare buffers for SampleName and ID1
    val sampleNameBuffer = NativeArray(context, 8192, Datatype.TILEDB_STRING_ASCII) // Adjust size
    val sampleNameOffsets = NativeArray(context, 1024, Datatype.TILEDB_UINT64)
    val id1Buffer = NativeArray(context, 8192, Datatype.TILEDB_STRING_ASCII) // Adjust size

    // Create the query
    val query = Query(array, QueryType.TILEDB_READ).apply {
        setLayout(Layout.TILEDB_UNORDERED)

        // Use prepared buffers for RefRange
        setDataBuffer("RefRange", refRangeBuffers.first)
        setOffsetsBuffer("RefRange", refRangeBuffers.second)

        // Set buffers for SampleName and ID1
        setDataBuffer("SampleName", sampleNameBuffer)
        setOffsetsBuffer("SampleName", sampleNameOffsets)
        setDataBuffer("ID1", id1Buffer)
    }

    // Submit the query
    query.submit()

    // Extract RefRange data
    val refRangeOffsetsArrayResult = refRangeBuffers.second.toJavaArray() as LongArray
    val refRangeRawData = String(refRangeBuffers.first.toJavaArray() as ByteArray)
    val refRanges = refRangeOffsetsArrayResult.take(query.resultBufferElements()["RefRange"]?.first?.toInt() ?: 0)
        .mapIndexed { index, offset ->
            val end = if (index < refRangeOffsetsArrayResult.size - 1) refRangeOffsetsArrayResult[index + 1].toInt()
            else refRangeRawData.length
            refRangeRawData.substring(offset.toInt(), end).trimEnd('\u0000')
        }

    // Extract SampleName data
    val sampleNameOffsetsArray = sampleNameOffsets.toJavaArray() as LongArray
    val sampleNameRawData = String(sampleNameBuffer.toJavaArray() as ByteArray)
    val sampleNames = sampleNameOffsetsArray.take(query.resultBufferElements()["SampleName"]?.first?.toInt() ?: 0)
        .mapIndexedNotNull { index, offset ->
            val end = if (index < sampleNameOffsetsArray.size - 1) {
                sampleNameOffsetsArray[index + 1].toInt()
            } else {
                sampleNameRawData.length // last string ends at end of raw data
            }

            // Check bounds and extract the substring if valid
            if (offset.toInt() >= 0 && offset.toInt() < end && end <= sampleNameRawData.length) {
                sampleNameRawData.substring(offset.toInt(), end).trimEnd('\u0000')
            } else {
                null // Skip invalid offsets
            }
        }

    // Extract ID1 data (fixed-length)
    val id1RawData = String(id1Buffer.toJavaArray() as ByteArray)
    val fixedLength = 32 // Fixed length for each ID1 value (MD5 hash)
    val id1s = (0 until (id1RawData.length / fixedLength)).map { index ->
        val start = index * fixedLength
        val end = start + fixedLength
        id1RawData.substring(start, end).trimEnd('\u0000')
    }


    // Combine results into a list of maps
    val results = refRanges.zip(sampleNames.zip(id1s)) { refRange, (sampleName, id1) ->
        mapOf(
            "RefRange" to refRange,
            "SampleName" to sampleName,
            "ID1" to id1
        )
    }
    val filteredResults = results.filter {it["RefRange"] in refRanges}
    println("queryByRefRanges: Size of results = ${results.size} and size of filteredResults = ${filteredResults.size}")

    // Close resources
    query.close()
    array.close()
    context.close()

    return filteredResults
}
