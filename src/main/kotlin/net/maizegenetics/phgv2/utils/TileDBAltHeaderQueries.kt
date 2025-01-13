package net.maizegenetics.phgv2.utils

import htsjdk.variant.vcf.VCFFileReader
import io.tiledb.java.api.*
import io.tiledb.java.api.Array
import io.tiledb.java.api.Constants.TILEDB_VAR_NUM
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * This class contains methods to query the alternate header of a TileDB VCF array.
 */

/**
 * This function queries the TileDB array for the sampleName and Id associated with refernece ranges
 * in the array name provided.  THis was written assuming the query was for the alt_headers_array -
 * but the general principle should be adapted to both arrays.
 * It uses a stream to continually read results while the query statue is INCOMPLETE.
 * The function returns a list of maps, where each map contains the fields SampleName and ID.
 *
 * Streaming works with smallseq, so probably really isn't streaming.  Needs testing with
 * a larger dataset.
 *
 * This was tested using the alt_header_array, not the hvcf_variants_array
 * THis works when I grab ALL data, then filter at the end.
 *
 * It returns a Map<String, List<Map<String,String>>> where the key is the RefRange and the value is a list of maps
 * THe list of Maps contains the fields SampleName and ID.  There will be an entry on the list for each sampleName,
 * e.g SampleName=LineA, ID=12345; SampleName=LineB, ID=67890, etc.
 *
 */

fun queryWithStreaming_sampleNameIdByRefRange(arrayName: String, refRangeList: List<String>): Map<String,List<Map<String,String>>> {
    val context = Context()

    // This allows for association of all the SampleName/Id pairs with a RefRange
    val resultMap = mutableMapOf<String, List<Map<String,String>>>()

    // Open the array in read mode
    val array = Array(context, arrayName, QueryType.TILEDB_READ)

    // Prepare buffers for RefRange (filtering dimension)
    val refRangeBuffers = TiledbCoreHvcfUtils.prepareVariableBuffer(context, refRangeList)
    println("queryWithStreaming: refRangeList size = ${refRangeList.size},  refRangeBuffers.first size = ${refRangeBuffers.first.size}, second size = ${refRangeBuffers.second.size}")
    // print the refRangeBuffers.second as a long array
    val refRangeOffsetsArray = refRangeBuffers.second.toJavaArray() as LongArray
    println("queryWithStreaming: refRangeOffsetsArray = ${refRangeOffsetsArray.contentToString()}")

    // Set reasonable initial buffer sizes
    val sampleNameBufferSize = 8192
    val sampleNameOffsetsSize = 1024
    val idBufferSize = 8192
    val idOffsetsSize = 1024

    // Allocate initial buffers
    val sampleNameBuffer = io.tiledb.java.api.NativeArray(context, sampleNameBufferSize, Datatype.TILEDB_STRING_ASCII)
    val sampleNameOffsets = io.tiledb.java.api.NativeArray(context, sampleNameOffsetsSize, Datatype.TILEDB_UINT64)
    val idBuffer = io.tiledb.java.api.NativeArray(context, idBufferSize, Datatype.TILEDB_STRING_ASCII)
    val idOffsets = io.tiledb.java.api.NativeArray(context, idOffsetsSize, Datatype.TILEDB_UINT64)

    // Create the query
    val query = Query(array, QueryType.TILEDB_READ).apply {
        setLayout(Layout.TILEDB_UNORDERED) // GLOBAL_ORDER didn't get all the data, and associated it incorrectly

        // Set RefRange Buffers  - these are the buffers from prepareVariableBuffer(refRangeList)
        setDataBuffer("RefRange", refRangeBuffers.first)
        setOffsetsBuffer("RefRange", refRangeBuffers.second)

        // Set SampleName buffers
        setDataBuffer("SampleName", sampleNameBuffer)
        setOffsetsBuffer("SampleName", sampleNameOffsets)

        // Set ID buffers
        setDataBuffer("ID", idBuffer)
        setOffsetsBuffer("ID", idOffsets)
    }

    // Add ranges for RefRange dimension. This will allow tiledb to do the "slicing" for us
    // This doesn't work due to conda issues.  It seems we need a core tiledb of 12.2.0 or later
    // I tried adding and creating a new environment but could not get the env to build on the iMac.
    // SO for now, we will just get all the data and filter it at the end.
//    refRangeList.forEach { refRange ->
//        query.addRange(0, refRange, refRange, null) // Constrain to specific RefRange values
//    }

    // Iterate until query is complete
    do {
        query.submit()

        // Use resultBufferElements to determine the number of valid entries
        val numSampleNames = query.resultBufferElements()["SampleName"]?.first?.toInt() ?: 0
        val numIDs = query.resultBufferElements()["ID"]?.first?.toInt() ?: 0
        val numRefRanges = query.resultBufferElements()["RefRange"]?.first?.toInt() ?: 0
        println("queryWithStreaming the query.submit results show : num SampleNames=$numSampleNames, num IDs=$numIDs, num RefRanges=$numRefRanges")


        // Extract RefRange data
        val refRangeOffsetsArrayResult = refRangeBuffers.second.toJavaArray() as LongArray
        val refRangeRawData = String(refRangeBuffers.first.toJavaArray() as ByteArray)
        val refRanges = refRangeOffsetsArrayResult.take(query.resultBufferElements()["RefRange"]?.first?.toInt() ?: 0)
            .mapIndexedNotNull { index, offset ->
                val end = if (index < refRangeOffsetsArrayResult.size - 1) refRangeOffsetsArrayResult[index + 1].toInt()
                else refRangeRawData.length
                refRangeRawData.substring(offset.toInt(), end).trimEnd('\u0000')
                if (offset.toInt() < end) {
                    refRangeRawData.substring(offset.toInt(), end).trimEnd('\u0000')
                } else {
                    null
                }
            }

        // Extract SampleName data
        val sampleNameOffsetsArray = sampleNameOffsets.toJavaArray() as LongArray
        val sampleNameRawData = String(sampleNameBuffer.toJavaArray() as ByteArray)
        val sampleNames = sampleNameOffsetsArray.take(numSampleNames).mapIndexedNotNull { index, offset ->
            val end = if (index < sampleNameOffsetsArray.size - 1) sampleNameOffsetsArray[index + 1].toInt() else sampleNameRawData.length
            if (offset.toInt() >=0 && offset.toInt() < end && end <= sampleNameRawData.length) {
                sampleNameRawData.substring(offset.toInt(), end).trimEnd('\u0000')
            } else null
        }

        // Extract ID data
        val idOffsetsArray = idOffsets.toJavaArray() as LongArray
        val idRawData = String(idBuffer.toJavaArray() as ByteArray)
        val ids = idOffsetsArray.take(numIDs).mapIndexedNotNull { index, offset ->
            val end = if (index < idOffsetsArray.size - 1) idOffsetsArray[index + 1].toInt() else idRawData.length
            if (offset.toInt() >=0 && offset.toInt() < end && end <= idRawData.length) {
                idRawData.substring(offset.toInt(), end).trimEnd('\u0000')
            } else null
        }

        // Associate the data, store to the resultMap
        refRanges.forEachIndexed { index, refRange ->
            val sampleName = sampleNames.getOrNull(index) ?: ""
            val id = ids.getOrNull(index) ?: ""
            resultMap[refRange] = resultMap.getOrDefault(refRange,emptyList()) + mapOf("SampleName" to sampleName, "ID" to id)
        }

    } while (query.queryStatus == QueryStatus.TILEDB_INCOMPLETE)

    // Close resources
    query.close()
    array.close()
    context.close()

    // Above got ALL The ref ranges.  TileDB only slices if you add ranges to the query.
    // (See commented out code above for query.addRange() )
    // We need a newer tiledb version for that, and I couldn't get it to build on the iMac.
    // Since I have limited time before leaving, I am filtering after retrieving the data.
    // See the commented out code above. (LCJ)

    //Filter the results to just what we requested
    val filteredResults = resultMap.filterKeys { it in refRangeList }
    println("queryWithStreaming: number of  resultMap entries = ${resultMap.keys.size}, num filtered results = ${filteredResults.keys.size}")
    //println("\nHere are all the values from the resultMap before filtering:")
    resultMap.forEach { (key, value) -> println("Key: $key, Value: $value") }

    return filteredResults
}

/**
 * Given a list of refRanges, query the array for all IDs associated with each reference range
 *
 * This function uses streaming ro continually read tiledb results until the query status is TILEDB_COMPLETE
 * (ie, no longer TILEDB_INCOMPLETE).  It returns a map of refRanges to a list of IDs.
 */
fun queryIDsByRefRange(arrayName: String, refRangesToQuery: List<String>): Map<String, List<String>> {
    val context = Context()

    // This one is not a list.  We have a list of refRanges to query, and we want to return a map of refRanges to IDs
    // SO the return is a Map<String, List<String>>.  The key is the refRange, and the value is a list of hapIDs that
    // are associated with each range.
    val result = mutableMapOf<String, List<String>>()

    // Open the array in read mode
    val array = Array(context, arrayName, QueryType.TILEDB_READ)

    // Prepare buffers for reading IDs and RefRanges
    val idBuffer = NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII) // Adjust size as needed
    val idOffsetsBuffer = NativeArray(context, 512, Datatype.TILEDB_UINT64)

    val refRangeBuffer = NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII) // Adjust size as needed
    val refRangeOffsetsBuffer = NativeArray(context, 512, Datatype.TILEDB_UINT64)

    // Create a query object
    val query = Query(array, QueryType.TILEDB_READ).apply {
        setLayout(Layout.TILEDB_GLOBAL_ORDER)

        // Set buffers for ID and RefRange
        setDataBuffer("ID", idBuffer)
        setOffsetsBuffer("ID", idOffsetsBuffer)

        setDataBuffer("RefRange", refRangeBuffer)
        setOffsetsBuffer("RefRange", refRangeOffsetsBuffer)
    }

    do {
        query.submit()

        // Extract RefRange data
        val refRangeOffsetsArray = refRangeOffsetsBuffer.toJavaArray() as LongArray
        val refRangeRawData = String(refRangeBuffer.toJavaArray() as ByteArray)
        val refRanges = refRangeOffsetsArray.mapIndexed { index, offset ->
            val end = if (index < refRangeOffsetsArray.size - 1) refRangeOffsetsArray[index + 1].toInt() else refRangeRawData.length
            if (offset.toInt() < end) refRangeRawData.substring(offset.toInt(), end).trimEnd('\u0000') else ""
        }.filter { it.isNotEmpty() } // Remove empty strings

        // Extract ID data
        val idOffsetsArray = idOffsetsBuffer.toJavaArray() as LongArray
        val idRawData = String(idBuffer.toJavaArray() as ByteArray)
        val ids = idOffsetsArray.mapIndexed { index, offset ->
            val end = if (index < idOffsetsArray.size - 1) idOffsetsArray[index + 1].toInt() else idRawData.length
            if (offset.toInt() < end) idRawData.substring(offset.toInt(), end).trimEnd('\u0000') else ""
        }.filter { it.isNotEmpty() } // Remove empty strings

        // Filter IDs by RefRange
        refRangesToQuery.forEach { refRangeQuery ->
            val matchingIDs = ids.zip(refRanges)
                .filter { it.second == refRangeQuery }
                .map { it.first }
            result[refRangeQuery] = matchingIDs
        }

    } while (query.queryStatus == QueryStatus.TILEDB_INCOMPLETE)

    // Close resources
    query.close()
    array.close()
    context.close()

    return result
}