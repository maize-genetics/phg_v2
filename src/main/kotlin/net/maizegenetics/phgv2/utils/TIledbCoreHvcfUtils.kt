package net.maizegenetics.phgv2.utils

import io.tiledb.java.api.*
import io.tiledb.java.api.Array
import io.tiledb.java.api.Constants.TILEDB_VAR_NUM
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * Functions to test java tiledb core api.  This hopefully runs on any platform.
 * The tiledb-java library is platform independent.
 *
 * When you give it an array name, it will create the array in the current directory.
 * And all the data will be written to the array.
 * On bl01, I was in the folder : /workdir/lcj34/phg_v2/tiledbTesting and I told it
 * to create an array named "alt_header_array".  It created a folder named "alt_header_array"
 * and there is data now in it.
 */

private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.utils.TiledbCoreHvcfUtils")
fun parseTiledbAltHeaders(vcfFile: String): List<Map<String, String>> {
    val altData = mutableListOf<Map<String, String>>()

    File(vcfFile).useLines { lines ->
        lines.forEach { line ->
            if (line.startsWith("##ALT=<")) {
                // Subset the line to start after ##ALT=< and end before >
                val lineData = line.substring(7, line.length - 1)
                // Split the lineData on commas
                val lineDataSplit = lineData.split(",")

                // Create a map for storing the selected fields
                val entry = mutableMapOf<String, String>()

                // Extract only the desired fields
                lineDataSplit.forEach { data ->
                    val dataSplit = data.split("=")
                    if (dataSplit.size == 2) { // Ensure valid key=value pairs
                        val key = dataSplit[0].trim()
                        val value = dataSplit[1].removeSurrounding("\"").trim()
                        if (key in listOf("ID", "SampleName", "Regions", "RefChecksum", "RefRange")) {
                            entry[key] = value
                        }
                    }
                }

                // Ensure RefChecksum exists with a default empty value if missing
                // older files stored RefChecksum as the RefRange, and did not use the contig:start-end
                // For this to work, we really need the newer format where RefRange is contig:start-end
                // and RefChecksum is its own field.
                entry.putIfAbsent("RefChecksum", "")

                // Add the entry to the altData list
                altData.add(entry)
            }
        }
    }
    return altData
}

// If there is no path for the arrayName, the array will be created
// in the current directory.
fun createTileDBArray2Dimension(arrayName: String) {
    val context = Context()

    // Define dimensions
    val dimRefRange = Dimension(
        context,
        "RefRange",
        Datatype.TILEDB_STRING_ASCII, // Variable-length string
        null, // No domain for variable-length strings
        null  // No tile extent for variable-length strings
    )

    val dimSampleName = Dimension(
        context,
        "SampleName",
        Datatype.TILEDB_STRING_ASCII, // Variable-length string
        null, // No domain for variable-length strings
        null  // No tile extent for variable-length strings
    )

    // Define domain
    val domain = Domain(context).apply {
        addDimension(dimRefRange)
        addDimension(dimSampleName)
    }

      // Define attributes
    val attrRegions = Attribute(context, "Regions", Datatype.TILEDB_STRING_ASCII).apply {
        setCellValNum(TILEDB_VAR_NUM) // Set as variable-length
    }

    val attrRefChecksum = Attribute(context, "RefChecksum", Datatype.TILEDB_STRING_ASCII).apply {
        setCellValNum(TILEDB_VAR_NUM) // Set as variable-length
    }

    val attrID = Attribute(context, "ID", Datatype.TILEDB_STRING_ASCII).apply {
        setCellValNum(TILEDB_VAR_NUM) // Set as variable-length
    }

    // Define schema
    val schema = ArraySchema(context, ArrayType.TILEDB_SPARSE).apply {
        setDomain(domain)
        addAttribute(attrRegions)
        addAttribute(attrRefChecksum)
        addAttribute(attrID)
    }

    // Create the array
    Array.create(arrayName, schema)

    // Close resources
    schema.close()
    domain.close()
    context.close()

}

// Prepare data and offsets for variable-length fields
// This executes 2 main tasks:
//  1. joins all strings into 1 large string creating a concatenated string buffer
//  2. Calculates the starting index of each string within the concatenated string.
//     This is important for variable-length dimensions or attributes.
fun prepareVariableBuffer(
    context: Context,
    data: List<String>
): Pair<io.tiledb.java.api.NativeArray, io.tiledb.java.api.NativeArray> {
    val concatenated = data.joinToString(separator = "")
    val offsets = data.runningFold(0) { acc, item -> acc + item.length }
        .dropLast(1)
        .map { it.toLong() }
        .toLongArray()

    // WIth variable length attributes, we need both a data buffer and an offsets buffer
    val dataBuffer = io.tiledb.java.api.NativeArray(context, concatenated, Datatype.TILEDB_STRING_ASCII)
    val offsetsBuffer = io.tiledb.java.api.NativeArray(context, offsets, Datatype.TILEDB_UINT64)

    return Pair(dataBuffer, offsetsBuffer)
}

// Function takes an array name and a list of maps, where each map contains
// the fields ID, SampleName, Regions, RefChecksum, and RefRange.
// The function writes the data to the TileDB array.
fun writeAltDataToTileDB(arrayName: String, altData: List<Map<String, String>>) {
    val context = Context()

    // Open the array in write mode
    val array = Array(context, arrayName, QueryType.TILEDB_WRITE)

    // Extract data
    val ids = altData.map { it["ID"].orEmpty() }
    val sampleNames = altData.map { it["SampleName"].orEmpty() }
    val regions = altData.map { it["Regions"].orEmpty() }
    val refChecksums = altData.map { it["RefChecksum"].orEmpty() }
    val refRanges = altData.map { it["RefRange"].orEmpty() }

    // Debug: Print extracted data
//    println("\nLCJ - ids: $ids")
//    println("LCJ - sampleNames: $sampleNames")
//    println("LCJ - regions: $regions")
//    println("LCJ - refChecksums: $refChecksums")
//    println("LCJ - refRanges: $refRanges")

    // Prepare data and offsets for all attributes
    val idBuffers = prepareVariableBuffer(context, ids)
    val sampleNameBuffers = prepareVariableBuffer(context, sampleNames)
    val regionBuffers = prepareVariableBuffer(context, regions)
    val refChecksumBuffers = prepareVariableBuffer(context, refChecksums)
    val refRangeBuffers = prepareVariableBuffer(context, refRanges)

//    println("LCJ - Concatenated SampleName Data: ${sampleNames.joinToString(separator = "")}")
//    println("LCJ - Offsets for SampleName: ${sampleNames.runningFold(0) { acc, s -> acc + s.length }}")
//
//    println("LCJ - Concatenated Regions Data: ${regions.joinToString(separator = "")}")
//    println("LCJ - Offsets for Regions: ${regions.runningFold(0) { acc, s -> acc + s.length }}")
//    println("\nBefore writing to DB:\nConcatenated RefRange Data: ${refRanges.joinToString("")}")
//    println("Before writing to DB:\nRefRange Offsets: ${refRanges.runningFold(0) { acc, s -> acc + s.length }}")

    println("\nwriteALtDataToTIleDB: Sending query now\n")
    // Prepare query
    val query = Query(array, QueryType.TILEDB_WRITE).apply {
        setLayout(Layout.TILEDB_UNORDERED)

        // Set buffers for variable-length attributes
        setDataBuffer("ID", idBuffers.first)
        setOffsetsBuffer("ID", idBuffers.second)

        setDataBuffer("SampleName", sampleNameBuffers.first)
        setOffsetsBuffer("SampleName", sampleNameBuffers.second)

        setDataBuffer("Regions", regionBuffers.first)
        setOffsetsBuffer("Regions", regionBuffers.second)

        setDataBuffer("RefChecksum", refChecksumBuffers.first)
        setOffsetsBuffer("RefChecksum", refChecksumBuffers.second)

        setDataBuffer("RefRange", refRangeBuffers.first)
        setOffsetsBuffer("RefRange", refRangeBuffers.second)
    }

    // Submit and finalize query
    query.submit()
    query.close()
    array.close()
    context.close()
}

// Given a set of IDs, this function return the sampleName and Regions for each ID
// This is just a test of functionality.
fun readSampleNameRegionsForID(arrayName: String, idList: List<String>): List<Map<String, String>> {
    val context = Context()

    // returning a Kotlin  List<Map<String, String>>, and it supports
    // Easy iteration over entries.
    // Direct access to fields by their keys.
    // WOuld we want a data class instead?
    val result = mutableListOf<Map<String, String>>()

    // Open the array in read mode
    val array = Array(context, arrayName, QueryType.TILEDB_READ)

    // Prepare ID buffers
    val idBuffers = prepareVariableBuffer(context, idList)

    // Prepare buffers for SampleName and Regions
    val sampleNameBuffer = io.tiledb.java.api.NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII) // Adjust size as needed
    val sampleNameOffsets = io.tiledb.java.api.NativeArray(context, 512, Datatype.TILEDB_UINT64)       // Adjust size as needed

    val regionsBuffer = io.tiledb.java.api.NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII)    // Adjust size as needed
    val regionsOffsets = io.tiledb.java.api.NativeArray(context, 512, Datatype.TILEDB_UINT64)          // Adjust size as needed

    // Prepare query
    val query = Query(array, QueryType.TILEDB_READ).apply {
        setLayout(Layout.TILEDB_UNORDERED)

        // Set ID buffers
        setDataBuffer("ID", idBuffers.first)
        setOffsetsBuffer("ID", idBuffers.second)

        // Set attribute buffers
        setDataBuffer("SampleName", sampleNameBuffer)
        setOffsetsBuffer("SampleName", sampleNameOffsets)

        setDataBuffer("Regions", regionsBuffer)
        setOffsetsBuffer("Regions", regionsOffsets)
    }

    // Submit query
    query.submit()

    // Process results
    val resultBufferElements = query.resultBufferElements()

    // Extract results for SampleName
    val sampleNameOffsetsArray = sampleNameOffsets.toJavaArray() as LongArray
    val sampleNameRawData = String(sampleNameBuffer.toJavaArray() as ByteArray)
//    println("LCJ - Raw SampleName Data: $sampleNameRawData")
//    println("LCJ - SampleName Offsets: ${sampleNameOffsetsArray.contentToString()}")


    val sampleNames = sampleNameOffsetsArray.mapIndexed { index, offset ->
        val end = if (index < sampleNameOffsetsArray.size - 1) sampleNameOffsetsArray[index + 1].toInt() else sampleNameRawData.length
        sampleNameRawData.substring(offset.toInt(), end).trimEnd('\u0000')
    }.filter { it.isNotEmpty() } // Remove empty strings

    // Extract results for Regions
    val regionsOffsetsArray = regionsOffsets.toJavaArray() as LongArray
    val regionsRawData = String(regionsBuffer.toJavaArray() as ByteArray)
//    println("LCJ - Raw Regions Data: $regionsRawData")
//    println("LCJ - Regions Offsets: ${regionsOffsetsArray.contentToString()}")

    val regions = regionsOffsetsArray.mapIndexed { index, offset ->
        val end = if (index < regionsOffsetsArray.size - 1) regionsOffsetsArray[index + 1].toInt() else regionsRawData.length
        regionsRawData.substring(offset.toInt(), end).trimEnd('\u0000')
    }.filter { it.isNotEmpty() } // Remove empty strings

    // Combine results
    for (idx in sampleNames.indices) {
        result.add(
            mapOf(
                "ID" to idList[idx],
                "SampleName" to sampleNames[idx],
                "Regions" to regions[idx]
            )
        )
    }

    // Close resources
    query.close()
    array.close()
    context.close()

    return result
}

// Given a list of Reference Ranges, query the array for all IDs associated with each Reference Range
fun queryIDsByRefRange(arrayName: String, refRangesToQuery: List<String>): Map<String, List<String>> {
    val context = Context()

    // This one is not a list.  We have a list of refRanges to query, and we want to return a map of refRanges to IDs
    // SO the return is a Map<String, List<String>>.  The key is the refRange, and the value is a list of hapIDs that
    // are associated with each range.
    val result = mutableMapOf<String, List<String>>()

    // Open the array in read mode
    val array = Array(context, arrayName, QueryType.TILEDB_READ)

    // Prepare buffers for reading IDs and RefRanges
    val idBuffer = io.tiledb.java.api.NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII) // Adjust size as needed
    val idOffsetsBuffer = io.tiledb.java.api.NativeArray(context, 512, Datatype.TILEDB_UINT64)

    val refRangeBuffer = io.tiledb.java.api.NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII) // Adjust size as needed
    val refRangeOffsetsBuffer = io.tiledb.java.api.NativeArray(context, 512, Datatype.TILEDB_UINT64)

    // Create a query object
    val query = Query(array, QueryType.TILEDB_READ).apply {
        setLayout(Layout.TILEDB_UNORDERED)

        // Set buffers for ID and RefRange
        setDataBuffer("ID", idBuffer)
        setOffsetsBuffer("ID", idOffsetsBuffer)

        setDataBuffer("RefRange", refRangeBuffer)
        setOffsetsBuffer("RefRange", refRangeOffsetsBuffer)
    }

    // Submit the query
    query.submit()

    // Process results
    val resultBufferElements = query.resultBufferElements()

    // NativeArray from TileDB doesn't have a toByteArray() method and doesn't support
    // conversion using Kotlin's toByteARra90.  So we move to Java to do the conversion
    val idOffsetsArray = idOffsetsBuffer.toJavaArray() as LongArray
    println("LCJ queryIdByRefRange - idOffsetsArray: ${idOffsetsArray.contentToString()}")
    val idRawData = String(idBuffer.toJavaArray() as ByteArray)

    val ids = idOffsetsArray.mapIndexed { index, offset ->
        val end = if (index < idOffsetsArray.size - 1) idOffsetsArray[index + 1].toInt() else idRawData.length
        //idRawData.substring(offset.toInt(), end).trimEnd('\u0000')
        if (offset.toInt() < end) idRawData.substring(offset.toInt(), end).trimEnd('\u0000') else ""
    }.filter { it.isNotEmpty() } // Remove empty strings

    val refRangeOffsetsArray = refRangeOffsetsBuffer.toJavaArray() as LongArray
    val refRangeRawData = String(refRangeBuffer.toJavaArray() as ByteArray)
//    println("LCJ - refRangeRawData: $refRangeRawData")
//    println("LCJ - refRangeOffsetsArray: ${refRangeOffsetsArray.contentToString()}")
    val refRanges = refRangeOffsetsArray.mapIndexed { index, offset ->
        val end = if (index < refRangeOffsetsArray.size - 1) refRangeOffsetsArray[index + 1].toInt() else refRangeRawData.length
        if (offset.toInt() < end) refRangeRawData.substring(offset.toInt(), end).trimEnd('\u0000') else ""
    }.filter { it.isNotEmpty() } // Remove empty strings

    // Filter IDs by RefRange
    refRangesToQuery.forEach { refRangeQuery ->
        val matchingIDs = ids.zip(refRanges)
            .filter { it.second == refRangeQuery }
            .map { it.first }
        result[refRangeQuery] = matchingIDs
    }

    // Close resources
    query.close()
    array.close()
    context.close()

    return result
}






