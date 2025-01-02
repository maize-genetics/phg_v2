package net.maizegenetics.phgv2.utils

import htsjdk.variant.vcf.VCFFileReader
import io.tiledb.java.api.*
import io.tiledb.java.api.Array
import io.tiledb.java.api.Constants.TILEDB_VAR_NUM
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * Functions to test java tiledb core api.  This hopefully runs on any platform.
 * The tiledb-java library is platform independent.
 *
 * Creating the array: When given  an array name, tiledb will create the array in the current directory
 * if there has been no path included.
 * And all the data will be written to the array.
 * On bl01, I was in the folder : /workdir/lcj34/phg_v2/tiledbTesting and I told it
 * to create an array named "alt_header_array".  It created a folder named "alt_header_array"
 * and there is data now in it.
 *
 * TODO:  must parse the body to get gt and sample names, this has tobe reworked to include body parseing/GT values
 * and to create the SamleGamete - which will still be a string, but with _0 _1 etc appended to the sample name.
 */

private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.utils.TiledbCoreHvcfUtils")


fun parseTiledbAltHeaders(reader: VCFFileReader): List<Map<String,String>> {
    val altData = mutableListOf<Map<String, String>>()

    val altHeaders = parseALTHeader(reader.header)
    // Need to put the data into the map above
    altHeaders.forEach{ header ->
        val entry = mutableMapOf<String, String>()
        entry["ID"] = header.value.id
        entry["SampleName"] = header.value.sampleName()
        entry["Regions"] = header.value.regions.map { it.first.toString() + "-" + it.second.position.toString() }.joinToString(",")
        entry["RefRange"] = header.value.refRange
        entry["RefChecksum"] = header.value.refChecksum
        altData.add(entry)
    }
    return altData
}
// Rather than procsesing the files with 2 methods, which would
// require reading the file twice, we'll read the file once and
// process the header and body in one pass.  Below might need to be
// re-written using htsjdk to read the file and parse the header.
fun parseTiledbAltHeadersORIG(vcfFile: String): List<Map<String, String>> {
    val altData = mutableListOf<Map<String, String>>()

    // Regex pattern to match key-value pairs, accounting for quoted values
    val regex = Regex("""(\w+)=("(.*?)"|[^,]+)""")

    File(vcfFile).useLines { lines ->
        lines.forEach { line ->
            if (line.startsWith("##ALT=<")) {
                // Subset the line to start after ##ALT=< and end before >
                val lineData = line.substring(7, line.length - 1)

                // Create a map for storing the selected fields
                val entry = mutableMapOf<String, String>()

                // Extract key-value pairs using the regex
                regex.findAll(lineData).forEach { matchResult ->
                    val key = matchResult.groupValues[1].trim()
                    val value = matchResult.groupValues[2].removeSurrounding("\"").trim()
                    if (key in listOf("ID", "SampleName", "Regions", "RefChecksum", "RefRange")) {
                        entry[key] = value
                    }
                }

                // Ensure RefChecksum exists with a default empty value if missing
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
fun createTileDBCoreArrays(dbPath: String) {
    val context = Context()

    // Create the alt header array
    var arrayNameAlt = dbPath + "/alt_header_array"
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

    // Cannot set fixed-length string for dimensions - they must be variable length.  ONly attributes
    // canbe fixed length.
    val dimID = Dimension(context, "ID",Datatype.TILEDB_STRING_ASCII, null, null)

    // CUrrently ordering domains as RefRange, SampleName and ID in that order.
    // ID should definitlye be last, but is a toss up between RefRange and SampleName
    // for which should be first.
    // Define domain
    val domain = Domain(context).apply {
        addDimension(dimRefRange)
        addDimension(dimSampleName)
        addDimension(dimID)
    }

      // Define attributes
    val attrRegions = Attribute(context, "Regions", Datatype.TILEDB_STRING_ASCII).apply {
        setCellValNum(TILEDB_VAR_NUM) // Set as variable-length
    }

    val attrRefChecksum = Attribute(context, "RefChecksum", Datatype.TILEDB_STRING_ASCII).apply {
        setCellValNum(TILEDB_VAR_NUM) // Set as variable-length
    }

    // ID is now a dimension, not an attribute
//    val attrID = Attribute(context, "ID", Datatype.TILEDB_STRING_ASCII).apply {
//        setCellValNum(TILEDB_VAR_NUM) // Set as variable-length - which is wrong, these are set at 32 bytes for md5 hash
//    }

    // Define schema
    val schema = ArraySchema(context, ArrayType.TILEDB_SPARSE).apply {
        setDomain(domain)
        addAttribute(attrRegions)
        addAttribute(attrRefChecksum)
        //addAttribute(attrID)  // attrID is now a dimension, not an attribute
    }

    // Create the array
    Array.create(arrayNameAlt, schema)

    // Close resources
    schema.close()
    domain.close()
    //context.close()  WIll use for 2nd array - don't close it yet.

    // Create the variants array
    // Dimensions are sampleName and refRange:  these are both strings
    // Attributes are ID1 and ID2 - these are currently string, will someday be uint64
    val arrayNameVariants = dbPath + "/hvcf_variants_array"

    val domain2 = Domain(context).apply {
        addDimension(Dimension(context, "RefRange", Datatype.TILEDB_STRING_ASCII, null, null))
        addDimension(Dimension(context, "SampleName", Datatype.TILEDB_STRING_ASCII, null, null))
    }

    // Setting fixed size based on MD5 hash.  If IDs change to uint64 this will have to change as well
    val attrID1 = Attribute(context, "ID1", Datatype.TILEDB_STRING_ASCII).apply {
        setCellValNum(32) // Define as fixed-length with 32 characters
    }
    val attrID2 = Attribute(context, "ID2", Datatype.TILEDB_STRING_ASCII).apply {
        setCellValNum(32) // Define as fixed-length with 32 characters
    }

    val schema2 = ArraySchema(context, ArrayType.TILEDB_SPARSE).apply {
        setDomain(domain2) // same domain as above
        addAttribute(attrID1)
        addAttribute(attrID2)
    }

    // Create the array
    Array.create(arrayNameVariants, schema2)

    // Close resources
    schema2.close()
    domain2.close()
    context.close()

}

// Prepare data and offsets for variable-length fields
// This executes 2 main tasks:
//  1. joins all strings into 1 large string creating a concatenated string buffer
//  2. Calculates the starting index of each string within the concatenated string.
//     This is important for variable-length dimensions or attributes.
// This function should be used when loading data to a tiledb, or for creating
// the buffers for a specific list of values for a dimension (not attribute) that is being queried.
// For example:  if you want to query all the samplenames and IDs (hapids) for a select set
// of RefRanges, you would use this function to create the buffers for the RefRanges.
fun prepareVariableBuffer(
    context: Context,
    data: List<String>
): Pair<io.tiledb.java.api.NativeArray, io.tiledb.java.api.NativeArray> {
    val concatenated = data.joinToString(separator = "")

    // Create a running total of the lengths of the strings in the data list, dropping
    // the last element as that is not needed, then convert to a long array.
    val offsets = data.runningFold(0) { acc, item -> acc + item.length }
        .dropLast(1)
        .map { it.toLong() } // tiledb wants this data as TILEDB_UINT64, which is Kotlin long
        .toLongArray()

    // WIth variable length attributes, we need both a data buffer and an offsets buffer
    val dataBuffer = io.tiledb.java.api.NativeArray(context, concatenated, Datatype.TILEDB_STRING_ASCII)
    val offsetsBuffer = io.tiledb.java.api.NativeArray(context, offsets, Datatype.TILEDB_UINT64)

    return Pair(dataBuffer, offsetsBuffer)
}


// Function takes an array name and a list of maps, where each map contains
// the fields ID, SampleName, Regions, RefChecksum, and RefRange.
// The function writes the data to the TileDB array.
// This does not deal with the contexts/variants of the hvcf file, only the alt headers.
fun writeAltDataToTileDB(arrayName: String, altData: List<Map<String, String>>) {
    val context = Context()

    // Open the array in write mode
    println("writeAltDataToTIleDB - arrayName: $arrayName")
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

/**
 * This function writes hvcf variants data to the hvcf_variants_array
 */
fun writeVariantsDataToTileDB(variantArrayName:String, combinedHvcfVariantData:List<Map<String, String>>) {
    val context = Context()

    // Open the array in write mode
    val array = Array(context, variantArrayName, QueryType.TILEDB_WRITE)

    // Extract data
    val samples = combinedHvcfVariantData.map { it["SampleName"].orEmpty() }
    val refRanges = combinedHvcfVariantData.map { it["RefRange"].orEmpty() }

    // Prepare data and offsets for all attributes
    val sampleBuffers = prepareVariableBuffer(context, samples) // SampleNames are variable length
    val refRangeBuffers = prepareVariableBuffer(context, refRanges) // RefRanges are variable length

    // This value must change if we go to uint64 for the IDs - 32 is when using MD5 hash
    val fixedLength = 32 // fixed length buffer while we use the MD5 hash, which is 32 bytes
    val id1s = combinedHvcfVariantData.map { it["ID1"].orEmpty().padEnd(fixedLength) }
    val id1Buffer = NativeArray(context, id1s.joinToString(""), Datatype.TILEDB_STRING_ASCII)

    val id2s = combinedHvcfVariantData.map { it["ID2"].orEmpty().padEnd(fixedLength) }
    val id2Buffer = NativeArray(context, id2s.joinToString(""), Datatype.TILEDB_STRING_ASCII)

    // Prepare query
    val query = Query(array, QueryType.TILEDB_WRITE).apply {
        setLayout(Layout.TILEDB_UNORDERED)

        // Set buffers for variable-length attributes
        setDataBuffer("SampleName", sampleBuffers.first)
        setOffsetsBuffer("SampleName", sampleBuffers.second)

        setDataBuffer("RefRange", refRangeBuffers.first)
        setOffsetsBuffer("RefRange", refRangeBuffers.second)

        // Set buffers for fixed length attributes
        setDataBuffer("ID1", id1Buffer)
        setDataBuffer("ID2", id2Buffer)
    }

    // Submit and finalize query
    query.submit()
    query.close()
    array.close()
    context.close()
}

/**
 * This function queries the TileDB array for the sampleName and Id associated with refernece ranges
 * in the array name provided.  THis was written assuming the query wasfor the alt_headers_array -
 * but the general principle shoudl be adapted to both arrays.
 * It uses a dynamic buffer size to handle incomplete queries, continuing to query until the query is complete.
 * The function returns a list of maps, where each map contains the fields SampleName and ID.
 *
 * THIS HAS NOT BEEN TESTED - it is here because we'll probably need something like it when testing at scale
 * as we won't know the sizeof the data we are querying.
 */
fun queryWithDynamicBuffers(arrayName: String, refRangeList: List<String>): List<Map<String, String>> {
    val context = Context()
    val result = mutableListOf<Map<String, String>>()

    // Open the array in read mode
    val array = Array(context, arrayName, QueryType.TILEDB_READ)

    // Prepare buffers for RefRange (filtering dimension)
    val refRangeBuffers = prepareVariableBuffer(context, refRangeList)

    // Initial buffer sizes
    var sampleNameBufferSize = 4096
    var sampleNameOffsetsSize = 512
    var idBufferSize = 4096
    var idOffsetsSize = 512

    // Initialize query buffers
    var sampleNameBuffer = io.tiledb.java.api.NativeArray(context, sampleNameBufferSize, Datatype.TILEDB_STRING_ASCII)
    var sampleNameOffsets = io.tiledb.java.api.NativeArray(context, sampleNameOffsetsSize, Datatype.TILEDB_UINT64)

    var idBuffer = io.tiledb.java.api.NativeArray(context, idBufferSize, Datatype.TILEDB_STRING_ASCII)
    var idOffsets = io.tiledb.java.api.NativeArray(context, idOffsetsSize, Datatype.TILEDB_UINT64)

    val query = Query(array, QueryType.TILEDB_READ).apply {
        setLayout(Layout.TILEDB_UNORDERED)

        // Set RefRange buffers
        setDataBuffer("RefRange", refRangeBuffers.first)
        setOffsetsBuffer("RefRange", refRangeBuffers.second)

        // Set SampleName buffers
        setDataBuffer("SampleName", sampleNameBuffer)
        setOffsetsBuffer("SampleName", sampleNameOffsets)

        // Set ID buffers
        setDataBuffer("ID", idBuffer)
        setOffsetsBuffer("ID", idOffsets)
    }

    // Iterate until query is complete
    do {
        query.submit()

        if (query.queryStatus == QueryStatus.TILEDB_INCOMPLETE) {
            println("Query incomplete. Expanding buffer sizes.")

            // Double buffer sizes
            sampleNameBufferSize *= 2
            sampleNameOffsetsSize *= 2
            idBufferSize *= 2
            idOffsetsSize *= 2

            // Reallocate buffers
            sampleNameBuffer = io.tiledb.java.api.NativeArray(context, sampleNameBufferSize, Datatype.TILEDB_STRING_ASCII)
            sampleNameOffsets = io.tiledb.java.api.NativeArray(context, sampleNameOffsetsSize, Datatype.TILEDB_UINT64)
            idBuffer = io.tiledb.java.api.NativeArray(context, idBufferSize, Datatype.TILEDB_STRING_ASCII)
            idOffsets = io.tiledb.java.api.NativeArray(context, idOffsetsSize, Datatype.TILEDB_UINT64)

            // Reset buffers in the query
            query.setDataBuffer("SampleName", sampleNameBuffer)
            query.setOffsetsBuffer("SampleName", sampleNameOffsets)
            query.setDataBuffer("ID", idBuffer)
            query.setOffsetsBuffer("ID", idOffsets)
        }

    } while (query.queryStatus == QueryStatus.TILEDB_INCOMPLETE)

    // Use resultBufferElements to determine the number of valid entries
    val numSampleNames = query.resultBufferElements()["SampleName"]?.first?.toInt() ?: 0
    val numIDs = query.resultBufferElements()["ID"]?.first?.toInt() ?: 0

    // Extract SampleName data
    val sampleNameOffsetsArray = sampleNameOffsets.toJavaArray() as LongArray
    val sampleNameRawData = String(sampleNameBuffer.toJavaArray() as ByteArray)
    val sampleNames = sampleNameOffsetsArray.take(numSampleNames).mapIndexedNotNull { index, offset ->
        val end = if (index < sampleNameOffsetsArray.size - 1) sampleNameOffsetsArray[index + 1].toInt() else sampleNameRawData.length
        if (offset.toInt() < end) sampleNameRawData.substring(offset.toInt(), end).trimEnd('\u0000') else null
    }

    // Extract ID data
    val idOffsetsArray = idOffsets.toJavaArray() as LongArray
    val idRawData = String(idBuffer.toJavaArray() as ByteArray)
    val ids = idOffsetsArray.take(numIDs).mapIndexedNotNull { index, offset ->
        val end = if (index < idOffsetsArray.size - 1) idOffsetsArray[index + 1].toInt() else idRawData.length
        if (offset.toInt() < end) idRawData.substring(offset.toInt(), end).trimEnd('\u0000') else null
    }

    // Combine results
    sampleNames.zip(ids).forEach { (sampleName, id) ->
        result.add(mapOf("SampleName" to sampleName, "ID" to id))
    }

    // Close resources
    query.close()
    array.close()
    context.close()

    return result
}



/**
 * Given a list of reference ranges, return the sampleNames and Ids (hapids) associated with each range.
 *
 */
fun querySampleNamesAndIDsByRefRange(arrayName: String, refRangeList: List<String>): Map<String, List<Map<String, String>>> {
    val context = Context()
    val result = mutableMapOf<String, List<Map<String, String>>>()

    // Open the array in read mode
    val array = Array(context, arrayName, QueryType.TILEDB_READ)

    // Prepare buffers for RefRange
    val refRangeBuffers = prepareVariableBuffer(context, refRangeList)

    // Prepare buffers for SampleName and ID
    val sampleNameBuffer = io.tiledb.java.api.NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII)
    val sampleNameOffsets = io.tiledb.java.api.NativeArray(context, 512, Datatype.TILEDB_UINT64)

    val idBuffer = io.tiledb.java.api.NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII)
    val idOffsets = io.tiledb.java.api.NativeArray(context, 512, Datatype.TILEDB_UINT64)

    // Create and configure the query
    val query = Query(array, QueryType.TILEDB_READ).apply {
        setLayout(Layout.TILEDB_UNORDERED)

        // Set buffers for RefRange (dimension filter)
        setDataBuffer("RefRange", refRangeBuffers.first)
        setOffsetsBuffer("RefRange", refRangeBuffers.second)

        // Set buffers for SampleName (dimension retrieval)
        setDataBuffer("SampleName", sampleNameBuffer)
        setOffsetsBuffer("SampleName", sampleNameOffsets)

        // Set buffers for ID (attribute retrieval)
        setDataBuffer("ID", idBuffer)
        setOffsetsBuffer("ID", idOffsets)
    }

    // Submit the query
    query.submit()

    // Use resultBufferElements to determine the number of valid entries
    val numSampleNames = query.resultBufferElements()["SampleName"]?.first?.toInt() ?: 0
    val numIDs = query.resultBufferElements()["ID"]?.first?.toInt() ?: 0
    val numRefRanges = query.resultBufferElements()["RefRange"]?.first?.toInt() ?: 0
    println("LCJ - Num Results: SampleNames=$numSampleNames, IDs=$numIDs, RefRanges=$numRefRanges")

    // Extract SampleName data
    val sampleNameOffsetsArray = sampleNameOffsets.toJavaArray() as LongArray
    val sampleNameRawData = String(sampleNameBuffer.toJavaArray() as ByteArray)
    println("Offsets: ${sampleNameOffsetsArray.contentToString()}")
    println("Raw Data: $sampleNameRawData")
    println("Num SampleNames: $numSampleNames")

    val sampleNames = sampleNameOffsetsArray.take(numSampleNames).mapIndexedNotNull { index, offset ->
        val end = if (index < numSampleNames - 1) {
            sampleNameOffsetsArray[index + 1].toInt() // Next offset as end
        } else {
            sampleNameRawData.length // Last string ends at the end of raw data
        }

        if (offset.toInt() >= 0 && offset.toInt() < end && end <= sampleNameRawData.length) {
            sampleNameRawData.substring(offset.toInt(), end).trimEnd('\u0000')
        } else {
            null // Skip invalid offsets
        }
    }
    println("Parsed SampleNames: $sampleNames")

    // Extract ID data
    val idOffsetsArray = idOffsets.toJavaArray() as LongArray
    val idRawData = String(idBuffer.toJavaArray() as ByteArray)

    val ids = idOffsetsArray.take(numIDs).mapIndexedNotNull { index, offset ->
        val end = if (index < numIDs - 1) {
            idOffsetsArray[index + 1].toInt() // Next offset as end
        } else {
            idRawData.length // Last string ends at the end of raw data
        }

        if (offset.toInt() >= 0 && offset.toInt() < end && end <= idRawData.length) {
            idRawData.substring(offset.toInt(), end).trimEnd('\u0000')
        } else {
            null // Skip invalid offsets
        }
    }
    println("Parsed IDs: $ids")

    // Combine results grouped by RefRange
    val refRangeOffsetsArray = refRangeBuffers.second.toJavaArray() as LongArray
    val refRangeRawData = String(refRangeBuffers.first.toJavaArray() as ByteArray)
    val refRanges = refRangeOffsetsArray.take(numRefRanges).mapIndexedNotNull { index, offset ->
        val end = if (index < refRangeOffsetsArray.size - 1) refRangeOffsetsArray[index + 1].toInt() else refRangeRawData.length
        if (offset.toInt() < end) refRangeRawData.substring(offset.toInt(), end).trimEnd('\u0000') else null
    }

    refRanges.forEachIndexed { index, refRange ->
        val sampleName = sampleNames.getOrNull(index) ?: ""
        val id = ids.getOrNull(index) ?: ""
        result[refRange] = result.getOrDefault(refRange, emptyList()) + mapOf(
            "SampleName" to sampleName,
            "ID" to id
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
    //println("LCJ queryIdByRefRange - idOffsetsArray: ${idOffsetsArray.contentToString()}")
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






