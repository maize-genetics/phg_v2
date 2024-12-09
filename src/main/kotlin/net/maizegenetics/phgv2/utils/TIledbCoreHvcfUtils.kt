package net.maizegenetics.phgv2.utils

import io.tiledb.java.api.*
import io.tiledb.java.api.Array
import io.tiledb.java.api.Constants.TILEDB_VAR_NUM
import java.io.File

/**
 * Functions to test java tiledb core api.  This hopefully runs on any platform.
 * The tiledb-java library is platform independent.
 * When you give it an array name, it will create the array in the current directory.
 * And all the data will be written to the array.
 * On bl01, I was in the folder : /workdir/lcj34/phg_v2/tiledbTesting and I told it
 * to create an array named "alt_header_array".  It created a folder named "alt_header_array"
 * and there is data now in it.
 * TODO: test the array can be created. Then write code to test that I can pull data from it
 */
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
                entry.putIfAbsent("RefChecksum", "")

                // Add the entry to the altData list
                altData.add(entry)
            }
        }
    }
    return altData
}

fun createTileDBArray(arrayName: String) {
    val context = Context()

    // Define dimensions (ID as a string dimension)
    val dimID = Dimension(
        context,
        "ID",
        Datatype.TILEDB_STRING_ASCII, // String data type
        null,                         // No domain for string dimensions
        null                          // No tile extent for string dimensions
    )

    // Define domain
    val domain = Domain(context).addDimension(dimID)

    // Define attributes (all variable-length)
    val attrSampleName = Attribute(context, "SampleName", Datatype.TILEDB_STRING_ASCII).apply {
        setCellValNum(TILEDB_VAR_NUM) // Set as variable-length
    }
    val attrRegions = Attribute(context, "Regions", Datatype.TILEDB_STRING_ASCII).apply {
        setCellValNum(TILEDB_VAR_NUM) // Set as variable-length
    }
    val attrRefChecksum = Attribute(context, "RefChecksum", Datatype.TILEDB_STRING_ASCII).apply {
        setCellValNum(TILEDB_VAR_NUM) // Set as variable-length
    }
    val attrRefRange = Attribute(context, "RefRange", Datatype.TILEDB_STRING_ASCII).apply {
        setCellValNum(TILEDB_VAR_NUM) // Set as variable-length
    }

    // Create schema
    val schema = ArraySchema(context, ArrayType.TILEDB_SPARSE).apply {
        setDomain(domain)
        addAttribute(attrSampleName)
        addAttribute(attrRegions)
        addAttribute(attrRefChecksum)
        addAttribute(attrRefRange)
    }

    // Use ArraySchema to create the array
    Array.create(arrayName, schema)

    // Close resources
    schema.close()
    domain.close()
    context.close()
}

fun createTileDBArrayORIG(arrayName: String) {
    val context = Context()

    // Define dimensions (ID as a string dimension)
    val dimID = Dimension(
        context,
        "ID",
        Datatype.TILEDB_STRING_ASCII, // String data type
        null,                         // No domain for string dimensions
        null                          // No tile extent for string dimensions
    )

    // Define domain
    val domain = Domain(context).addDimension(dimID)

    // Define attributes
    val attrSampleName = Attribute(context, "SampleName", Datatype.TILEDB_STRING_ASCII)
    val attrRegions = Attribute(context, "Regions", Datatype.TILEDB_STRING_ASCII)
    val attrRefChecksum = Attribute(context, "RefChecksum", Datatype.TILEDB_STRING_ASCII)
    val attrRefRange = Attribute(context, "RefRange", Datatype.TILEDB_STRING_ASCII)

    // Create schema - note:  DENSE arrays do not support dimension datatype 'STRING_ASCII'
    // What does that mean for us?
    val schema = ArraySchema(context, ArrayType.TILEDB_SPARSE).apply {
        setDomain(domain)
        addAttribute(attrSampleName)
        addAttribute(attrRegions)
        addAttribute(attrRefChecksum)
        addAttribute(attrRefRange)
    }

    // Use ArraySchema to create the array
    //schema.create(arrayName)
    Array.create(arrayName, schema)

    // Close resources
    schema.close()
    domain.close()
    context.close()
}


// Prepare data and offsets for variable-length fields
fun prepareVariableBuffer(
    context: Context,
    data: List<String>
): Pair<io.tiledb.java.api.NativeArray, io.tiledb.java.api.NativeArray> {
    val concatenated = data.joinToString(separator = "")
    val offsets = data.runningFold(0) { acc, item -> acc + item.length }
        .dropLast(1)
        .map { it.toLong() }
        .toLongArray()

    val dataBuffer = io.tiledb.java.api.NativeArray(context, concatenated, Datatype.TILEDB_STRING_ASCII)
    val offsetsBuffer = io.tiledb.java.api.NativeArray(context, offsets, Datatype.TILEDB_UINT64)

    return Pair(dataBuffer, offsetsBuffer)
}

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
    println("\nLCJ - ids: $ids")
    println("LCJ - sampleNames: $sampleNames")
    println("LCJ - regions: $regions")
    println("LCJ - refChecksums: $refChecksums")
    println("LCJ - refRanges: $refRanges")

    // Prepare data and offsets for all attributes
    val idBuffers = prepareVariableBuffer(context, ids)
    val sampleNameBuffers = prepareVariableBuffer(context, sampleNames)
    val regionBuffers = prepareVariableBuffer(context, regions)
    val refChecksumBuffers = prepareVariableBuffer(context, refChecksums)
    val refRangeBuffers = prepareVariableBuffer(context, refRanges)

    println("LCJ - Concatenated SampleName Data: ${sampleNames.joinToString(separator = "")}")
    println("LCJ - Offsets for SampleName: ${sampleNames.runningFold(0) { acc, s -> acc + s.length }}")

    println("LCJ - Concatenated Regions Data: ${regions.joinToString(separator = "")}")
    println("LCJ - Offsets for Regions: ${regions.runningFold(0) { acc, s -> acc + s.length }}")


    println("\nwriteALtDataToTIleDBSending query now\n")
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

fun readDataFromTileDB(arrayName: String, idList: List<String>): List<Map<String, String>> {
    val context = Context()
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
//    println("LCJ - Buffer Elements: $resultBufferElements")

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

//    println("LCJ - Processed SampleNames: $sampleNames")
//    println("LCJ - Processed Regions: $regions")

    // Combine results
    for (i in sampleNames.indices) {
        result.add(
            mapOf(
                "ID" to idList[i],
                "SampleName" to sampleNames[i],
                "Regions" to regions[i]
            )
        )
    }

    // Debug final results
    println("LCJ - Final Results: $result")

    // Close resources
    query.close()
    array.close()
    context.close()

    return result
}


fun readDataFromTileDB_PREVIOUS(arrayName: String, idList: List<String>): List<Map<String, String>> {
    val context = Context()
    val result = mutableListOf<Map<String, String>>()

    // Open the array in read mode
    val array = Array(context, arrayName, QueryType.TILEDB_READ)

    // Prepare ID buffers
    val idBuffer = idList.joinToString(separator = "")
    val idOffsets = idList.runningFold(0) { acc, item -> acc + item.length }.dropLast(1).map { it.toLong() }.toLongArray()

    val idDataBuffer = io.tiledb.java.api.NativeArray(context, idBuffer, Datatype.TILEDB_STRING_ASCII)
    val idOffsetsBuffer = io.tiledb.java.api.NativeArray(context, idOffsets, Datatype.TILEDB_UINT64)

    // Prepare buffers for SampleName and Regions
    val sampleNameBuffer = io.tiledb.java.api.NativeArray(context, 1024, Datatype.TILEDB_STRING_ASCII) // Adjust size as needed
    val regionsBuffer = io.tiledb.java.api.NativeArray(context, 1024, Datatype.TILEDB_STRING_ASCII)   // Adjust size as needed

    // Prepare query
    val query = Query(array, QueryType.TILEDB_READ).apply {
        setLayout(Layout.TILEDB_UNORDERED)

        // Set ID buffers
        setDataBuffer("ID", idDataBuffer)
        setOffsetsBuffer("ID", idOffsetsBuffer)

        // Set attribute buffers
        setDataBuffer("SampleName", sampleNameBuffer)
        setDataBuffer("Regions", regionsBuffer)
    }

    // Submit query
    query.submit()


    // Process results
    val resultBufferElements = query.resultBufferElements()
    println("\nLCJ - resultBufferElements: $resultBufferElements\n")
    // Debug: Inspect raw buffers
    val rawSampleNameData = String((sampleNameBuffer.toJavaArray() as ByteArray))
    val rawRegionsData = String((regionsBuffer.toJavaArray() as ByteArray))

    println("\nLCJ - Raw SampleName Data: $rawSampleNameData\n")
    println("\nLCJ - Raw Regions Data: $rawRegionsData\n")

    val sampleNameLength = resultBufferElements["SampleName"]?.second?.toInt() ?: 0
    val regionsLength = resultBufferElements["Regions"]?.second?.toInt() ?: 0

    println("\nLCJ - sampleNameLength: $sampleNameLength , regionsLength: $regionsLength\n")

    // Extract results
    val sampleNames = String((sampleNameBuffer.toJavaArray() as ByteArray), 0, sampleNameLength).split("\u0000").filter { it.isNotEmpty() }
    val regions = String((regionsBuffer.toJavaArray() as ByteArray), 0, regionsLength).split("\u0000").filter { it.isNotEmpty() }

    // Combine results
    for (i in sampleNames.indices) {
        result.add(
            mapOf(
                "ID" to idList[i],
                "SampleName" to sampleNames[i],
                "Regions" to regions[i]
            )
        )
    }

    // Close resources
    query.close()
    array.close()
    context.close()

    return result
}




