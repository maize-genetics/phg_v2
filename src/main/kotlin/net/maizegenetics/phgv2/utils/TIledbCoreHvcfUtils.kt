package net.maizegenetics.phgv2.utils

import htsjdk.variant.vcf.VCFFileReader
import io.tiledb.java.api.*
import io.tiledb.java.api.Array
import io.tiledb.java.api.Constants.TILEDB_VAR_NUM
import org.apache.logging.log4j.LogManager
import java.io.File
import io.tiledb.java.api.TileDBObject

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

enum class HvcfAttributes {
    REFCHECKSUM, REGIONS, ID1, ID2
}

// ID is only a dimension of the alt header array, not the variants array
enum class HvcfDimensions {
    REF_RANGE, SAMPLE_NAME, ID
}

// THis returns a list of maps, where each map contains the fields ID, SampleName, Regions, RefChecksum, and RefRange.
// with the values for each field from the ALT headers in the VCF file.
// Each item on the list is the data from a single ALT header line.
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

fun parseTiledbVariantData(vcfReader:VCFFileReader): List<Map<String,String>> {
    val variantData = mutableListOf<Map<String, String>>()

    val sampleName = vcfReader.header.genotypeSamples.firstOrNull()

    // Iterate through VariantContext objects in the VCF
    vcfReader.iterator().forEach { variantContext ->
        val genotype = variantContext.getGenotype(sampleName)

        if (genotype.isNoCall) {
            // In the imputed hvcf files, all ref ranges are present, but there will be
            // no data in the GT field if that range is not represented in the sample.
            println("LoadHvcf: No GT call for sample $sampleName, variant at ${variantContext.contig}:${variantContext.start}")
            return@forEach
        }

        val altAlleles = variantContext.alternateAlleles
        val gtField = genotype.genotypeString // need this to determine how many alleles there are
        val alleles = gtField.split("[/|]".toRegex()) // Split on "/" or "|" for diploid

        // create the refRange mapping
        val chr = variantContext.contig
        val start = variantContext.start
        val end = variantContext.end
        val refRange = "${chr}:${start}-${end}"

        // Create the entry.  Only 1 if is haploid, 2 if is diploid
        val entry = mutableMapOf<String, String>()
        entry["RefRange"] = refRange
        entry["ID1"] = altAlleles[0].displayString
        entry["SampleName"] = "${sampleName}"
        entry["ID2"] = altAlleles[0].displayString // Will be changed below if there are 2 alleles

        // DO we want this to be 0, or do we want this to be
        // the same value as the ID1 field?
        // If there is only 1 allele, then the ID2 field should be the same as the ID1 field
        // when it is diploid.  But what about haploid - how do we know it is hapoid if we
        // fill in the ID2 field?
        if (alleles.size == 2 && altAlleles.size == 2) {
            entry["ID2"] = altAlleles[1].displayString // If there are 2 alt alleles, set ID2 to the second allele
        }
        variantData.add(entry)
    }
    return variantData
}

// If there is no path for the arrayName, the array will be created
// in the current directory.
// In both arrays I am setting the capacity to 5000.  This is based on a chatGPT suggestion
// for 70K ranges and 100 samples.  This may need to be adjusted.  It affects the sub-arrays
// and query efficiency.
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
        setCapacity(5000)
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
        setCapacity(5000) // picked based on 70K ranges, 100 samples - a chatGPT suggestion
    }

    // Create the array
    Array.create(arrayNameVariants, schema2)

    // Close resources
    schema2.close()
    domain2.close()
    context.close()

}

/**
 * THis function verifies the tiledb array exists in the folder provided.
 * Currently it is only checking for 1 of the arrays, and assume if 1 is there
 * and is good they both are.
 */
fun verifyHvcfArray(dbPath:String): Boolean {
    // Verify the folder exists and the hvcf_array exists in that folder
    // as a tiledb array
    val arrayName = "${dbPath}/alt_header_array"
    if (!File(dbPath).exists() || !File(dbPath).isDirectory()) {
        myLogger.warn("Folder $dbPath does not exist - creating.")
        File(dbPath).mkdirs()
        createTileDBCoreArrays(arrayName)
    } else {
        myLogger.info("Folder $dbPath exists, trying to open array name $arrayName")
        try {
            val array = Array(Context(), arrayName, QueryType.TILEDB_READ)
        } catch (exc: TileDBError) {
            myLogger.error("Array $arrayName does not exist or exists but is not a tiledb array.\n" +
                    "Please check your arrayName for accuracy and run  Initdb to set up the arrays if necessary.")
            throw exc
        }
    }
    return true

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
    //println("prepareVariableBuffer: size of data list = ${data.size}")

    // Create a running total of the lengths of the strings in the data list, dropping
    // the last element as that is not needed, then convert to a long array.
    val offsets = data.runningFold(0) { acc, item -> acc + item.length }
        .dropLast(1)
        .map { it.toLong() } // tiledb wants this data as TILEDB_UINT64, which is Kotlin long
        .toLongArray()

    //println("prepareVariableBuffer: size of offsets array = ${offsets.size}")
    // WIth variable length attributes, we need both a data buffer and an offsets buffer
    val dataBuffer = io.tiledb.java.api.NativeArray(context, concatenated, Datatype.TILEDB_STRING_ASCII)
    val offsetsBuffer = io.tiledb.java.api.NativeArray(context, offsets, Datatype.TILEDB_UINT64)

    return Pair(dataBuffer, offsetsBuffer)
}


// Function takes an array name and a list of maps, where each map contains
// the fields ID, SampleName, Regions, RefChecksum, and RefRange.
// The function writes the data to the TileDB array.
// This does not deal with the contexts/variants of the hvcf file, only the alt headers.
// TODO - this should be batched to handle large numbers of alt headers.
fun writeAltDataToTileDB(arrayName: String, altData: List<Map<String, String>>) {
    val context = Context()

    // Open the tiledb array in write mode
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

    // Prepare query
    val query = Query(array, QueryType.TILEDB_WRITE).apply {
        // We need to write with TILEDB_UNORDERED as the ALT headers are sorted by hapid, and
        // our tiledb dimensions are sorted by refRange and SampleName.  If we try to load
        // in TILEDB_GLOBAL_ORDER it fails because the refRanges are not in order.
        // But when reading the data, we can us TILEDB_GLOBAL_ORDER to retrieve the data in the order
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

    // THis is to help with fragmenting, and repositioning the tiles for efficiency
    val consolidateContext = Context()
    Array.consolidate(consolidateContext, arrayName) // Static method for consolidation
    consolidateContext.close()

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
        setDataBuffer("RefRange", refRangeBuffers.first)
        setOffsetsBuffer("RefRange", refRangeBuffers.second)

        setDataBuffer("SampleName", sampleBuffers.first)
        setOffsetsBuffer("SampleName", sampleBuffers.second)

        // Set buffers for fixed length attributes (fixed length attributes do not need offsets)
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
 */

fun queryWithStreaming_sampleNameIdByRefRange(arrayName: String, refRangeList: List<String>): Map<String,List<Map<String,String>>> {
    val context = Context()
    //val combinedResult = mutableListOf<Map<String, String>>()
    val resultMap = mutableMapOf<String, List<Map<String,String>>>()

    // Open the array in read mode
    val array = Array(context, arrayName, QueryType.TILEDB_READ)

    // Prepare buffers for RefRange (filtering dimension)
    val refRangeBuffers = prepareVariableBuffer(context, refRangeList)
    println("queryWithStreaming: refRangeList size = ${refRangeList.size},  refRangeBuffers.first size = ${refRangeBuffers.first.size}, second size = ${refRangeBuffers.second.size}")
    // print the refRangeBuffers.second as a long array
    val refRangeOffsetsArray = refRangeBuffers.second.toJavaArray() as LongArray
    println("queryWithStreaming: refRangeOffsetsArray = ${refRangeOffsetsArray.contentToString()}")

    // Set reasonable initial buffer sizes
    val sampleNameBufferSize = 8192
    val sampleNameOffsetsSize = 1024
    val idBufferSize = 8192
    val idOffsetsSize = 1024
    val refRangeBufferSize = 8192
    val refRangeOffsetsSize = 1024

    // Allocate initial buffers
    val sampleNameBuffer = io.tiledb.java.api.NativeArray(context, sampleNameBufferSize, Datatype.TILEDB_STRING_ASCII)
    val sampleNameOffsets = io.tiledb.java.api.NativeArray(context, sampleNameOffsetsSize, Datatype.TILEDB_UINT64)
    val idBuffer = io.tiledb.java.api.NativeArray(context, idBufferSize, Datatype.TILEDB_STRING_ASCII)
    val idOffsets = io.tiledb.java.api.NativeArray(context, idOffsetsSize, Datatype.TILEDB_UINT64)
    val refRangeBuffer = io.tiledb.java.api.NativeArray(context, refRangeBufferSize, Datatype.TILEDB_STRING_ASCII)
    val refRangeOffsets = io.tiledb.java.api.NativeArray(context, refRangeOffsetsSize, Datatype.TILEDB_UINT64)

    // Create the query
    val query = Query(array, QueryType.TILEDB_READ).apply {
        setLayout(Layout.TILEDB_UNORDERED) // GLOBAL_ORDER didn't get all the data, and associated it incorrectly

        // Set RefRange buffers
//        setDataBuffer("RefRange", refRangeBuffer)
//        setOffsetsBuffer("RefRange", refRangeOffsets)

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
//        val refRangeOffsetsArray = refRangeOffsets.toJavaArray() as LongArray
//        val refRangeRawData = String(refRangeBuffer.toJavaArray() as ByteArray)
//        val refRanges = refRangeOffsetsArray.take(numRefRanges).mapIndexedNotNull { index, offset ->
//            val end = if (index < refRangeOffsetsArray.size - 1) refRangeOffsetsArray[index + 1].toInt() else refRangeRawData.length
//            if (offset.toInt() < end) refRangeRawData.substring(offset.toInt(), end).trimEnd('\u0000') else null
//        }

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

        // Combine results into a list of map
//        val results = refRanges.zip(sampleNames.zip(ids)) { refRange, (sampleName, id) ->
//            mapOf(
//                "RefRange" to refRange,
//                "SampleName" to sampleName,
//                "ID" to id
//            )
//        }
//        combinedResult.addAll(results)

        // TRy this instead of above: No - doesn't help.  STill not getting the correct
        // results.  Not sure if the problem is reading the data or storing the data.
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
    // Above got ALL The ref ranges.  I don't know why since we gave it a buffer
    // with just a couple.  So now filter the results to just what we requested
    //val filteredResults = combinedResult.filter { it["RefRange"] in refRangeList }
    val filteredResults = resultMap.filterKeys { it in refRangeList }
    println("queryWithStreaming: number of  resultMap entries = ${resultMap.keys.size}, num filtered results = ${filteredResults.keys.size}")
    println("\nHere are all the values from the resultMap before filtering:")
    resultMap.forEach { (key, value) -> println("Key: $key, Value: $value") }
    //println("\nHere are all the values from the original combinedResult map before filtering:")
    // print the values in combinedResult.  These turn out to be the same as the resultMap
    //combinedResult.forEach { println(it) }

    return filteredResults
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
    println("querySampleNamesAndIDsByRefRange: RefRangeBuffers.first = ${refRangeBuffers.first}")
    println("querySampleNamesAndIDsByRefRange: RefRangeBuffers.second = ${refRangeBuffers.second}")

    // Prepare buffers for SampleName and ID
    val sampleNameBuffer = io.tiledb.java.api.NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII)
    val sampleNameOffsets = io.tiledb.java.api.NativeArray(context, 512, Datatype.TILEDB_UINT64)

    val idBuffer = io.tiledb.java.api.NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII)
    val idOffsets = io.tiledb.java.api.NativeArray(context, 512, Datatype.TILEDB_UINT64)

    // Create and configure the query
    val query = Query(array, QueryType.TILEDB_READ).apply {
        setLayout(Layout.TILEDB_GLOBAL_ORDER)

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

    // LCJ - raw data is NOT coming back correctly - I'm missing 2:1-1000
    println("RefRangeOffsets: ${refRangeOffsetsArray.contentToString()}")
    println("RefRangeRaw Data: $refRangeRawData")
    val refRanges = refRangeOffsetsArray.take(numRefRanges).mapIndexedNotNull { index, offset ->
        val end = if (index < refRangeOffsetsArray.size - 1) refRangeOffsetsArray[index + 1].toInt() else refRangeRawData.length
        if (offset.toInt() < end) refRangeRawData.substring(offset.toInt(), end).trimEnd('\u0000') else null
    }

    // This only shows 2 ref ranges:  they are both 1:1-1000, but one for LineA and one for LineB.
    // What happened to the other ref ranges?  I sent in 1:1-1000","2:1-1000
    println("querySampleNamesAndIDsByRefRange: Parsed RefRanges: $refRanges")
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






