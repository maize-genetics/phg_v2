package net.maizegenetics.phgv2.utils

import io.tiledb.java.api.*
import io.tiledb.java.api.Array
import org.jetbrains.kotlinx.dataframe.DataFrame
import org.jetbrains.kotlinx.dataframe.api.*
import io.tiledb.java.api.Query
import org.apache.logging.log4j.LogManager

/**
 * This class contains queries to be run against a TileDB core hvcf_variant_arrays
 * Each test takes an array name - that array should be an hvcf_variant_array, populated
 * with hvcf file data created from both the assembly alignment pipeline and the imputation pipeline.
 * THere is a separate class for the alt_header_array queries
 *
 * TODO Refactor some of these functions as they could likely be simplified
 */
class TileDBCoreVariantQueries {
    companion object {
        val myLogger = LogManager.getLogger(TileDBCoreVariantQueries::class.java)


        /**
         * Convert the tiledb data to a Kotlin DataFrame for easier processing
         * and writing to a file.  THe format here may not be what we want - I'm thinking
         * there needs to be more processing into columns to get this accurate.
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
         * query distinct sample names from the tiledb array/  Could  use this query for
         * either array, as SampleName is a dimension named in both the alt_header_array and the
         * hvcf_variants_array.
         */
        fun queryDistinctSampleNames(arrayName: String): Set<String> {
            val context = Context()
            val allSamples = mutableSetOf<String>()

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

            do {
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

                allSamples.addAll(sampleNames)
            } while (query.queryStatus == QueryStatus.TILEDB_INCOMPLETE)

            // Close resources
            query.close()
            array.close()
            context.close()

            return allSamples
        }

        /**
         * this function takes an array name and queries for a list of distinct ref ranges
         * found in that array.  It could be run on either the alt_header_array or the hvcf_variants_array
         * as both have RefRange as a dimension.  This is a simple query that returns a set of strings.
         */
        fun queryDistinctRefRanges(arrayName: String): Set<String>{
            val context = Context()

            // Open the TileDB array in read mode
            myLogger.info("queryDistinctRefRanges: arrayName = $arrayName")
            val array = Array(context, arrayName, QueryType.TILEDB_READ)

            val allRanges = mutableSetOf<String>()
            // Prepare buffers for SampleName (variable-length strings)
            val refRangeBuffer = NativeArray(context, 4096, Datatype.TILEDB_STRING_ASCII) // Adjust size as needed
            val refRangeOffsets = NativeArray(context, 512, Datatype.TILEDB_UINT64) // Adjust size as needed

            // Set up the query
            val query = Query(array, QueryType.TILEDB_READ).apply {
                setLayout(Layout.TILEDB_GLOBAL_ORDER) // Read in global order
                setDataBuffer("RefRange", refRangeBuffer)
                setOffsetsBuffer("RefRange", refRangeOffsets)
            }

            do {

                query.submit()

                // Extract the number of results
                val numRefRanges = query.resultBufferElements()["RefRange"]?.first?.toInt() ?: 0

                // Extract Range data
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
                allRanges.addAll(refRanges)
            } while (query.queryStatus == QueryStatus.TILEDB_INCOMPLETE)

            // Close resources
            query.close()
            array.close()
            context.close()

            return allRanges

        }

        /**
         * This returns the SampleName, ID1 and ID2 from the hvcf_variants_array for a list of RefRanges
         * The results is a map of RefRange to a list of maps where each map contains the fields SampleName, ID1, and ID2.
         *
         * It uses streaming to get all values.  This needs to be tested on a larger dataset than the test data.
         * TODO: test on BL01 with the phgv2_commonSept2024 data loaded to tiledb core arrays.
         */
        fun queryVariantArrayByRefRange(
            arrayName: String,
            refRangeList: List<String>?
        ): Map<String, List<Map<String, String>>> {

            val resultMap = mutableMapOf<String, List<Map<String, String>>>()
            val context = Context()

            // Open the TileDB array in read mode
            val array = Array(context, arrayName, QueryType.TILEDB_READ)

            if (refRangeList != null) {
                myLogger.info("queryByRefRanges: refRangeList size = ${refRangeList.size}")
                refRangeList.forEach { myLogger.info("queryByRefRanges: refRange = $it") }
            }


            // Prepare buffers for RefRange using prepareVariableBuffer
            val refRangeBuffers = if (refRangeList.isNullOrEmpty()) {
                val buffer = NativeArray(context, 8192, Datatype.TILEDB_STRING_ASCII)
                val offset = NativeArray(context, 1024, Datatype.TILEDB_UINT64)
                Pair(buffer,offset)
            } else {
                TileDBCoreHvcfUtils.prepareVariableBuffer(context, refRangeList)
            }
            val refRangeOffsetsArray = refRangeBuffers.second.toJavaArray() as LongArray
            myLogger.info("variantArray: refRangeOffsetsArray = ${refRangeOffsetsArray.contentToString()}")

            // Allocate buffers
            val sampleNameBufferSize = 8192
            val sampleNameOffsetsSize = 1024
            val idBufferSize = 8192

            val sampleNameBuffer = NativeArray(context, sampleNameBufferSize, Datatype.TILEDB_STRING_ASCII)
            val sampleNameOffsets = NativeArray(context, sampleNameOffsetsSize, Datatype.TILEDB_UINT64)
            val id1Buffer = NativeArray(context, idBufferSize, Datatype.TILEDB_STRING_ASCII)
            val id2Buffer = NativeArray(context, idBufferSize, Datatype.TILEDB_STRING_ASCII)

            // Create the query
            // This gets ALL data. To have tiledb do the slicing we need to do an "addRange" to the query.
            // That code appears to be available in a Java API version that is not yet available on Maven
            // So all data is returned and then filtered in the code below.
            val query = Query(array, QueryType.TILEDB_READ).apply {
                setLayout(Layout.TILEDB_UNORDERED)

                setDataBuffer("RefRange", refRangeBuffers.first)
                setOffsetsBuffer("RefRange", refRangeBuffers.second)

                setDataBuffer("SampleName", sampleNameBuffer)
                setOffsetsBuffer("SampleName", sampleNameOffsets)
                setDataBuffer("ID1", id1Buffer) // no offsets buffers for fixed size attributes
                setDataBuffer("ID2", id2Buffer)
            }

            // Iterate until the query is complete
            do {
                query.submit()

                // Determine the number of valid entries
                val numSampleNames = query.resultBufferElements()["SampleName"]?.first?.toInt() ?: 0
                val numRefRanges = query.resultBufferElements()["RefRange"]?.first?.toInt() ?: 0

                // Extract data from buffers
                val refRangeRawData = String(refRangeBuffers.first.toJavaArray() as ByteArray)
                val refRangeOffsetsArrayResult = refRangeBuffers.second.toJavaArray() as LongArray
                val refRanges = refRangeOffsetsArrayResult.take(numRefRanges).mapIndexedNotNull { index, offset ->
                    val end = if (index < refRangeOffsetsArrayResult.size - 1) refRangeOffsetsArrayResult[index + 1].toInt()
                    else refRangeRawData.length
                    if (offset.toInt() < end) refRangeRawData.substring(offset.toInt(), end).trimEnd('\u0000') else null
                }

                val sampleNameRawData = String(sampleNameBuffer.toJavaArray() as ByteArray)
                val sampleNameOffsetsArray = sampleNameOffsets.toJavaArray() as LongArray
                val sampleNames = sampleNameOffsetsArray.take(numSampleNames).mapIndexedNotNull { index, offset ->
                    val end = if (index < sampleNameOffsetsArray.size - 1) sampleNameOffsetsArray[index + 1].toInt()
                    else sampleNameRawData.length
                    if (offset.toInt() < end) sampleNameRawData.substring(offset.toInt(), end).trimEnd('\u0000') else null
                }

                val id1RawData = String(id1Buffer.toJavaArray() as ByteArray)
                val id1s = id1RawData.chunked(32) // Fixed size of 32 bytes

                val id2RawData = String(id2Buffer.toJavaArray() as ByteArray)
                val id2s = id2RawData.chunked(32) // Fixed size of 32 bytes

                // Associate the data
                refRanges.forEachIndexed { index, refRange ->
                    val sampleName = sampleNames.getOrNull(index) ?: ""
                    val id1 = id1s.getOrNull(index) ?: ""
                    val id2 = id2s.getOrNull(index) ?: ""

                    resultMap[refRange] = resultMap.getOrDefault(refRange, emptyList()) + mapOf(
                        "SampleName" to sampleName,
                        "ID1" to id1,
                        "ID2" to id2
                    )
                }

            } while (query.queryStatus == QueryStatus.TILEDB_INCOMPLETE)

            // Close resources
            query.close()
            array.close()
            context.close()

            // Filter and return results
            if (refRangeList.isNullOrEmpty()) {
                return resultMap
            } else {
                return resultMap.filterKeys { it in refRangeList }
            }

        }

    }
}


