package net.maizegenetics.phgv2.cli

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.enum
import net.maizegenetics.phgv2.utils.TileDBCoreHvcfUtils
import net.maizegenetics.phgv2.utils.TileDBCoreVariantQueries
import org.apache.logging.log4j.LogManager
import java.io.File

//Enum to hold the Query Types currently supported
enum class QueryType {
    DISTINCT_SAMPLES,
    DISTINCT_RANGES
}

//Enum to hold the Array Types currently supported
enum class ArrayType {
    ALT_HEADER,
    VARIANTS
}

/**
 * This class processing user queries for data stored in the tiledb core arrays.
 * There are 2 arrays:  the alt_header_array and the hvcf_variants_array
 * HVCF from imputed data will only have data stored to the hvcf_variants_array as those
 * files do not contain ALT header lines.
 *
 * Users may specify a list of ranges, a list of samples, and/or a list of attributes of
 * data they wish returned.
 */
class QueryHvcfArrays: CliktCommand(help = "Query tiledb core arrays for hvcf file data") {

    private val myLogger = LogManager.getLogger(QueryHvcfArrays::class.java)

    val dbPath by option(help = "Folder name under which TileDB datasets will be created. If this folder does not exist, it will be created.")
        .default("")
    val refRangeFile by option(help = "Full path to BED file formatted list of reference ranges to query.  Default is all")
        .default("")
    val queryType by option(help = "Type of query to perform.  Choices are:  distinctSamples, distinctRanges")
        .enum<QueryType>()
        .required()
    val arrayType by option(help = "Type of array to query.  Choices are:  altHeader, variants.  Default is variants")
        .enum<ArrayType>()
        .default(ArrayType.VARIANTS)

    //Leaving this commented as we will likely need this in the future
//    val sampleNames by option(help = "Full path to text file with sampleNames, 1 per line. Default is all")
//        .default("")

    val outputFile by option(help = "Full path to Output file for results.")
        .required()

    override fun run() {
        logCommand(this)

        // If the dbPath is blank, then use the current working directory
        val dbPath = if (dbPath.isBlank()) {
            System.getProperty("user.dir")
        } else {
            dbPath
        }

        // call method to handle query

        processUserQuery(dbPath, refRangeFile, queryType, arrayType, outputFile)
    }


    // This is initial query.  Only handling distinctSamples and distinctRanges
    // Should be expanded to handle other queries, and perhaps break processing
    // into separate functions.
    fun processUserQuery(dbPath:String, refRangeFile:String, queryType:QueryType, arrayType:ArrayType, outputFile:String) {
        myLogger.info("Processing user query: $queryType")
        // Read in the refRangeFile, sampleNames, and ids files if they are not blank
        // THis is a BED file - entries need to be translated to the format chrom:start-end
        val refRangeList = if (refRangeFile.isNotBlank()) {
            bufferedReader(refRangeFile).readLines().map { line ->
                val tokens = line.split("\t")
                val chrom = tokens[0]
                val start = tokens[1].toInt() + 1 // BEd is 0-based, our refRanges are 1-based.
                val end = tokens[2].toInt()
                "$chrom:$start-$end"
            }
        } else {
            null
        }

        // Verify the tiledbURI - an exception is thrown from verifyURI if the URI is not valid
        myLogger.info("QueryHvcfArrays: verifying array")
        val goodArray = TileDBCoreHvcfUtils.verifyHvcfArray(dbPath)
        myLogger.info("QueryHvcfArrays: goodArray: $goodArray")

        // Define the array to query
        // Using an enum means we don't need a else catch all
        val queryArray = when(arrayType) {
            ArrayType.ALT_HEADER -> {
                "$dbPath/alt_header_array"
            }
            ArrayType.VARIANTS -> {
                "$dbPath/hvcf_variants_array"
            }
        }
        // query based on queryType
        // Using an enum means we don't need a else catch all
        when (queryType) {
            QueryType.DISTINCT_SAMPLES -> {
                val distinctSamples = TileDBCoreVariantQueries.queryDistinctSampleNames(queryArray)
                // Write the results from both the altHeaderArray and the variantsArray to the output file
                File(outputFile).writeText(distinctSamples.joinToString("\n"))
            }
            QueryType.DISTINCT_RANGES -> {
                // The ranges shoudl be identical in the altHeaderArray and the variantsArray
                // so only 1 query is needed
                val distinctRanges = TileDBCoreVariantQueries.queryDistinctRefRanges(queryArray)
                myLogger.info("Distinct refRanges: $distinctRanges\n")
                // Write the results to the output file
                File(outputFile).writeText(distinctRanges.joinToString("\n"))
            }
        }
        myLogger.info(" query written to file $outputFile")
    }
}