package net.maizegenetics.phgv2.brapi.service

import com.google.common.cache.CacheBuilder
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.brapi.model.TokenPagination
import net.maizegenetics.phgv2.brapi.model.Variant
import net.maizegenetics.phgv2.brapi.utilities.BrAPIConfig
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.*


/**
 * Class to handle all the Variant creation.
 */

private val variantCache = CacheBuilder.newBuilder()
    .maximumSize(100)
    .build<String, List<ReferenceRange>>()

class VariantsService {
    private val myLogger = LogManager.getLogger(VariantsService::class.java)

    init {

        // Search the tiledbURI/reference directory for bed files.
        // there should only be 1.  If there are more, we will use the first one.
        // return the name of that file
        val bedFileList = File("${BrAPIConfig.tiledbURI}/reference/").walk().filter { it.name.endsWith(".bed") }.toList()
        if ( bedFileList.size < 1) {
            // The bedfile is copied to tiledbURI/references when CreateRefVcf is run.
            // When running specific juint tests (e.g. ServerInfoTest) which do not run CreateRefVcf
            // as a prerequisite, compilation will fail on this init because there is no bedfile.
            // This conditional takes care of that problem.
            // Need to think through if this is a general problem or just a testing problem.
            myLogger.error("VariantsService:init - no bed files found in ${BrAPIConfig.tiledbURI}/reference/")
            println("VariantsService:init - no bed files found in ${BrAPIConfig.tiledbURI}/reference/")

        } else {
            val bedFile = bedFileList[0]
            val referenceRanges = createRefRangesForCache(bedFile) as ArrayList<ReferenceRange>
            // This list of reference ranges as stored with key "haplotype"
            variantCache.put("all", referenceRanges)
        }

    }

    // Use the tiledbURI defined in SamplesService to get the bed file.  From this,
    // we will create the reference range objects
    fun createRefRangesForCache(bedFile:File): List<ReferenceRange> {

        // Bed files may be white-space or tab-delimited.  They could contain header lines that
        // begin with either "track" or "browser".  We will ignore these lines.
        // Read the lines from file bedFile , filter out the header lines, create a ReferenceRange
        // object from each of the other lines and put onto a list.
        val referenceRanges = bedFile.readLines().filter { !it.startsWith("track") && !it.startsWith("browser") }
            .map { range ->
                val fields = range.split("\\s+".toRegex())
                val referenceRange = ReferenceRange(
                    contig = fields[0],
                    start = fields[1].toInt() + 1, // because bed is 0-based, but ref ranges are 1-based
                    end = fields[2].toInt(), // bed file end is exclusive, no need to add 1
                )
                referenceRange
            }
        return referenceRanges

    }

    // Function to query the cache for a set of reference ranges.
    // If the cache doesn't exist, create it.
    private fun getReferenceRanges(groupName:String): List<ReferenceRange>? {
        var referenceRanges = variantCache.getIfPresent(groupName)
        if (referenceRanges == null) {
            val bedFileList = File("${BrAPIConfig.tiledbURI}/reference/").walk().filter { it.name.endsWith(".bed") }.toList()
            if ( bedFileList.size < 1) {
                // The bedfile is copied to tiledbURI/references when CreateRefVcf is run.
                // When running specific juint tests (e.g. ServerInfoTest) which do not run CreateRefVcf
                // as a prerequisite, compilation will fail on this init because there is no bedfile.
                // This conditional takes care of that problem.
                // Need to think through if this is a general problem or just a testing problem.
                myLogger.error("VariantsService:init - no bed files found in ${BrAPIConfig.tiledbURI}/reference/")
                myLogger.error("VariantsService:init - please ensure CreateRefVcf has been run and the bed file is in ${BrAPIConfig.tiledbURI}/reference/")

            } else {
                val bedFile = bedFileList[0]
                referenceRanges = createRefRangesForCache(bedFile) as ArrayList<ReferenceRange>
                // This list of reference ranges as stored with key "haplotype"
                variantCache.put("all", referenceRanges)
            }
        }
        return referenceRanges
    }

    // This method will return a list of variants for a given page.  The page is defined by the
    // next index into the ReferenceRange list.  The pageSize is the number of variants to return.
    // Currently, there is only a single group named "all".  This will be used to get the
    // reference ranges from the cache.  having the groupName parameter allows for this to change.
    fun generateVariantsListFromCache(currentPageToken:Int, pageSize:Int, groupName:String = "all"): Pair<TokenPagination,List<Variant>> {
        //check for cache - if it doesn't exist, create it
        // remove println when myLogger is working
        //var referenceRanges: ArrayList<ReferenceRange>? = null
        val referenceRanges = getReferenceRanges(groupName)

        // if no reference ranges were found, return empty list
        if (referenceRanges == null) {
            return Pair(TokenPagination(pageSize = pageSize, nextPageToken = null, currentPageToken = currentPageToken.toString(), totalCount = 0), emptyList())
        }
        // This calculation rounds up
        val totalPages = (referenceRanges.size + pageSize - 1)/pageSize
        if (currentPageToken >= referenceRanges.size) {
            return Pair(TokenPagination(pageSize = pageSize, nextPageToken = null, currentPageToken = currentPageToken.toString(), totalCount = totalPages), emptyList())
        }

        // Return paginated list of variants.
        // convert the ArrayList to a list then sort it
        // With this, we can use the index into the list as the page tokens
        val sortedRangeList = referenceRanges.toList()
        Collections.sort(sortedRangeList)

        val startId = currentPageToken-1 // because tokens are 1 based, but list is 0-based
        val endId = if(startId + pageSize  < referenceRanges.size) startId + pageSize else referenceRanges.size  // endId is exclusive

        val totalCount = referenceRanges.size

        // Pull from the sortedRangeList just those list entries that are between the startId and endId,
        // and write them to a new list
        println("\nLCJ VariantsService:generateVariantsListFromCache - startId: $startId, endId: $endId \n")
        // Note: subList the endId is exclusive, accounted for when creating the endId above
        val variants = sortedRangeList.subList(startId,endId).map { range ->
            Variant(
                referenceName = range.contig,
                start = range.start,
                end = range.end,
                variantDbId = "${range.contig}:${range.start}-${range.end}",
                variantType = "REF_RANGE",
                referenceBases = "",
                alternateBases = emptyList(),
                filtersApplied = false,
                filtersFailed = emptyList(),
                filtersPassed = true,
                variantNames = emptyList(),
                variantSetDbId = emptyList(),
                additionalInfo = emptyMap(),
                ciend = emptyList(),
                cipos = emptyList(),
                created = null,
                svlen = range.end - range.start + 1,
                updated = null
            )
        }

        // add 1 to get to the next page token.  If the endId is the last reference range,
        var nextPageToken: String? = (endId+1).toString()
        if (endId > referenceRanges.size) {
            nextPageToken = null
        }

        var pagination = TokenPagination(pageSize=pageSize, nextPageToken=nextPageToken, currentPageToken=currentPageToken.toString(), totalCount=totalCount, totalPages=totalPages)
        return Pair<TokenPagination, List<Variant>>(pagination,variants)

    }

    // Find the data for a specific Variant (ie referenceRange)
    fun generateVariantFromID(variantDbId:String, pageToken:Int, pageSize:Int, groupName:String = "all"):Variant? {
        var referenceRanges = getReferenceRanges(groupName)

        // if no reference ranges were found, return empty list
        if (referenceRanges == null) {
            //  Return a single variant, or null?
            return null

        } else {
            val contig= variantDbId.split(":")[0]
            val start = variantDbId.split(":")[1].split("-")[0].toInt()
            val end = variantDbId.split(":")[1].split("-")[1].toInt()

            // Filter the reference ranges to find one that matches the variantDbId
            // if one is found, create a Variant object from it
            // if not found, return null
            val variant = referenceRanges.filter{it.contig == contig && it.start == start && it.end == end}.map { range ->
                Variant(
                    referenceName = range.contig,
                    start = range.start,
                    end = range.end,
                    variantDbId = "${range.contig}:${range.start}-${range.end}",
                    variantType = "REF_RANGE",
                    referenceBases = "",
                    alternateBases = emptyList(),
                    filtersApplied = false,
                    filtersFailed = emptyList(),
                    filtersPassed = true,
                    variantNames = emptyList(),
                    variantSetDbId = emptyList(),
                    additionalInfo = emptyMap(),
                    ciend = emptyList(),
                    cipos = emptyList(),
                    created = null,
                    svlen = range.end - range.start + 1,
                    updated = null
                )
            }

            if (variant.size == 0) {
                return null
            } else {
                return variant[0]
            }
        }
    }
}