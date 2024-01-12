package net.maizegenetics.phgv2.brapi.service

import com.typesafe.config.ConfigFactory
import io.ktor.server.config.*
import model.TokenPagination
import model.Variant
import net.maizegenetics.phgv2.api.ReferenceRange
import org.apache.logging.log4j.LogManager
import org.ehcache.config.builders.CacheConfigurationBuilder
import org.ehcache.config.builders.CacheManagerBuilder
import org.ehcache.config.builders.ResourcePoolsBuilder
import org.ehcache.core.internal.statistics.DefaultStatisticsService
import org.ehcache.core.spi.service.StatisticsService
import java.io.File
import java.util.*
import kotlin.collections.ArrayList
import kotlin.collections.HashSet

val statisticsVariantsService: StatisticsService = DefaultStatisticsService()
private val config = HoconApplicationConfig(ConfigFactory.load())


/**
 * Class to handle all the Variant creation.
 */

// NOTE: when using the variantCache to pull reference range data,
// use variantDbId-1.  The variantDbIds are 1-based:  actual reference range ids
// These are stored in order, in a list that is 0-based.
private val cacheManager = CacheManagerBuilder.newCacheManagerBuilder()
    .using(statisticsVariantsService)
    .withCache(
        "variantCache",
        CacheConfigurationBuilder.newCacheConfigurationBuilder(
            String::class.java,
            ArrayList::class.java,
            ResourcePoolsBuilder.heap(1).build()
        )
    )
    .build(true)

private val variantCache = cacheManager.getCache("variantCache",String::class.java, ArrayList::class.java)

class VariantsService {
    private val myLogger = LogManager.getLogger(VariantsService::class.java)

    init {
        println("VariantsService:init -  heapSIze before create variants cache:")
        //Sizeof.printMemoryUse()  // copy TASSEL function to phgv2?

        // Search the tiledbURI/reference directory for bed files.
        // there should only be 1.  If there are more, we will use the first one.
        // return the name of that file
        val bedFile = File("${tiledbURI}/reference/").walk().filter { it.name.endsWith(".bed") }.toList()[0]
        val referenceRanges = createRefRangesForCache(bedFile) as ArrayList<ReferenceRange>
        if (referenceRanges != null) {
            // This list of reference ranges as stored with key "haplotype"
            variantCache.put("all",referenceRanges)
            val ehCacheStat = statisticsVariantsService.getCacheStatistics("variantCache");
            // using println until we get myLogger working
            println("variant (reference Range) cache size: ${ehCacheStat.getTierStatistics().get("OnHeap")?.getMappings()}")
            myLogger.info("variant (reference Range) cache size: ${ehCacheStat.getTierStatistics().get("OnHeap")?.getMappings()}")
        }
        println("VariantsService:init -  heapSize after create variants cache:")
        //Sizeof.printMemoryUse()
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
                    start = fields[1].toInt(),
                    end = fields[2].toInt(),
                )
                referenceRange
            }
        return referenceRanges

    }

    // This method will return a list of variants for a given page.  The page is defined by the
    // next index into the ReferenceRange list.  The pageSize is the number of variants to return.
    // Currently, there is only a single group named "all".  This will be used to get the
    // reference ranges from the cache.  having the groupName parameter allows for this to change.
    fun generateVariantsListFromCache(currentPageToken:Int, pageSize:Int, groupName:String): Pair<TokenPagination,List<Variant>> {
        //check for cache - if it doesn't exist, create it
        // remove println when myLogger is working
        //var referenceRanges: ArrayList<ReferenceRange>? = null
        var referenceRanges: List<ReferenceRange>? = null
        if (variantCache.containsKey(groupName)) {
            myLogger.info("generateVariantsListFromCache - getting reference/variants from cache")
            println("generateVariantsListFromCache - getting reference/variants from cache")
            referenceRanges = variantCache[groupName] as ArrayList<ReferenceRange>

        } else {
            myLogger.info("generateVariantsListFromCache - no variantCache - creating it, heapSIze before create cache:")
            println("generateVariantsListFromCache - no variantCache - creating it, heapSIze before create cache:")
            //Sizeof.printMemoryUse()
            val bedFile = File("${tiledbURI}/reference/").walk().filter { it.name.endsWith(".bed") }.toList()[0]
            referenceRanges = createRefRangesForCache(bedFile) as ArrayList<ReferenceRange>

            variantCache.put(groupName, referenceRanges)
            myLogger.info("generateVariantsListFromCache: heapSize after create variant cache:")
            println("generateVariantsListFromCache: heapSize after create variant cache:")
            //Sizeof.printMemoryUse()
        }

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
        val endId = startId + pageSize -1; // endId is inclusive
        val totalCount = referenceRanges.size

        // Now we pull from the sortedRangeList just those list entries that are between the startId and endId,
        // and write them to a new list
        // LCJ - manually check this - we don't need a variant, we just want the list entries.
        val variants = sortedRangeList.subList(startId,endId+1).map { range ->
            Variant(
                referenceName = range.contig,
                start = range.start,
                end = range.end,
                variantDbId = "${range.contig}:${range.start}:${range.end}",
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

        // Plus 2:  Because we dropped 1 to go to 0-based for the referenceRanges list.
        // Now we need to add that back to get to 1-based, plus add 1 more to get to the
        // first variant beyond what we've already processed.
        var nextPageToken: String? = (endId+2).toString()
        if (endId > referenceRanges.size) {
            nextPageToken = null
        }

        var pagination = TokenPagination(pageSize=pageSize, nextPageToken=nextPageToken, currentPageToken=currentPageToken.toString(), totalCount=totalCount, totalPages=totalPages)
        return Pair<TokenPagination, List<Variant>>(pagination,variants)

    }
}