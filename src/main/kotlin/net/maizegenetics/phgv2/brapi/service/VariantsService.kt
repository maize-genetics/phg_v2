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
        // object from each of hte other lines and put onto a list.
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

    fun generateVariantsListFromCache(currentPageToken:Int, pageSize:Int): Pair<TokenPagination,List<Variant>> {

        // LCJ - fix this !! is dummy code
        return Pair(TokenPagination(pageSize = pageSize, nextPageToken = null, currentPageToken = currentPageToken.toString(), totalCount = 0), emptyList())
    }
}