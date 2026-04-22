package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import org.apache.logging.log4j.LogManager
import java.io.File

class MappingCountTableQc: CliktCommand("Read mapping count QC") {

    val myLogger = LogManager.getLogger(ReadMappingCountQc::class.java)

    val hvcfDir by option(help = "Directory with Haplotype VCF files")
        .required()

    val readMappingDir by option(help = "Read mapping directory")
        .required()

    val outputDir by option(help = "Output directory")
        .required()


    /**
     * Entry point for the CLI command. Delegates to [buildMappingTables] using
     * the parsed option values.
     */
    override fun run() {
        myLogger.info("Read mapping count QC")

        buildMappingTables(hvcfDir, File(readMappingDir), File(outputDir))
    }

    /**
     * Builds a [HaplotypeGraph] from [hvcfDir] and pre-computes the two lookup
     * maps needed for count conversion, then processes every `*_readMapping.txt`
     * file found (recursively) under [readMappingDir], writing one counts table
     * per file to [outputDir].
     */
    fun buildMappingTables(hvcfDir: String, readMappingDir: File, outputDir: File) {
        myLogger.info("Building mapping tables")
        myLogger.info("Building Haplotype Graph")
        val graph = HaplotypeGraph(hvcfDir)

        val hapIdAndRangeToSampleGametes = graph.hapIdAndRefRangeToSampleGametes()
        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()


        //Walk through each readMapping file in ReadMappingDir
        readMappingDir.walkTopDown().filter { it.isFile }
            .filter { it.name.endsWith("_readMapping.txt") }
            .forEach { processSingleReadMappingFile(graph,it, outputDir, hapIdAndRangeToSampleGametes, hapIdToRefRangeMap) }
    }

    /**
     * Processes a single readMapping file end-to-end:
     * 1. Determines the output file path via [buildCountOutputFile].
     * 2. Parses the readMapping file with [AlignmentUtils.importReadMapping].
     * 3. Converts the raw hapId set → count entries into per-(range, sampleGamete)
     *    counts via [convertReadMappingToHapIdCounts].
     * 4. Writes the resulting count table via [writeOutCounts].
     *
     * @param graph the full haplotype graph built from the HVCF directory
     * @param readMappingFile the `*_readMapping.txt` file to process
     * @param outputDir directory where the counts file will be written
     * @param hapIdAndRangeToSampleGametes pre-built map of (ReferenceRange, hapId) → SampleGametes
     * @param hapIdToRefRangeMap pre-built map of hapId → list of ReferenceRanges it appears in
     */
    fun processSingleReadMappingFile(graph: HaplotypeGraph, readMappingFile: File, outputDir: File,
                                     hapIdAndRangeToSampleGametes: Map<Pair<ReferenceRange,String>,List<SampleGamete>>,
                                     hapIdToRefRangeMap: Map<String,List<ReferenceRange>>) {
        myLogger.info("Processing read mapping file: ${readMappingFile.name}")

        val outputFile = buildCountOutputFile(outputDir.absolutePath, readMappingFile)

        val readMapping = AlignmentUtils.importReadMapping(readMappingFile.absolutePath)

        val hapIdCounts = convertReadMappingToHapIdCounts(readMapping, hapIdAndRangeToSampleGametes, hapIdToRefRangeMap)

        writeOutCounts(hapIdCounts, graph,outputFile)
    }

    /**
     * Derives the output file path for a given read mapping file by replacing
     * the `_readMapping` suffix in the base name with `_mappingCounts` and
     * placing the result in [outputDir].
     *
     * For example, `LineA_1_readMapping.txt` → `<outputDir>/LineA_1_mappingCounts.txt`.
     */
    fun buildCountOutputFile(outputDir: String, readMappingFile: File): File {
        //Remove the _readMapping.txt and add _mappingCounts.txt
        val outputFileName = readMappingFile.nameWithoutExtension.replace("_readMapping", "_mappingCounts") + ".txt"
        return File(outputDir, outputFileName)
    }

    /**
     * Converts a raw readMapping (hapId set → read count) into a per-
     * (ReferenceRange, SampleGamete) count map suitable for the output table.
     *
     * For each entry in [readMappings]:
     * - The consensus [ReferenceRange] for the hapId set is resolved via
     *   [findRefRangeForHapIdSet]. Sets that span multiple ranges (i.e. whose
     *   hapIds disagree on range) cannot be attributed to a single range and
     *   are skipped with a warning.
     * - Every [SampleGamete] that carries any hapId in the set at that range
     *   has its count incremented by the read count for that entry.
     *
     * @param readMappings map of hapId set → number of reads that mapped to it
     * @param hapIdAndRangeToSampleGametes lookup of (ReferenceRange, hapId) → SampleGametes
     * @param hapIdToRefRangeMap lookup of hapId → ReferenceRanges it belongs to
     * @return map of (ReferenceRange, SampleGamete) → total read count
     */
    fun convertReadMappingToHapIdCounts(readMappings: Map<List<String>,Int>,
                                        hapIdAndRangeToSampleGametes: Map<Pair<ReferenceRange,String>,List<SampleGamete>>,
                                        hapIdToRefRangeMap: Map<String,List<ReferenceRange>>): Map<Pair<ReferenceRange, SampleGamete>, Int> {
        val finalCountMap = mutableMapOf<Pair<ReferenceRange, SampleGamete>, Int>()

        //go through each hapIdSet from the readMappings.
        // Determine which reference range it belongs to.
        // Then walk through the hapIds in the set and increase the count for that RefRange+hapId pair
        for((hapIdSet, count) in readMappings) {
            val topRefRangeHit = findRefRangeForHapIdSet(hapIdSet, hapIdToRefRangeMap)
            if(topRefRangeHit.contig == "UNKNOWN") {
                myLogger.warn("Haplotype $hapIdSet with count $count unable to resolve reference range.  " +
                        "Make sure the HVCFs used match what were used for building the ropebwt index")
                continue
            }

            //build the list of sample gametes to increase
            hapIdSet.flatMap { hapIdAndRangeToSampleGametes[Pair(topRefRangeHit, it)]!! }
                .forEach {
                    //increment the count map
                    finalCountMap[Pair(topRefRangeHit, it)] = finalCountMap.getOrDefault(Pair(topRefRangeHit, it), 0) + count
                }

        }
        return finalCountMap
    }

    /**
     * Determines the single [ReferenceRange] that all hapIds in [hapIdSet]
     * belong to.
     *
     * Each hapId is looked up in [hapIdToRefRange] and its associated ranges
     * are tallied. If every hapId in the set maps to the same range (i.e. the
     * highest-frequency range is hit exactly `hapIdSet.size` times), that range
     * is returned. Otherwise — meaning the hapIds span more than one reference
     * range — `ReferenceRange("UNKNOWN", -1, -1)` is returned to signal that
     * the set is ambiguous and should be skipped.
     *
     * @param hapIdSet list of haplotype checksums from a single readMapping entry
     * @param hapIdToRefRange map of hapId → list of ReferenceRanges it appears in
     * @return the consensus [ReferenceRange], or an UNKNOWN sentinel if ambiguous
     */
    fun findRefRangeForHapIdSet(hapIdSet: List<String>, hapIdToRefRange: Map<String,List<ReferenceRange>>): ReferenceRange {
        //Simple algorithm is just to walk through and count how many times each refRange is hit.
        // Take the max and check to see if the count equals hapIdSet size
        val refRangeCounts = hapIdSet.flatMap { hapId ->
            if(!hapIdToRefRange.containsKey(hapId)) {
               myLogger.info("Haplotype $hapId not in hapIdToRefRange map")
            }
            hapIdToRefRange[hapId]!!
        }.groupingBy { it }.eachCount()

        //get the max
        val mostHitRefRange = refRangeCounts.maxByOrNull { it.value }!!
        return if(mostHitRefRange.value != hapIdSet.size) {
            ReferenceRange("UNKNOWN", -1,-1)
        } else {
            mostHitRefRange.key
        }
    }

    /**
     * Writes the final mapping count table to [outputFile] as a tab-delimited
     * text file.
     *
     * The header row is `RefRange` followed by every [SampleGamete] in the
     * graph (sorted). Each subsequent row corresponds to one [ReferenceRange]
     * (in sorted order) and contains the read count attributed to each
     * sampleGamete at that range, defaulting to 0 when no reads mapped there.
     *
     * @param hapIdCounts map of (ReferenceRange, SampleGamete) → total read count,
     *   as produced by [convertReadMappingToHapIdCounts]
     * @param graph the haplotype graph used to enumerate all ranges and sampleGametes
     * @param outputFile destination file; created or overwritten
     */
    fun writeOutCounts(hapIdCounts: Map<Pair<ReferenceRange, SampleGamete>, Int>, graph: HaplotypeGraph, outputFile: File) {
        val sampleGametes = graph.sampleGametesInGraph()
        val ranges = graph.ranges().sorted()
        bufferedWriter(outputFile.absolutePath).use { writer ->
            //writing out refRange and all the sampleGametes
            writer.write("RefRange\t${sampleGametes.joinToString("\t")}\n")
            for(range in ranges) {
                val countsForRange = sampleGametes.map { sampleGamete ->
                    hapIdCounts.getOrDefault(Pair(range, sampleGamete), 0)
                }
                writer.write("$range\t${countsForRange.joinToString("\t")}\n")
            }
        }
    }


}