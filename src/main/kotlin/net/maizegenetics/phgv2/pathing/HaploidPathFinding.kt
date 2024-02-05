package net.maizegenetics.phgv2.pathing

import biokotlin.seqIO.NucSeqIO
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.boolean
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.int
import htsjdk.variant.vcf.VCFAltHeaderLine
import htsjdk.variant.vcf.VCFHeaderVersion
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.Job
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.channels.ReceiveChannel
import kotlinx.coroutines.channels.SendChannel
import kotlinx.coroutines.launch
import kotlinx.coroutines.runBlocking
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.*
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * Finds the most likely path through a HaplotypeGraph based on read mapping data. Input: a keyfile with taxa names
 * and read mapping file paths and a haplotype graph. Read mappings must be in the form of
 * haplotype list -> count of reads mapping to the haplotypes. The haplotype graph must be the one used to map reads or
 * a subset of that restricted to specific taxa and/or ReferenceRanges.
 *
 * Steps:
 * 1. Read the keyfile
 * 1. Build the HaplotypeGraph
 * 2. For each entry in the keyfile (multithreaded)
 *      a. Get the read mapping data
 *      b. Use Viterbi algorithm to find the best path
 *      c. Store the path (write hvcf)
 *
 *
 */
class HaploidPathFinding : CliktCommand(help = "Impute haploid paths") {
    val pathKeyfile by option(help = "tab-delimited file with first two columns: SampleName, ReadMappingFiles. ReadMappingFiles " +
            "must be the full path to a read mapping file or a comma separated list of file paths.")
        .required()
        .validate { require(File(it).exists()) {"$it is not a valid file"} }

    val hvcfDir by option(help = "The directory containing the hvcf files used to build a HaplotypeGraph for path finding.")
        .required()
        .validate { require(File(it).isDirectory) {"$it is not a valid directory."} }

    val referenceGenome by option(help = "path to reference genome (fasta or fastq")
        .required()
        .validate { require(File(it).exists()) {"$it is not a valid file"} }

    val outputDir by option(help = "The directory where the output hvcfs will be written. The output file names will be <sampleName>.h.vcf.")
        .required()
        .validate { require(File(it).isDirectory) {"$it is not a valid directory."} }

    val probCorrect by option(help = "The probability that a mapped read was mapped correctly")
        .double()
        .default(0.99)
        .validate { require(it in 0.5..1.0) {"prob-correct must be between 0.5 and 1.0"} }

    val probSameGamete by option(help = "The probability of transitioning to the same gamete in the next reference range")
        .double()
        .default(0.99)

    val minGametes by option(help = "The minimum number of gametes with a haplotype in a reference range. " +
            "Reference ranges with fewer gametes will not be imputed.")
        .int()
        .default(2)
        .validate { require(it > -1) {"min-gametes must be a positive integer"} }

    val minReads by option(help = "The minimum number of reads per ReferenceRange. Reference ranges with fewer reads will not be imputed")
        .int()
        .default(0)
        .validate { require(it > -1) {"min-reads must be a positive integer."} }

    val maxReadsPerKb by option(help = "ReferenceRanges with more than max-reads-per-kb will not be imputed.")
        .int()
        .default(1000)
        .validate { require(it > -1) {"max-reads-per-kb must be a positive integer."} }

    val useLikelyAncestors by option(help="Use only the most likely ancestors of each sample for path finding.")
        .boolean()
        .default(false)

    val maxAncestors by option(help = "If use-likely-ancestors = true, use at most max-ancestors.")
        .int()
        .default(Int.MAX_VALUE)

    val minCoverage by option(help = "If use-likely-ancestors = true, use the fewest number of ancestors that together have this proportion of mappable reads.")
        .double()
        .default(1.0)
        .validate { require(it in 0.5..1.0) {"min-coverage must be between 0.5 and 1.0"} }

    val likelyAncestorFile by option(help="If useLikelyAncestors is true, a record of the ancestors used for each sample will be written to this file.")
        .default("")

    val threads by option(help = "number of threads used to find paths.").int().default(3)


    private val myLogger = LogManager.getLogger(HaploidPathFinding::class.java)

    private fun buildHaplotypeGraph(): HaplotypeGraph {
        val listOfHvcfFilenames = File(hvcfDir).listFiles().filter { it.name.endsWith(".h.vcf") || it.name.endsWith(".h.vcf.gz")}.map { it.path }
        return HaplotypeGraph(listOfHvcfFilenames)
    }

    override fun run() {
        //TODO determine how samples with no haplotype are being handled
        processKeyFile()
    }

    private fun processKeyFile() {
        //processes the key file to produce a map of sample name -> list of read mapping files for that sample
        myLogger.info("Processing key file $pathKeyfile")
        val keyFileLines = getBufferedReader(pathKeyfile).use { it.readLines()}
            .map { it.split("\t") }
            .filter { it.size > 1 }
        val headerList = keyFileLines[0]

        val sampleNameIsFirstColumn = headerList[0].equals("SampleName", true)
        val readMappingIsSecondColumn =  headerList[1].equals("ReadMappingFiles", true)
        require(sampleNameIsFirstColumn && readMappingIsSecondColumn) {"The first column heading of the keyfile " +
                "must be SampleName and the second column heading must be ReadMappingFiles."}

        //check for duplicate sample names (not allowed)
        val sampleNameSet = keyFileLines.drop(1).map { it[0] }.toSet()
        check(sampleNameSet.size == keyFileLines.size - 1) {"key file contains duplicate sample names"}

        //map of sampleName -> list of read mapping file names
        //Todo check for existence of output file for each sample and skip if it exists
        val sampleToFiles = keyFileLines.drop(1).associate { Pair(it[0], it[1]) }
            .mapValues { it.value.split(",").map{filename -> filename.trim()} }
        processReadMappings(sampleToFiles)
    }

    data class ReadMappingResult(val name: String, val readMappingCounts: Map<List<String>, Int>) //placeholder
    data class Path(val name: String, val hapidList: List<String>, val graph: HaplotypeGraph, val parentsUsed: List<MostLikelyParents.ParentStats> )

    private fun processReadMappings(sampleToFiles: Map<String, List<String>>) = runBlocking {
        myLogger.info("processing read mappings.")
        val myGraph = buildHaplotypeGraph()

        //create a channel for read mappings
        val readMappingChannel = Channel<ReadMappingResult>(10)

        //create a channel for paths
        val pathChannel = Channel<Path>(10)

        //create worker threads to process entries from the read mapping channel
        val jobList = mutableListOf<Job>()
        repeat(threads) {
            jobList.add(launch(Dispatchers.Default) {
                imputePath(myGraph, readMappingChannel, pathChannel)
            })
        }

        //create a coroutine to store paths
        launch(Dispatchers.IO) {
            savePath(pathChannel)
        }

        //load read mappings for each sample into a channel
        for (sampleFileList in sampleToFiles) {
            val listOfReadMaps = sampleFileList.value.map { filename -> importReadMapping(filename) }
            val readMappingsForSample = mergeReadMappings(listOfReadMaps)
            myLogger.info("submitting read mapping for $sampleFileList")
            readMappingChannel.send(ReadMappingResult(sampleFileList.key, readMappingsForSample))
        }
        readMappingChannel.close()

        //now wait for all the impute jobs to finish, then close the pathChannel
        for (imputeJob in jobList) imputeJob.join()
        pathChannel.close()

        //block should not finish until savePath is finished processing results
    }

    private suspend fun imputePath(graph: HaplotypeGraph, readMappings: ReceiveChannel<ReadMappingResult>, resultChannel: SendChannel<Path>) {
        val pathFinder = PathFinderWithViterbiHMM(graph = graph,
            probCorrect = probCorrect,
            sameGameteProbability = probSameGamete,
            minGametesPerRange = minGametes,
            minReadsPerRange = minReads,
            maxReadsPerKB =  maxReadsPerKb,
            useLikelyParents = useLikelyAncestors,
            maxParents = maxAncestors,
            minCoverage = minCoverage,
            isHaploidPath = true
            )

        for (result in readMappings) {
            val mappingsByRefrange = readMappingByRange(result.readMappingCounts, graph)
//            val hapidList = pathFinder.findBestHaploidPath(mappingsByRefrange)
            val hapidList = Pair(listOf<String>(), listOf<MostLikelyParents.ParentStats>())
            resultChannel.send(Path(result.name, hapidList.first, graph, hapidList.second))
        }
    }

    private suspend fun savePath(pathChannel : ReceiveChannel<Path>) {
        for (path in pathChannel) {
            if (path.hapidList.isNullOrEmpty())
                myLogger.error("No hVCF was written for ${path.name} because no haplotypes were imputed.")
            else {
                writeHvcf(path)
                if (path.parentsUsed.size > 0) appendParentStats(path)
            }
        }
    }

    private fun writeHvcf(myPath: Path) {

        //use method exportVariantContext() in VariantLoadingUtils
        //use the alt header lines from the graph filtered to the hapids in the VariantContexts
        //create a list of VariantContexts from myPath: use method createHvcfRecord() in VariantLoadingUtils


        //start the altHeaderList with the Reference haplotypes
        val referenceSequence = NucSeqIO(referenceGenome).readAll()
        val altHeadersSample = mutableListOf<AltHeaderMetaData>()
        for (hapid in myPath.hapidList) {
            val sampleAlt = myPath.graph.altHeader(hapid)
            check(sampleAlt != null) {"There is no ReferenceRange for sample haplotype $hapid."}
            altHeadersSample.add(sampleAlt)
        }

        val hapidToRefRange = myPath.graph.hapIdToRefRangeMap()
        val variantContexts = altHeadersSample.map { sampleAlt ->
            val refRange = hapidToRefRange[sampleAlt.id]!!
            val startPos = Position(refRange.contig, refRange.start)
            val endPos = Position(refRange.contig, refRange.end)

            //the ref allele needs to b the nucleotide at start Pos, not the haplotype id (hash)
            //Biokotlin positions are 0-based, so use refRange.start - 1
            val refAllele = referenceSequence[refRange.contig]!!.sequence[refRange.start - 1].name
            createHVCFRecord(myPath.name, startPos ,endPos, Pair(refAllele, sampleAlt.id))
        }

        //exportVariantContext()
        val hvcfFileName = if (outputDir.endsWith("/")) "${outputDir}${myPath.name}.h.vcf"
        else "${outputDir}/${myPath.name}.h.vcf"
        val headerSet = altHeadersSample.map { altHeaderMetadataToVCFHeaderLine(it) }.toSet()

        exportVariantContext(myPath.name, variantContexts, hvcfFileName, referenceSequence, headerSet)

    }

    private fun appendParentStats(path: Path) {
        val fileExists = File(likelyAncestorFile).exists()
        getBufferedWriter(likelyAncestorFile, append = true).use { myWriter ->
            if (!fileExists) myWriter.write("sample\tancestor\tgameteId\treads\tcoverage\n")
            for (stats in path.parentsUsed) {
                myWriter.write("${path.name}\t${stats.parent.name}\t${stats.parent.gameteId}\t${stats.readCount}\t${stats.coverage}\n")
            }
        }
    }

    companion object {

        /**
         * Takes read mapping counts, [readMappings], identifies the ReferenceRange for each hapid set,
         * then outputs the result as a map of ReferenceRange to read mapping counts for that ReferenceRange.
         * Here read mapping counts is a map of hapid set (List<String>) to a count.
         */
        fun readMappingByRange(readMappings: Map<List<String>, Int>, graph: HaplotypeGraph): Map<ReferenceRange, Map<List<String>, Int>> {
            val hapidToRefrangeMap = graph.ranges().map { refrange ->
                graph.hapIdToSampleGametes(refrange).keys.map { Pair(it, refrange) }
            }.flatten().toMap()

            val mappingByRange = mutableMapOf<ReferenceRange, MutableMap<List<String>, Int>>()
            for ((haplist, count) in readMappings.entries) {
                val refrange = hapidToRefrangeMap[haplist.first()]!!
                mappingByRange.getOrPut(refrange) { mutableMapOf() }[haplist] = count
            }
            return mappingByRange
        }
        fun altHeaderMetadataToVCFHeaderLine(altHeaderData: AltHeaderMetaData): VCFAltHeaderLine {

            return VCFAltHeaderLine(
                "<ID=${altHeaderData.id}, " +
                        "Description=\"haplotype data for line: ${altHeaderData.sampleName}\">," +
                        "Source=\"${altHeaderData.source}\",SampleName=\"${altHeaderData.sampleName}\"," +
                        "Regions=\"${altHeaderData.regions.joinToString(",") { "${it.first.contig}:${it.first.position}-${it.second.position}" }}\"," +
                        "Checksum=\"Md5\",RefRange=\"${altHeaderData.refRange}\">",
                VCFHeaderVersion.VCF4_2
            )
        }

    }
}

