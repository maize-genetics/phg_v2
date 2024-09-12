package net.maizegenetics.phgv2.pathing

import biokotlin.seqIO.NucSeqIO
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.boolean
import com.github.ajalt.clikt.parameters.types.choice
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.int
import htsjdk.variant.variantcontext.VariantContext
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.Job
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.channels.ReceiveChannel
import kotlinx.coroutines.channels.SendChannel
import kotlinx.coroutines.launch
import kotlinx.coroutines.runBlocking
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.pathing.AlignmentUtils.Companion.importReadMapping
import net.maizegenetics.phgv2.pathing.AlignmentUtils.Companion.mergeReadMappings
import net.maizegenetics.phgv2.utils.*
import org.apache.logging.log4j.LogManager
import java.io.File
import java.nio.file.Files
import kotlin.io.path.absolutePathString

sealed class PathInputFile {
    abstract fun getReadFiles() : List<KeyFileData>
    data class KeyFile(val keyFile: String): PathInputFile() {
        @Override
        override fun getReadFiles(): List<KeyFileData> {
            check(File(keyFile).exists()) { "Key file $keyFile does not exist." }
            val linesWithHeader = File(keyFile).bufferedReader().readLines()
            val header = linesWithHeader.first().split("\t")
            //convert the header into a map of column name to column index
            val headerMap = header.mapIndexed { index, s -> s to index }.toMap()
            check(headerMap.containsKey("sampleName")) { "Key file $keyFile must have a column named sampleName." }
            check(headerMap.containsKey("filename")) { "Key file $keyFile must have a column named filename." }
            return linesWithHeader.drop(1).map{lines -> lines.split("\t")}.map { linesSplit ->
                KeyFileData(linesSplit[headerMap["sampleName"]!!], linesSplit[headerMap["filename"]!!], if(headerMap.containsKey("filename2") && linesSplit.indices.contains(headerMap["filename2"]!!)) linesSplit[headerMap["filename2"]!!] else "")
            }
        }
    }

    data class ReadFiles(val readFiles: String): PathInputFile() {
        @Override
        override fun getReadFiles(): List<KeyFileData> {
            check(readFiles.isNotEmpty()) { "--read-files must have at least one file." }
            val isReadMapping = readFiles.contains("_readMapping.txt")
            return if (isReadMapping) {
                val sampleName = File(readFiles.split(",")[0]).name.removeSuffix("_readMapping.txt")
                listOf(KeyFileData(sampleName, readFiles, ""))
            } else {
                val fileNames = readFiles.split(",")
                check(fileNames.size <= 2) { "--read-files must have 1 or 2 files separated by commas.  You provided: ${fileNames.size}" }
                val sampleName = File(fileNames[0]).name.removeSuffix(".gz").removeSuffix(".fq").removeSuffix(".fastq")
                listOf(KeyFileData(sampleName,fileNames.first(), if(fileNames.size==1) "" else fileNames.last()))
            }
        }

    }
}


/**
 * This version of FindPaths uses .h.vcf files to build the HaplotypeGraph used for path finding
 * and writes the imputed paths as .h.vcf files to an output directory. It does not use a tiledb for anything.
 * That is expected to change at some point.
 *
 * The keyfile can list read mapping files or fastq files and can take a comma separated list of files to
 * combine for path finding. For paired-end reads filename and filename2 can also be comma separated lists.
 *
 * --read-files expects either a comma-separated list of read mapping files, one fastq file for single end reads,
 * or a pair of fastq files for paired end reads.
 */
class FindPaths: CliktCommand(help = "Impute best path(s) using read mappings.")  {

    val readInputFiles: PathInputFile by mutuallyExclusiveOptions<PathInputFile>(
        option("--path-keyfile", help = "Name of tab-delimited key file.  Columns for samplename and filename" +
                " are required. Files must be either read mapping files (ending in _readMapping.txt) or fastq files. " +
                "If using paired end fastqs, a filename2 column can be included. A value must be entered for " +
                "either --path-keyfile or --read-files.")
            .convert{ PathInputFile.KeyFile(it) },
        option("--read-files", help = "Comma separated list of fastq files for a single sample. " +
                "Either 1(for single end) or 2(for paired end) files can be input at a time this way.  Any more and " +
                "an error will be thrown. If listing files from read mapping, they must end in _readMapping.txt.")
            .convert{ PathInputFile.ReadFiles(it) }
    ).single().required()

    val hvcfDir by option(help = "The directory containing the hvcf files used to build a HaplotypeGraph for path finding. Required parameter.")
        .required()
        .validate { require(File(it).isDirectory) {"$it is not a valid directory. Required parameter."} }

    val referenceGenome by option(help = "path to reference genome (fasta or fastq). Required parameter.")
        .required()
        .validate { require(File(it).exists()) {"$it is not a valid file"} }

    val outputDir by option(help = "The directory where the output hvcfs will be written. The output file names will be <sampleName>.h.vcf. Required parameter.")
        .required()
        .validate { require(File(it).isDirectory) {"$it is not a valid directory."} }

    val outParentsDir by option(help="The directory where the imputed parents (ancestors) will be written for each sample. File names will be <sampleName>_imputed_parents.txt. If no directory name is supplied, the files will not be written.")
        .default("")

    val pathType by option(help = "The type of path to find. Must be lower case 'haploid' or 'diploid' (without quotes). 'haploid' infers a single path through the graph. 'diploid' infers a pair of paths. Required parameter.")
        .choice("haploid", "diploid")
        .required()

    val kmerIndex by option(help = "The name and path to the kmerIndex. Default = <hvcfDir>/kmerIndex.txt")
        .default("")

    val probCorrect by option(help = "The probability that a mapped read was mapped correctly.")
        .double()
        .default(0.99)
        .validate { require(it in 0.5..1.0) {"prob-correct must be between 0.5 and 1.0"} }

    val probSameGamete by option(help = "The probability of transitioning to the same gamete (sample) in the next reference range.")
        .double()
        .default(0.99)
        .validate { require(it in 0.01..1.0) {"prob-correct must be between 0.01 and 1.0"} }

    val minGametes by option(help = "The minimum number of gametes with a haplotype in a reference range. " +
            "Reference ranges with fewer gametes will not be imputed.")
        .int()
        .default(1)
        .validate { require(it > -1) {"min-gametes must be a positive integer"} }

    val minReads by option(help = "The minimum number of reads per ReferenceRange. Reference ranges with fewer reads " +
            "will not be imputed. If minReads = 0, all ReferenceRanges will be imputed.")
        .int()
        .default(0)
        .validate { require(it > -1) {"min-reads must be a positive integer."} }

    val inbreedingCoefficient by option(help = "The estimated coefficient of inbreeding for the samples being evaluated. " +
            "Only used for diploid path type.")
        .double()
        .default(0.0)
        .validate { require(it in 0.0..1.0) {"inbreeding-coefficient must be between 0.0 and 1.0"} }

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

    private val myLogger = LogManager.getLogger(FindPaths::class.java)

    override fun run() {
        val keyFileLines = readInputFiles.getReadFiles()
        require(keyFileLines.isNotEmpty()) {"Must provide either --path-keyfile or --read-files."}

        //create the outParentsDir, if it does not already exist
        if(outParentsDir.isNotBlank()) File(outParentsDir).mkdirs()

        var samplesToReadMappingFiles = when (getKeyFileType(keyFileLines)) {
            KeyfileType.READ ->
                keyFileLines.associate { Pair(it.sampleName, it.file1.split(",")) }

            KeyfileType.FASTQ -> mapReads(keyFileLines)
        }
        samplesToReadMappingFiles = checkForExistingOutput(samplesToReadMappingFiles)

        processReadMappings(samplesToReadMappingFiles)
    }

    private enum class KeyfileType {FASTQ, READ}
    private fun getKeyFileType(keyfileLines: List<KeyFileData>): KeyfileType {
        val fastqEndings = listOf(".fastq", ".fastq.gz", ".fq", ".fq.gz")
        val readMappingEndings = listOf("_readMapping.txt", "_readMapping.txt.gz")
        val firstFilename = keyfileLines.first().file1.split(",")[0]
        return when {
            fastqEndings.any { firstFilename.endsWith(it) } -> KeyfileType.FASTQ
            readMappingEndings.any { firstFilename.endsWith(it) } -> KeyfileType.READ
            else -> throw IllegalArgumentException("Key filenames must be either fastq or readMapping. Fastq files" +
                    " must end in one of $fastqEndings. ReadMapping files must end in one of " +
                    " $readMappingEndings")
        }
    }

    private fun mapReads(keyfileLines: List<KeyFileData>): Map<String, List<String>> {

        val outputDirPath = Files.createTempDirectory("readMapping")

        //write a key file to outputDirPath
        val readKeyfile = outputDirPath.resolve("readKeyFile.txt").absolutePathString()
        getBufferedWriter(readKeyfile).use { myWriter ->
            myWriter.write("sampleName\tfilename\tfilename2\n")
            for (keyLine in keyfileLines) {
                val filenames = keyLine.file1.split(",")
                val filenames2 = keyLine.file2.split(",")
                require(filenames2[0] == "" || filenames.size == filenames2.size) {"Keyfile problem: The number of filenames in filename2 and filename are different for sampleName = ${keyLine.sampleName}"}
                if (filenames2[0] == "") {
                    for (ndx in filenames.indices) myWriter.write("${keyLine.sampleName}_$ndx\t${filenames[ndx]}\t\n")
                } else {
                    for (ndx in filenames.indices) myWriter.write("${keyLine.sampleName}_$ndx\t${filenames[ndx]}\t${filenames2[ndx]}\n")
                }
            }
        }

        val argList = mutableListOf("--hvcf-dir", hvcfDir, "--key-file", readKeyfile,
            "--output-dir", outputDirPath.absolutePathString())
        if (kmerIndex.isNotBlank()) {
            argList.add("--kmer-index")
            argList.add(kmerIndex)
        }

        MapKmers().main(argList)

        //return a map of sampleName -> a file named <sampleName>_readMapping.txt in the outputDir
        return keyfileLines.associate { Pair(it.sampleName,
            listOf(outputDirPath.resolve("${it.sampleName}_readMapping.txt").absolutePathString())) }

    }

    private fun buildHaplotypeGraph(): HaplotypeGraph {
        val listOfHvcfFilenames = File(hvcfDir).listFiles()
            .filter {  it.name.endsWith("h.vcf") || it.name.endsWith("h.vcf.gz")  }.map { it.path }

        return HaplotypeGraph(listOfHvcfFilenames)
    }

    /**
     * Takes a map of sample -> file list as input. If the output directory contains an hvcf file for any sample,
     * that sample will be removed from the map and the trimmed map will be returned.
     */
    private fun checkForExistingOutput(sampleToFiles: Map<String, List<String>>): Map<String, List<String>> {

        val existingSamples = File(outputDir).listFiles().filter { it.name.endsWith(".h.vcf") || it.name.endsWith(".h.vcf.gz") }
            .map { file -> file.name.substringBefore(".h.vcf") }

        val newSamples = sampleToFiles.filter { (sampleName, _) -> !existingSamples.contains(sampleName) }

        if (newSamples.size < sampleToFiles.size)
            myLogger.info("Processing has been requested for ${sampleToFiles.size} samples. ${newSamples.size} will be imputed." +
                    "The others already exist in the output directory")

        return newSamples
    }

    private data class ReadMappingResult(val name: String, val readMappingCounts: Map<List<String>, Int>)
    private data class Path(val name: String, val hapidList: List<PathFinderWithViterbiHMM.PathNode>, val graph: HaplotypeGraph, val likelyParents: List<MostLikelyParents.ParentStats>)

    /**
     * Takes a map of sample name -> list of read mapping files. It processes the samples in parallel.
     * Multi-threading is across samples. The results are written to files by an additional thread.
     * The results are a h.vcf for each sample and, if the most likely parents are used, a single file of
     * the chosen ancestors (parents) for each sample.
     */
    private fun processReadMappings(sampleToFiles: Map<String, List<String>>) = runBlocking {
        myLogger.info("processing read mappings.")
        val myGraph = buildHaplotypeGraph()

        //check conditions for useLikelyAncestor
        if (useLikelyAncestors) {
            require(maxAncestors < myGraph.sampleGametesInGraph().size || minCoverage < 1.0)
            {"UseLikelyAncestors is true but likelyAncestors will not be checked because minCoverage is 1.0 " +
                    "and maxAncestors is >= the number of potential ancestors."}
            if (likelyAncestorFile.isBlank()) myLogger.warn("Likely ancestors will be determined for each sample " +
                    "but information about the ancestors used will not be written to a file because no file name was provided.")
        }

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
        val saveJob = launch(Dispatchers.IO) {
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

        saveJob.join()
        //block does not finish until savePath is finished processing results
    }

    private suspend fun imputePath(graph: HaplotypeGraph, readMappingChannel: ReceiveChannel<ReadMappingResult>, pathChannel: SendChannel<Path>) {
        val pathFinder = PathFinderWithViterbiHMM(graph = graph,
            isHaploidPath = pathType == "haploid",
            probCorrect = probCorrect,
            sameGameteProbability = probSameGamete,
            minGametesPerRange = minGametes,
            minReadsPerRange = minReads,
            maxReadsPerKB =  maxReadsPerKb,
            useLikelyParents = useLikelyAncestors,
            maxParents = maxAncestors,
            minCoverage = minCoverage,
            inbreedCoef = inbreedingCoefficient
        )

        for (result in readMappingChannel) {
            val mappingsByRefrange = readMappingByRange(result.readMappingCounts, graph)
            val pathResult = pathFinder.findBestPath(mappingsByRefrange)
            pathChannel.send(Path(result.name, pathResult.first, graph, pathResult.second))
        }

    }

    private suspend fun savePath(pathChannel : ReceiveChannel<Path>) {
        for (path in pathChannel) {
            if (path.hapidList.isNotEmpty()) {
                writeHvcf(path)
                if (useLikelyAncestors && likelyAncestorFile.isNotBlank()) appendParentStats(path)
            } else {
                myLogger.warn("Path was not imputed for ${path.name}")
            }
        }
    }

    private fun writeHvcf(myPath: Path) {
        //use method exportVariantContext() in VariantLoadingUtils
        //use the alt header lines from the graph for the hapids in myPath
        //create a list of VariantContexts from myPath: use method createHvcfRecord() in VariantLoadingUtils

        val referenceSequence = NucSeqIO(referenceGenome).readAll()
        val altHeadersSample = mutableListOf<AltHeaderMetaData>()

        val variantContextList = mutableListOf<VariantContext>()

        //process the path
        for ( node in myPath.hapidList) {
            val startPos = Position(node.refRange.contig, node.refRange.start)
            val endPos = Position(node.refRange.contig, node.refRange.end)
            val refAllele = referenceSequence[node.refRange.contig]!!.sequence[node.refRange.start - 1].name
            val hapids = node.sampleGametes.map { myPath.graph.sampleToHapId(node.refRange, it) }

            variantContextList.add(createDiploidHVCFRecord(myPath.name, startPos, endPos, hapids, refAllele))

            //add alt headers for the unique hapids
            hapids.filterNotNull().distinct().mapNotNull { myPath.graph.altHeader(it) }.forEach { altHeadersSample.add(it) }
        }

        //exportVariantContext()
        val hvcfFileName = if (outputDir.endsWith("/")) "${outputDir}${myPath.name}.h.vcf"
        else "${outputDir}/${myPath.name}.h.vcf"
        val headerSet = altHeadersSample.map { altHeaderMetadataToVCFHeaderLine(it) }.toSet()

        exportVariantContext(myPath.name, variantContextList, hvcfFileName, referenceSequence, headerSet)

        //write the parents (if requested)
        if(outParentsDir.isNotBlank()) {
            val outFile = File(outParentsDir).resolve("${myPath.name}_imputed_parents.txt")
            getBufferedWriter(outFile).use { myWriter ->
                myWriter.write("chrom\tstart\tend\tsample1\tsample2\n")
                for ( node in myPath.hapidList) {
                    myWriter.write("${node.refRange.contig}\t${node.refRange.start}\t${node.refRange.end}" +
                            "\t${node.sampleGametes.getOrElse(0){""}}\t${node.sampleGametes.getOrElse(1){""}}\n")
                }
            }
        }

    }

    /**
     * Appends taxa chosen by likely parents to a file along with incremental read counts and cumulative coverage.
     */
    private fun appendParentStats(path: Path) {
        val fileExists = File(likelyAncestorFile).exists()
        getBufferedWriter(likelyAncestorFile, append = true).use { myWriter ->
            if (!fileExists) myWriter.write("sample\tancestor\tgameteId\treads\tcoverage\n")
            for (stats in path.likelyParents) {
                myWriter.write("${path.name}\t${stats.parent.name}\t${stats.parent.gameteId}\t${stats.readCount}\t${stats.coverage}\n")
            }
        }
    }

    private fun readMappingByRange(readCounts: Map<List<String>, Int>, graph: HaplotypeGraph): Map<ReferenceRange, Map<List<String>, Int>> {
        val hapid2Refrange = graph.hapIdToRefRangeMap()
        //groups readCounts for the entire genome into separate maps for each reference range
        //since some hapids map to more than one reference range, assign each hapid set to the reference range with the most hapids in the set

        return readCounts.entries.groupBy { referenceRangeForHapidList(it.key, hapid2Refrange) }
            .mapValues { (_,mapEntries) -> mapEntries.associateBy({it.key}, {it.value}) }
    }

    private fun referenceRangeForHapidList(hapidList: List<String>, hapid2Refrange: Map<String, List<ReferenceRange>>): ReferenceRange {
        val referenceRangeList = hapidList.mapNotNull { hapid2Refrange[it] }.flatten()
        val referenceRangeCount = referenceRangeList.groupingBy { it }.eachCount()
        return referenceRangeCount.maxBy {it.value}.key
    }

}