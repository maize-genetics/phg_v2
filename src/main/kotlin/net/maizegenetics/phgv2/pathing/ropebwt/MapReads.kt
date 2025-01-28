package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.convert
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.pathing.AlignmentUtils
import net.maizegenetics.phgv2.pathing.KeyFileData
import net.maizegenetics.phgv2.pathing.ReadInputFile
import org.apache.logging.log4j.LogManager
import java.io.BufferedInputStream
import java.io.BufferedReader
import java.io.File
import java.io.InputStreamReader

/**
 * Some classes to hold the data from the ropebwt3 mem output
 */
data class MEM(val readName: String, val readStart: Int, val readEnd: Int, val numHits: Int, val listMemHits: List<MEMHit>)
data class MEMHit(val contig: String, val strand: String, val pos: Int)

/**
 * MapReads will map read files independently to a pangenome indexed by ropeBWT3
 * This will create a standard read mapping file that can be used for path finding
 * Note that this will take some time to run.  Our internal tests show that WGS files can take about 20 minutes to run.
 * Each file from a pair are processed independently.
 */
class MapReads : CliktCommand(help="BETA: Map reads to a pangenome using ropeBWT3") {

    private val myLogger = LogManager.getLogger(MapReads::class.java)

    val index by option(help = "The full path of the ropebwt3 index file.")
        .required()

    val readInputFiles: ReadInputFile by mutuallyExclusiveOptions<ReadInputFile>(
        option("--key-file", help = "Name of tab-delimited key file.  Columns for samplename and filename are required.  If using paired end fastqs, a filename2 column can be included. A value must be entered for either --key-file or --read-files.").convert{ ReadInputFile.KeyFile(it) },
        option("--read-files", help = "Comma separated list of fastq files for a single sample.  Either 1(for single end) or 2(for paired end) files can be input at a time this way.  Any more and an error will be thrown.").convert{ ReadInputFile.ReadFiles(it) }
    ).single().required()

    val outputDir by option("-o", "--output-dir", help = "Name for output ReadMapping file Directory (Required)")
        .required()

    val hvcfDir by option(help = "Directory for the haplotype VCF files. These hvcfs will be used to filter the hits to only those that hit a single referenceRange.")
        .required()

    val threads by option(help = "Number of threads to use.")
        .int()
        .default(5)

    val minMemLength by option(help = "Minimum length of a match to be considered a match.")
        .int()
        .default(148)

    val maxNumHits by option(help = "Number of hits to report.  Note ropebwt can hit more than --max-num-hits but any alignment hitting more haplotypes than this will be ignored.")
        .int()
        .default(50)

    val condaEnvPrefix by option (help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
        .default("")

    val maxStart by option(help = "Maximum start position for a read to be considered a match. Any alignments with a start above this will be ignored.")
        .int()
        .default(0)

    val minEnd by option(help = "Minimum end position for a read to be considered a match. Any alignments with an end below this will be ignored.")
        .int()
        .default(70)


    override fun run() {
        logCommand(this)

        myLogger.info("Building the Graph")

        val graph = HaplotypeGraph(hvcfDir)

        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()

        myLogger.info("Mapping reads to pangenome")

        mapAllReadFiles(index, readInputFiles.getReadFiles(), outputDir, threads, minMemLength, maxNumHits, condaEnvPrefix, hapIdToRefRangeMap, maxStart, minEnd)
    }

    /**
     * Function to map all the read files in the keyFileDataEntries to the index
     */
    fun mapAllReadFiles(index: String, keyFileDataEntries: List<KeyFileData>, outputDir: String, threads: Int,
                        minMemLength: Int, maxNumHits: Int, condaEnvPrefix: String,
                        hapIdToRefRangeMap : Map<String, List<ReferenceRange>>, maxStart: Int, minEnd: Int) {
        //Loop through the keyFileDataEntries and map the reads
        //If there is a second file it processes that and makes a separate readMapping file.
        val readNameToFileMap = mutableMapOf<String,MutableList<String>>()
        for(readFile in keyFileDataEntries) {
            val fileList = readNameToFileMap[readFile.sampleName]?: mutableListOf()
            val outputFile1 = "$outputDir/${readFile.sampleName}_1_readMapping.txt"
            mapSingleReadFile(index, readFile.sampleName, readFile.file1,outputFile1, threads, minMemLength, maxNumHits, condaEnvPrefix, hapIdToRefRangeMap, maxStart, minEnd)
            fileList.add(outputFile1)
            if(readFile.file2 != "") {
                val outputFile2 = "$outputDir/${readFile.sampleName}_2_readMapping.txt"
                mapSingleReadFile(index, readFile.sampleName, readFile.file2, outputFile2, threads, minMemLength, maxNumHits, condaEnvPrefix, hapIdToRefRangeMap, maxStart, minEnd)
                fileList.add(outputFile2)
            }
            readNameToFileMap[readFile.sampleName] = fileList
        }
        exportPathKeyFile(outputDir, readNameToFileMap)
    }

    /**
     * Function to map a single read file to the index and write the read mapping to the outputFile
     */
    fun mapSingleReadFile(index: String, sampleName: String,readFile: String, outputFile: String, threads: Int,
                          minMemLength: Int, maxNumHits: Int, condaEnvPrefix: String,
                          hapIdToRefRangeMap: Map<String,List<ReferenceRange>>, maxStart: Int, minEnd: Int) {
        myLogger.info("Mapping reads in $readFile to $index")

        val bedFileReader = setupMappingProcess(index, readFile, threads, minMemLength, maxNumHits, condaEnvPrefix)

        val readMapping = createReadMappingsForFileReader(bedFileReader, maxNumHits, hapIdToRefRangeMap, maxStart, minEnd)

        bedFileReader.close()

        myLogger.info("Writing read mapping to $outputFile")
        AlignmentUtils.exportReadMapping(outputFile,readMapping, sampleName,Pair(readFile,""))
    }



    /**
     * Function to setup the ropebwt3 mem process and pass a BufferedReader for use by the rest of the program
     * //time ../ropebwt3/ropebwt3 mem -t40 -l148 -p50 /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmd /workdir/zrm22/phgv2/ropeBWT/Reads/B97_HVMFTCCXX_L7_1.clean.fq.gz > B97_1_fullASM_pos_matches2NM.bed
     *
     */
    fun setupMappingProcess(index: String, readFile: String, threads: Int, minMemLength: Int, maxNumHits: Int, condaEnvPrefix: String): BufferedReader {
        val prefixArg = if(condaEnvPrefix.isNotBlank()) {
            Pair("-p",condaEnvPrefix)
        }
        else {
            Pair("-n", "phgv2-ropebwt-conda")
        }

        val ropebwt3Process = ProcessBuilder("conda","run",prefixArg.first,prefixArg.second,"ropebwt3", "mem", "-t$threads", "-l$minMemLength", "-p$maxNumHits", index, readFile)
            .redirectError(ProcessBuilder.Redirect.INHERIT)
            .start()

        return BufferedReader(InputStreamReader(ropebwt3Process.inputStream))
    }

    /**
     * Function to create a readMapping Map from a BufferedReader of a ropebwt3 mem output
     * This will work directly from a file but this command is not setup to do that yet.
     */
    fun createReadMappingsForFileReader(
        bedFileReader: BufferedReader,
        maxNumHits: Int,
        hapIdToRefRangeMap: Map<String, List<ReferenceRange>>,
        maxStart: Int,
        minEnd: Int
    ): MutableMap<List<String>, Int> {
        var currentLine = bedFileReader.readLine()
        val tempMems = mutableListOf<MEM>()
        val readMapping = mutableMapOf<List<String>, Int>()
        while (currentLine != null) {
            if(currentLine.isEmpty()) {
                currentLine = bedFileReader.readLine()
                continue
            }
            val alignmentParsed = parseStringIntoMem(currentLine)
            if (tempMems.isNotEmpty() && tempMems[0].readName != alignmentParsed.readName) {
                //write out the tempMems
                processMemsForRead(tempMems, readMapping, maxNumHits,hapIdToRefRangeMap, maxStart, minEnd)
                tempMems.clear()
            }
            tempMems.add(alignmentParsed)
            currentLine = bedFileReader.readLine()
        }

        processMemsForRead(tempMems, readMapping, maxNumHits, hapIdToRefRangeMap, maxStart, minEnd)
        return readMapping
    }

    /**
     * Function to parse the current alignment line from ropebwt3 mem into a usable object
     */
    fun parseStringIntoMem(string: String) : MEM {
        val split = string.split("\t")
        val readName = split[0]
        val readStart = split[1].toInt()
        val readEnd = split[2].toInt()
        val numHits = split[3].toInt()
        val listMemHits = split.subList(5, split.size).map {
            val hitSplit = it.split(":")
            MEMHit(hitSplit[0], hitSplit[1], hitSplit[2].toInt())
        }
        return MEM(readName, readStart, readEnd, numHits, listMemHits)
    }

    /**
     * Function to process the mems for a single read and add them to the readMapping
     * This will filter things based on the maxStart and minEnd positions, then  maxNumHits and retain the mems that are longest
     * Because MEMs are Maximal Exact Matches, the longest MEMs are the best hits
     */
    fun processMemsForRead(tempMems: List<MEM>, readMapping: MutableMap<List<String>, Int>, maxNumHits: Int,
                           hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, maxStart: Int, minEnd: Int) {
        val posFilteredMems = tempMems.filter { it.readStart <= maxStart && it.readEnd >= minEnd }
        if(posFilteredMems.isEmpty()) {
            return
        }
        //get the longest hits
        val maxLength = posFilteredMems.maxOf { it.readEnd - it.readStart }
        //remove any hits that are not the longest
        val bestHits = posFilteredMems.filter { it.readEnd - it.readStart == maxLength }

        val totalNumHits = bestHits.sumOf { it.numHits }

        if(totalNumHits <= maxNumHits) {
            //if the total number of hits is less than the maxNumHits, Filter the haps down to a single ref range and then add all those to the readMapping
            //First turn the hapIds into a set then back to a list and sort so it will be consistent
            val bestHapIds = bestHits.flatMap { it.listMemHits.map{hits -> hits.contig} }.toSet()

            val filteredBestHapIds = filterToOneReferenceRange(bestHapIds, hapIdToRefRangeMap)

            val hapIdsHit = filteredBestHapIds.sorted()
            readMapping[hapIdsHit] = (readMapping[hapIdsHit]?:0) + 1
        }
    }

    /**
     * Function to filter the hapIds to only those that hit a single reference range
     * This will pick the highest occurring reference range and if there is a tie it will choose the first one.
     */
    fun filterToOneReferenceRange(hapIds: Set<String>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>) : List<String> {
        //figure out which reference range is hit the most
        val refRangeCounts = hapIds.map { hapIdToRefRangeMap[it] }
            .filterNotNull()
            .flatten()
            .groupingBy { it }
            .eachCount()

        val maxRefRange = refRangeCounts.maxByOrNull { it.value }?.key
        //Then filter out any hapIds that don't have that reference range
        return hapIds.filter { hapIdToRefRangeMap[it]?.contains(maxRefRange)?:false }
    }

    /**
     * Function to export the path key file for the readMapping files that have been processed.
     */
    fun exportPathKeyFile(outputDir: String, readNameToFileMap: Map<String, List<String>>) {
        val pathKeyFile = "$outputDir/pathKeyFile.txt"
        myLogger.info("Writing pathKeyFile to $pathKeyFile")
        bufferedWriter(pathKeyFile).use {writer ->
            writer.write("sampleName\tfilename\n")
            for ((sampleName, fileList) in readNameToFileMap) {
                for (file in fileList) {
                    writer.write("$sampleName\t$file\n")
                }
            }
            writer.close()
        }
    }
}