package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.convert
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.pathing.AlignmentUtils
import net.maizegenetics.phgv2.pathing.KeyFileData
import net.maizegenetics.phgv2.pathing.ReadInputFile
import org.apache.logging.log4j.LogManager
import java.io.BufferedInputStream
import java.io.BufferedReader
import java.io.InputStreamReader

data class MEM(val readName: String, val readStart: Int, val readEnd: Int, val numHits: Int, val listMemHits: List<MEMHit>)
data class MEMHit(val contig: String, val strand: String, val pos: Int)


class MapReads : CliktCommand(help="Map reads to a pangenome using ropeBWT3") {

    private val myLogger = LogManager.getLogger(MapReads::class.java)

    val index by option(help = "The full path of the ropebwt3 index file.")
        .required()

    val readInputFiles: ReadInputFile by mutuallyExclusiveOptions<ReadInputFile>(
        option("--key-file", help = "Name of tab-delimited key file.  Columns for samplename and filename are required.  If using paired end fastqs, a filename2 column can be included. A value must be entered for either --key-file or --read-files.").convert{ ReadInputFile.KeyFile(it) },
        option("--read-files", help = "Comma separated list of fastq files for a single sample.  Either 1(for single end) or 2(for paired end) files can be input at a time this way.  Any more and an error will be thrown.").convert{ ReadInputFile.ReadFiles(it) }
    ).single().required()

    val outputDir by option("-o", "--output-dir", help = "Name for output ReadMapping file Directory (Required)")
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


    override fun run() {
        logCommand(this)
        myLogger.info("Mapping reads to pangenome")

        mapAllReadFiles(index, readInputFiles.getReadFiles(), outputDir, threads, minMemLength, maxNumHits, condaEnvPrefix)
    }

    fun mapAllReadFiles(index: String, keyFileDataEntries: List<KeyFileData>, outputDir: String, threads: Int, minMemLength: Int, maxNumHits: Int, condaEnvPrefix: String) {
        //time ../ropebwt3/ropebwt3 mem -t40 -l148 -p50 /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmd /workdir/zrm22/phgv2/ropeBWT/Reads/B97_HVMFTCCXX_L7_1.clean.fq.gz > B97_1_fullASM_pos_matches2NM.bed

        for(readFile in keyFileDataEntries) {
            val outputFile1 = "$outputDir/${readFile.sampleName}_1_readMapping.txt"
            mapSingleReadFile(index, readFile.sampleName, readFile.file1,outputFile1, threads, minMemLength, maxNumHits, condaEnvPrefix)
            if(readFile.file2 != "") {
                val outputFile2 = "$outputDir/${readFile.sampleName}_2_readMapping.txt"
                mapSingleReadFile(index, readFile.sampleName, readFile.file2, outputFile2, threads, minMemLength, maxNumHits, condaEnvPrefix)
            }
        }
    }

    fun mapSingleReadFile(index: String, sampleName: String,readFile: String, outputFile: String, threads: Int, minMemLength: Int, maxNumHits: Int, condaEnvPrefix: String) {
        myLogger.info("Mapping reads in $readFile to $index")

        val bedFileReader = setupMappingProcess(index, readFile, threads, minMemLength, maxNumHits, condaEnvPrefix)

        var currentLine = bedFileReader.readLine()
        val tempMems = mutableListOf<MEM>()
        val readMapping = mutableMapOf<List<String>, Int>()
        while(currentLine != null) {
            val alignmentParsed = parseStringIntoMem(currentLine)
            if(tempMems.isNotEmpty() && tempMems[0].readName != alignmentParsed.readName) {
                //write out the tempMems
                processMemsForRead(tempMems, readMapping, maxNumHits)
                tempMems.clear()
            }
            tempMems.add(alignmentParsed)
            currentLine = bedFileReader.readLine()
        }

        processMemsForRead(tempMems, readMapping, maxNumHits)

        bedFileReader.close()

        myLogger.info("Writing read mapping to $outputFile")
        AlignmentUtils.exportReadMapping(outputFile,readMapping, sampleName,Pair(readFile,""))
    }

    /**
     * Function to setup the ropebwt3 mem process
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

    fun processMemsForRead(tempMems: MutableList<MEM>, readMapping: MutableMap<List<String>, Int>, maxNumHits: Int) {
        //get the longest hits
        val maxLength = tempMems.maxOf { it.readEnd - it.readStart }
        //remove any hits that are not the longest
        val bestHits = tempMems.filter { it.readEnd - it.readStart == maxLength }

        val totalNumHits = bestHits.sumOf { it.numHits }

        if(totalNumHits <= maxNumHits) {
            //if the total number of hits is less than the maxNumHits, add all the hits to the readMapping
            //First turn the hapIds into a set then back to a list and sort so it will be consistent
            val hapIdsHit = bestHits.flatMap { it.listMemHits.map{hits -> hits.contig} }.toSet().toList().sorted()
            readMapping[hapIdsHit] = (readMapping[hapIdsHit]?:0) + 1
        }
    }
}