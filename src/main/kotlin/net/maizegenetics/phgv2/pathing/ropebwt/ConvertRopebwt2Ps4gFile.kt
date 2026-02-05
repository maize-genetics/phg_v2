package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.headerCommand
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.utils.Position
import org.apache.logging.log4j.LogManager
import java.io.File

sealed class BedInputFile {
    abstract fun getBedFiles(): List<String>
    data class BedFile(val bedFile: String) : BedInputFile() {
        @Override
        override fun getBedFiles(): List<String> {
            check(File(bedFile).exists()) { "Bed file $bedFile does not exist." }
            return listOf(bedFile)
        }
    }

    data class BedFileList(val bedFiles: String) : BedInputFile() {
        @Override
        override fun getBedFiles(): List<String> {
            check(bedFiles.isNotEmpty()) { "--bed-files must have at least one file." }
            val fileNames = bedFiles.split(",")
            fileNames.forEach {
                check(File(it).exists()) { "Bed file $it does not exist." }
            }
            return fileNames
        }
    }

    data class BedListFile(val bedListFile: String) : BedInputFile() {
        @Override
        override fun getBedFiles(): List<String> {
            check(File(bedListFile).exists()) { "Bed list file $bedListFile does not exist." }
            val lines = File(bedListFile).bufferedReader().readLines()
            lines.forEach {
                check(File(it).exists()) { "Bed file $it listed in $bedListFile does not exist." }
            }
            return lines
        }
    }
}


/**
 * This class will convert a RopebwtBed file to a PS4G file.  It will only work with ropebwt3 files where the reads are
 * aligned to the whole assembly chromosomes using the mem command.  MEMs are Maximal Exact Matches and are used to
 * determine what the optimal mapping is.  One downside to this approach is that if there is a SNP in the middle of the read, the mappings will be ignored.
 * We can run SW in the future to fix this but it is very slow.
 */
class ConvertRopebwt2Ps4gFile : CliktCommand(help = "Convert RopebwtBed to PS4G") {

    val myLogger = LogManager.getLogger(ConvertRopebwt2Ps4gFile::class.java)

    //val readInputFiles: ReadInputFile by mutuallyExclusiveOptions<ReadInputFile>(
    //        option("--key-file", help = "Name of tab-delimited key file.  Columns for samplename and filename are required.  If using paired end fastqs, a filename2 column can be included. A value must be entered for either --key-file or --read-files.").convert{ ReadInputFile.KeyFile(it) },
    //        option("--read-files", help = "Comma separated list of fastq files for a single sample.  Either 1(for single end) or 2(for paired end) files can be input at a time this way.  Any more and an error will be thrown.").convert{ ReadInputFile.ReadFiles(it) }
    //    ).single().required()

    val ropebwtBedFiles: BedInputFile by mutuallyExclusiveOptions<BedInputFile>(
            option("--ropebwt-bed", help = "RopebwtBed file.").convert{ BedInputFile.BedFile(it) },
            option("--ropebwt-bed-files", help = "Comma separated list of RopebwtBed files.").convert{ BedInputFile.BedFileList(it) },
            option("--ropebwt-bed-list-file", help = "File containing list of RopebwtBed files, one per line.").convert{ BedInputFile.BedListFile(it) }
        ).single().required()

//    val ropebwtBed by option(help = "RopebwtBed file")
//        .required()

    val outputDir by option(help = "Output directory or file name")
        .required()

   val splineKnotDir by option(help = "Spline Knot Directory")
        .required()

    val minMemLength by option(help = "Minimum length of a match to be considered a match.")
        .int()
        .default(148)

    val maxNumHits by option(help = "Number of hits to report.  Note ropebwt can hit more than --max-num-hits but any alignment hitting more haplotypes than this will be ignored.")
        .int()
        .default(50)

    val maxRange by option(help = "Maximum range (in bins) that a single alignment can hit on different haplotypes. Alignments that hit a wider range than this will be ignored.")
        .int()
        .default(50)
        .validate { require(it >= 0) {"max range must be non-negative"}  }

    val sortPositions by option(help = "Sort positions in the resulting PS4G file.")
        .flag(default = true)

    /**
     * Function to run the command.  This goes from a RopeBWT3 BED file for reads aligned against the whole chromosomes to a PS4G file.
     */
    override fun run() {
        logCommand(this)

        val command = headerCommand(this)

        myLogger.info("Convert Ropebwt to PS4G")

        myLogger.info("Loading Spline Knot File")
        val (splineKnots, chrIndexMap, gameteIndexMap) = SplineUtils.loadSplineKnotLookupFromDirectory(splineKnotDir)

        myLogger.info("Converting Spline Knots to Splines")

        val splineLookup = LinearLookupFunction(splineKnots, chrIndexMap)

        val sampleGameteIndexMap = gameteIndexMap.map { SampleGamete(it.key) to it.value}.toMap()

        val bedFilesToProcess = ropebwtBedFiles.getBedFiles()
        myLogger.info("Processing ${bedFilesToProcess.size} Ropebwt Bed Files")
        //check to make sure its not empty and the files exist
        check(bedFilesToProcess.isNotEmpty()) {"No Ropebwt Bed files provided"}
        bedFilesToProcess.forEach { ropebwtBed ->
            processSingleRopeBwtBed(ropebwtBed, outputDir,splineLookup, gameteIndexMap, sampleGameteIndexMap, command,minMemLength, maxNumHits, maxRange, sortPositions)
        }
    }

    private fun processSingleRopeBwtBed(
        ropebwtBed: String,
        outputDir: String,
        splineLookup: LinearLookupFunction,
        gameteIndexMap: Map<String, Int>,
        sampleGameteIndexMap: Map<SampleGamete, Int>,
        command: String,
        minMemLength: Int,
        maxNumHits: Int,
        maxRange: Int,
        sortPositions: Boolean = true
    ) {
        myLogger.info("Building PS4G Output File Name for file: $ropebwtBed")
        //build the output file name
        val outputFile = if (File(outputDir).isDirectory) {
            PS4GUtils.buildOutputFileName(ropebwtBed, outputDir)
        } else {
            outputDir
        }

        myLogger.info("Processing Ropebwt Bed File: $ropebwtBed")
        myLogger.info("Building PS4G Data")
        //build and write out the PS4G file
        val (ps4GData, sampleGameteCountMap) = buildPS4GData(
            ropebwtBed,
            splineLookup,
            gameteIndexMap,
            minMemLength,
            maxNumHits,
            maxRange,
            sortPositions
        )

        myLogger.info("Writing out PS4G File")
        PS4GUtils.writeOutPS4GFile(ps4GData, sampleGameteCountMap, sampleGameteIndexMap, outputFile, listOf(), command)
    }

    /**
     * Function to build the output PS4G file.
     */
    fun buildPS4GData(ropebwtBed: String,
                      splineLookup: LinearLookupFunction,
                      gameteToIdxMap: Map<String,Int>,
                      minMEMLength: Int, maxNumHits: Int,
                      maxRange: Int,
                      sortPositions: Boolean = true) : Pair<List<PS4GData>, Map<SampleGamete,Int>> {

        val gameteIdxToSampleGameteMap = gameteToIdxMap.map { it.value to SampleGamete(it.key) }.toMap()
        val bedFileReader = bufferedReader(ropebwtBed)
        var currentLine = bedFileReader.readLine()
        val tempMems = mutableListOf<MEM>()
        val countMap = mutableMapOf<Pair<Position, List<Int>>, Int>()
        val sampleGameteCountMap = mutableMapOf<SampleGamete,Int>()
        while (currentLine != null) {
            if(currentLine.isEmpty()) {
                currentLine = bedFileReader.readLine()
                continue
            }
            val alignmentParsed = RopeBWTUtils.parseStringIntoMem(currentLine)
            if (tempMems.isNotEmpty() && tempMems[0].readName != alignmentParsed.readName) {
                //write out the tempMems
                processTempMEMs(
                    tempMems,
                    splineLookup,
                    minMEMLength,
                    maxNumHits,
                    maxRange,
                    gameteToIdxMap,
                    countMap,
                    sampleGameteCountMap,
                    gameteIdxToSampleGameteMap
                )
                tempMems.clear()
            }
            tempMems.add(alignmentParsed)
            currentLine = bedFileReader.readLine()
        }

        processTempMEMs(
            tempMems,
            splineLookup,
            minMEMLength,
            maxNumHits,
            maxRange,
            gameteToIdxMap,
            countMap,
            sampleGameteCountMap,
            gameteIdxToSampleGameteMap
        )

        val ps4gDataList = PS4GUtils.convertCountMapToPS4GData(countMap, sortPositions)

        return Pair(ps4gDataList, sampleGameteCountMap)
    }

    /**
     * Function to process a set of MEMs for a given read.
     * This will do the full spline lookup and then create a consensus position, it will then increment the counter for
     * this gamete set/average position
     */
    fun processTempMEMs(
        tempMems: MutableList<MEM>,
        splineLookup: LinearLookupFunction,
        minMEMLength: Int,
        maxNumHits: Int,
        maxRange: Int,
        gameteToIdxMap: Map<String, Int>,
        countMap: MutableMap<Pair<Position, List<Int>>, Int>,
        sampleGameteCountMap: MutableMap<SampleGamete, Int>,
        gameteIdxToSampleGameteMap: Map<Int, SampleGamete>
    ) {
        val pairPosAndGameteSet =
            processMemsForRead(tempMems, splineLookup, minMEMLength, maxNumHits, maxRange, gameteToIdxMap)
        if (pairPosAndGameteSet.first.position != -1) {
            countMap[pairPosAndGameteSet] = countMap.getOrDefault(pairPosAndGameteSet, 0) + 1
            for (gameteIdx in pairPosAndGameteSet.second) { //Need to convert this to a set otherwise we get multiple counts for a given gamete
                sampleGameteCountMap[gameteIdxToSampleGameteMap[gameteIdx]!!] =
                    sampleGameteCountMap.getOrDefault(gameteIdxToSampleGameteMap[gameteIdx]!!, 0) + 1
            }
        }
    }

    /**
     * Function to process the MEMs collected for a given read.
     */
    fun processMemsForRead(tempMems: List<MEM>,
                           splineLookup: LinearLookupFunction,
                           minMEMLength: Int, maxNumHits: Int,
                           maxRange: Int,
                           gameteToIdxMap: Map<String, Int>): Pair<Position, List<Int>> {
        val bestHits = findBestMems(tempMems, minMEMLength, maxNumHits)
        if(bestHits.isEmpty()) {
            return Pair(Position("unknown",-1), listOf())
        }

        //convert MEMHits into Encoded Positions
        //TODO Rename as we are outside of int packing
        val referenceLookupPositions = lookupHitsToRefPosition(bestHits, splineLookup)

        //Create consensus Position
        val consensusPositionsAndGametes = createConsensusPositionAndGametes(referenceLookupPositions, gameteToIdxMap)

        //Filter on maximum range
        if(consensusPositionsAndGametes.third > maxRange) {
            return Pair(Position("unknown",-1), listOf())
        }

        return Pair(consensusPositionsAndGametes.first, consensusPositionsAndGametes.second)
    }

    /**
     * Function to do the Linear spline positional lookup from the RopeBWT3 MEMs
     * Returns a list of pairs of contig and Position for the reference lookups
     * Note these are still in binned positions (div 256)
     */
    fun lookupHitsToRefPosition(bestHits: List<MEMHit>, splineLookup: LinearLookupFunction) : List<Pair<String,Position>> {
        return bestHits.map { hit ->
            val position = Position(hit.contig, hit.pos)

            val lookupValue = splineLookup.value(position)

            Pair(hit.contig, lookupValue)
        }.filter { it.second != Position("unknown", 0) }
    }

    /**
     *  Function to find a consensus position for the gametes and output a Pair that can be used to increase counts
     */
    fun createConsensusPositionAndGametes(referenceLookupPositions: List<Pair<String,Position>>,gameteToIdxMap: Map<String, Int>) : Triple<Position, List<Int>, Int> {
        if(referenceLookupPositions.isEmpty()) {
            return Triple(Position("unknown",-1), listOf(), 0)
        }
        
        //determine chromosome majority
        val chrCounts = referenceLookupPositions.groupingBy { it.second.contig }.eachCount()
        val bestChromosome = chrCounts.maxBy { it.value }.key
        //remove hits that don't hit our chosen chromosome
        val bestHitsForChrom = referenceLookupPositions.filter { it.second.contig == bestChromosome }
        //compute the average position
        val averagePosition = bestHitsForChrom.sumOf { it.second.position} / bestHitsForChrom.size
        val rangePosition = bestHitsForChrom.maxOf{ it.second.position } - bestHitsForChrom.minOf{ it.second.position }

        //for now we just use the average position
        //Best chromosome is already in index form
        val binnedPosition = Position(bestChromosome, averagePosition)

        val gameteIndicesHit = bestHitsForChrom
            .map { it.first.split("_") }
            .map { it.last() } //needs to be last because there are scaffolds delimited by _
            .map {
                if(!gameteToIdxMap.containsKey(it)) {
                    myLogger.info("Gamete $it not found in. Chr: $bestChromosome, Position: $averagePosition")
                }
                gameteToIdxMap[it]!! }
            .toSortedSet()
            .toList()

        return Triple(binnedPosition, gameteIndicesHit, rangePosition)
    }

    /**
     * Function to remove MEMs that are not the longest and then return the best hits
     */
    fun findBestMems(tempMems: List<MEM>, minMEMLength: Int, maxNumHits: Int): List<MEMHit> {

        if(tempMems.isEmpty()) {
            return listOf()
        }
        //get the longest hits
        val maxLength = tempMems.maxOf { it.readEnd - it.readStart }
        if(maxLength < minMEMLength) {
            return listOf()
        }
        //remove any hits that are not the longest
        val bestHits = tempMems.filter { it.readEnd - it.readStart == maxLength }

        val totalNumHits = bestHits.sumOf { it.numHits }

        //If its below the maxNumHits then we can return all the hits as they are equally good
        return if(totalNumHits <= maxNumHits) {
            bestHits.flatMap { it.listMemHits }
        } else {
            listOf()
        }
    }
}