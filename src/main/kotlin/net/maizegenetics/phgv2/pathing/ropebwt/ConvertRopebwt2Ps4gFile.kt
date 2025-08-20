package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.headerCommand
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.utils.Position
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction
import org.apache.logging.log4j.LogManager

/**
 * This class will convert a RopebwtBed file to a PS4G file.  It will only work with ropebwt3 files where the reads are
 * aligned to the whole assembly chromosomes using the mem command.  MEMs are Maximal Exact Matches and are used to
 * determine what the optimal mapping is.  One downside to this approach is that if there is a SNP in the middle of the read, the mappings will be ignored.
 * We can run SW in the future to fix this but it is very slow.
 */
class ConvertRopebwt2Ps4gFile : CliktCommand(help = "Convert RopebwtBed to PS4G") {

    val myLogger = LogManager.getLogger(ConvertRopebwt2Ps4gFile::class.java)

    val ropebwtBed by option(help = "RopebwtBed file")
        .required()

    val outputDir by option(help = "Output directory")
        .required()

   val splineKnotDir by option(help = "Spline Knot Directory")
        .required()

    val minMemLength by option(help = "Minimum length of a match to be considered a match.")
        .int()
        .default(148)

    val maxNumHits by option(help = "Number of hits to report.  Note ropebwt can hit more than --max-num-hits but any alignment hitting more haplotypes than this will be ignored.")
        .int()
        .default(50)

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

        val linearSplineKnots = SplineUtils.convertKnotsToLinearSpline(splineKnots)
        val splineLookup = SplineUtils.convertKnotMapToLinearLookupFunction(linearSplineKnots)

        val sampleGameteIndexMap = gameteIndexMap.map { SampleGamete(it.key) to it.value}.toMap()

        myLogger.info("Building PS4G Output File Name")
        //build the output file name
        val outputFile = PS4GUtils.buildOutputFileName(ropebwtBed, outputDir)

        myLogger.info("Building PS4G Data")
        //build and write out the PS4G file
        val (ps4GData, sampleGameteCountMap) = buildPS4GData(ropebwtBed, splineLookup, chrIndexMap, gameteIndexMap,minMemLength, maxNumHits, sortPositions)

        myLogger.info("Writing out PS4G File")
        PS4GUtils.writeOutPS4GFile(ps4GData, sampleGameteCountMap, sampleGameteIndexMap, outputFile, listOf(), command)
    }

    /**
     * Function to build the output PS4G file.
     */
    fun buildPS4GData(ropebwtBed: String,
                      splineLookup: LinearLookupFunction,
                      chrIndexMap:Map<String,Int>,
                      gameteToIdxMap: Map<String,Int>,
                      minMEMLength: Int, maxNumHits: Int,
                      sortPositions: Boolean = true) : Pair<List<PS4GData>, Map<SampleGamete,Int>> {

        val gameteIdxToSampleGameteMap = gameteToIdxMap.map { it.value to SampleGamete(it.key) }.toMap()
        val bedFileReader = bufferedReader(ropebwtBed)
        var currentLine = bedFileReader.readLine()
        val tempMems = mutableListOf<MEM>()
        val countMap = mutableMapOf<Pair<Int, List<Int>>, Int>()
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
                    chrIndexMap,
                    minMEMLength,
                    maxNumHits,
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
            chrIndexMap,
            minMEMLength,
            maxNumHits,
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
        chrIndexMap: Map<String, Int>,
        minMEMLength: Int,
        maxNumHits: Int,
        gameteToIdxMap: Map<String, Int>,
        countMap: MutableMap<Pair<Int, List<Int>>, Int>,
        sampleGameteCountMap: MutableMap<SampleGamete, Int>,
        gameteIdxToSampleGameteMap: Map<Int, SampleGamete>
    ) {
        val pairPosAndGameteSet =
            processMemsForRead(tempMems, splineLookup, chrIndexMap, minMEMLength, maxNumHits, gameteToIdxMap)
        if (pairPosAndGameteSet.first != -1) {
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
                           chrIndexMap: Map<String, Int>,
                           minMEMLength: Int, maxNumHits: Int,
                           gameteToIdxMap: Map<String, Int>): Pair<Int, List<Int>> {
        val bestHits = findBestMems(tempMems, minMEMLength, maxNumHits)
        if(bestHits.isEmpty()) {
            return Pair(-1, listOf())
        }

        //convert MEMHits into Encoded Positions
        //TODO Rename as we are outside of int packing
        val encodedPositions = encodeHitsToPosition(bestHits, splineLookup)

        //Create consensus Position
        return createConsensusPositionAndGametes(encodedPositions, chrIndexMap, gameteToIdxMap)
    }

    /**
     * Function to do the spline positional lookup from the RopeBWT3 MEMs
     */
    fun encodeHitsToPosition(bestHits: List<MEMHit>, splineLookup: Map<String, PolynomialSplineFunction>) : List<Pair<String,Int>> {
        return bestHits.map { hit ->
            //hit.contig is chr_SampleGamete
            val spline = splineLookup[hit.contig] ?: return@map Pair(hit.contig, -1)
            val position = hit.pos
            if (spline.isValidPoint(position.toDouble())) {
                Pair(hit.contig, spline.value(position.toDouble()).toInt())
            } else {
                Pair(hit.contig, -1)
            }
        }.filter { it.second != -1 }
    }

    /**
     * Function to do the spline positional lookup from the RopeBWT3 MEMs
     */
    fun encodeHitsToPosition(bestHits: List<MEMHit>, splineLookup: LinearLookupFunction) : List<Pair<String,Int>> {
        return bestHits.map { hit ->
            val position = Position(hit.contig, hit.pos)

            val lookupValue = splineLookup.value(position)

            if (lookupValue == Position("unknown", 0)) {
                Pair(hit.contig, -1)
            }
            else {
                Pair(hit.contig, lookupValue.position)
            }
        }.filter { it.second != -1 }
    }

    /**
     *  Function to find a consensus position for the gametes and output a Pair that can be used to increase counts
     */
    fun createConsensusPositionAndGametes(decodedPositions: List<Pair<String,Int>>,chrIndexMap: Map<String,Int> ,gameteToIdxMap: Map<String, Int>) : Pair<Int, List<Int>> {
        if(decodedPositions.isEmpty()) {
            return Pair(-1, listOf())
        }

        //determine chromosome majority
        val chrCounts = decodedPositions.groupingBy { it.first.substringBeforeLast("_") }.eachCount()
        val bestChromosome = chrCounts.maxBy { it.value }.key
        //remove hits that don't hit our chosen chromosome
        val bestHitsForChrom = decodedPositions.filter { it.first.substringBeforeLast("_") == bestChromosome }
        //compute the average position
        val averagePosition = bestHitsForChrom.sumOf { it.second / bestHitsForChrom.size }
        //TODO future task remove hits that are too far from average...
        //for now we just use the average position
        //Best chromosome is already in index form
        val bestChromosomeIndex = chrIndexMap[bestChromosome.substringBeforeLast("_")] ?: throw IllegalArgumentException("Chromosome $bestChromosome not found in chrIndexMap")
        val encodedPosition = PS4GUtils.encodePositionFromIdxAndPos(bestChromosomeIndex, averagePosition)

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

        if(encodedPosition < 0) {
            myLogger.warn("Encoded position is negative: $encodedPosition for chromosome $bestChromosome at position $averagePosition")
        }
        return Pair(encodedPosition,gameteIndicesHit )
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