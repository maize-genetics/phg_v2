package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.headerCommand
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.parseALTHeader
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction
import org.apache.logging.log4j.LogManager
import java.io.File

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

    val hvcfDir by option(help = "Directory containing the hvcf files")
        .required()

    val minMemLength by option(help = "Minimum length of a match to be considered a match.")
        .int()
        .default(148)

    val maxNumHits by option(help = "Number of hits to report.  Note ropebwt can hit more than --max-num-hits but any alignment hitting more haplotypes than this will be ignored.")
        .int()
        .default(50)

    /**
     * Function to run the command.  This goes from a RopeBWT3 BED file for reads aligned against the whole chromosomes to a PS4G file.
     */
    override fun run() {
        logCommand(this)

        val command = headerCommand(this)

        myLogger.info("Convert Ropebwt to PS4G")

        myLogger.info("Building Spline Lookup")
        val (splineLookup, chrIndexMap, gameteIndexMap) = buildSplineLookup(hvcfDir)


        val sampleGameteIndexMap = gameteIndexMap.map { SampleGamete(it.key) to it.value}.toMap()

        myLogger.info("Building PS4G Output File Name")
        //build the output file name
        val outputFile = PS4GUtils.buildOutputFileName(ropebwtBed, outputDir)

        myLogger.info("Building PS4G Data")
        //build and write out the PS4G file
        val (ps4GData, sampleGameteCountMap) = buildPS4GData(ropebwtBed, splineLookup, chrIndexMap, gameteIndexMap,minMemLength, maxNumHits)

        myLogger.info("Writing out PS4G File")
        PS4GUtils.writeOutPS4GFile(ps4GData, sampleGameteCountMap, sampleGameteIndexMap, outputFile, listOf(), command)
    }

    /**
     * Function to build spline lookups from the hvcf files in the directory.
     * This will create Cubic splines based on the PolynomialSplineFunction from the Apache Commons Math3 library.
     */
    fun buildSplineLookup(hvcfDir: String) : Triple<Map<String, PolynomialSplineFunction>, Map<String,Int>, Map<String,Int>> {
        val hvcfFiles = File(hvcfDir).listFiles()
        val splineMap = mutableMapOf<String, PolynomialSplineFunction>()
        val chrIndexMap = mutableMapOf<String,Int>()
        val gameteIndexMap = mutableMapOf<String,Int>()
        for (hvcfFile in hvcfFiles!!) {
            processHvcfFileIntoSplines(hvcfFile, splineMap, chrIndexMap, gameteIndexMap)
        }
        return Triple(splineMap, chrIndexMap, gameteIndexMap)
    }

    /**
     * Function to process a single HVCF file into the spline map for the PS4G file.
     */
    fun processHvcfFileIntoSplines(
        hvcfFile: File?,
        splineMap: MutableMap<String, PolynomialSplineFunction>,
        chrIndexMap : MutableMap<String,Int>,
        gameteIndexMap: MutableMap<String, Int>
    ) {
        VCFFileReader(hvcfFile, false).use { reader ->
            val header = reader.header
            val headerParsed = parseALTHeader(header)
            val splineBuilder = AkimaSplineInterpolator()
            var currentChrom = ""
            val sampleName = reader.fileHeader.sampleNamesInOrder[0]
            checkMapAndAddToIndex(gameteIndexMap, sampleName)

            //This is ASM,REF
            val listOfPoints = mutableListOf<Pair<Double, Double>>()

            val iterator = reader.iterator()
            while (iterator.hasNext()) {
                val currentVariant = iterator.next()
                val chrom = currentVariant.contig
                //Check to see if we need to add a new entry in our chrIndexMap
                checkMapAndAddToIndex(chrIndexMap, chrom)
                if (chrom != currentChrom) {
                    if (listOfPoints.isEmpty()) {
                        currentChrom = chrom
                    } else {
                        buildSpline(listOfPoints, splineBuilder, splineMap, currentChrom, sampleName)
                        currentChrom = chrom
                    }
                }

                val stPosition = currentVariant.start
                val endPosition = currentVariant.end

                val stPositionEncoded = PS4GUtils.encodePosition(Position(chrom, stPosition), chrIndexMap)
                val endPositionEncoded = PS4GUtils.encodePosition(Position(chrom, endPosition), chrIndexMap)

                val genotype = currentVariant.genotypes
                    .get(0)
                    .alleles
                    .first()
                    .displayString
                    .replace("<", "")
                    .replace(">", "")

                headerParsed[genotype]?.let { altHeaderMetaData ->
                    val regions = altHeaderMetaData.regions
                    val asmStart = regions.first().first.position
                    val asmEnd = regions.last().second.position

                    listOfPoints.add(Pair(asmStart.toDouble(), stPositionEncoded.toDouble()))
                    listOfPoints.add(Pair(asmEnd.toDouble(), endPositionEncoded.toDouble()))
                }

            }
            iterator.close()

            if (listOfPoints.isNotEmpty()) {
                buildSpline(listOfPoints, splineBuilder, splineMap, currentChrom, sampleName)
            }
        }
    }

    /**
     * Simple function to check to see if a value already exists in our map and adds to it if not.
     */
    fun checkMapAndAddToIndex(
        stringToIndexMap: MutableMap<String, Int>,
        sampleName: String
    ) {
        if (!stringToIndexMap.containsKey(sampleName)) {
            stringToIndexMap[sampleName] = stringToIndexMap.size
        }
    }

    /**
     * Fnction to build the splines based on a list of Pair<Double, Double>
     */
    fun buildSpline(
        listOfPoints: MutableList<Pair<Double, Double>>,
        splineBuilder: AkimaSplineInterpolator,
        splineMap: MutableMap<String, PolynomialSplineFunction>,
        currentChrom: String,
        sampleName: String?
    ) {
        val sortedPoints = listOfPoints.sortedBy { it.first }.distinctBy { it.first }
        val asmArray = sortedPoints.map { it.first }.toDoubleArray()
        val refArray = sortedPoints.map { it.second }.toDoubleArray()
        val splineFunction = splineBuilder.interpolate(asmArray, refArray)
        splineMap["${currentChrom}_${sampleName}"] = splineFunction
        listOfPoints.clear()
    }

    /**
     * Function to build the output PS4G file.
     */
    fun buildPS4GData(ropebwtBed: String,  splineLookup: Map<String, PolynomialSplineFunction>, chrIndexMap:Map<String,Int>,
                      gameteToIdxMap: Map<String,Int>,
                      minMEMLength: Int, maxNumHits: Int) : Pair<List<PS4GData>, Map<SampleGamete,Int>> {

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

        val ps4gDataList = convertCountMapToPS4GData(countMap)


        return Pair(ps4gDataList, sampleGameteCountMap)
    }

    /**
     * Function to process a set of MEMs for a given read.
     * This will do the full spline lookup and then create a consensus position, it will then increment the counter for
     * this gamete set/average position
     */
    fun processTempMEMs(
        tempMems: MutableList<MEM>,
        splineLookup: Map<String, PolynomialSplineFunction>,
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
    fun processMemsForRead(tempMems: List<MEM>, splineLookup: Map<String, PolynomialSplineFunction>,
                           chrIndexMap: Map<String, Int>, minMEMLength: Int, maxNumHits: Int,
                           gameteToIdxMap: Map<String, Int>): Pair<Int, List<Int>> {
        val bestHits = findBestMems(tempMems, minMEMLength, maxNumHits)
        if(bestHits.isEmpty()) {
            return Pair(-1, listOf())
        }

        //convert MEMHits into Encoded Positions
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
     *  Function to find a consensus position for the gametes and output a Pair that can be used to increase counts
     */
    fun createConsensusPositionAndGametes(encodedPositions: List<Pair<String,Int>>, chrIndexMap: Map<String,Int>, gameteToIdxMap: Map<String, Int>) : Pair<Int, List<Int>> {
        if(encodedPositions.isEmpty()) {
            return Pair(-1, listOf())
        }
        val decodePositions = encodedPositions.map { Pair(it.first,PS4GUtils.decodePosition(it.second)) }

        //determine chromosome majority
        val bestChromosome = decodePositions.map { it.second }.maxOf { it.contig }
        //remove hits that don't hit our chosen chromosome
        val bestHitsForChrom = decodePositions.filter { it.second.contig == bestChromosome }
        //compute the average position
        val averagePosition = bestHitsForChrom.sumOf { it.second.position } / bestHitsForChrom.size
        //TODO future task remove hits that are too far from average...
        //for now we just use the average position
        //Best chromosome is already in index form
        val encodedPosition = PS4GUtils.encodePositionFromIdxAndPos(bestChromosome.toInt(), averagePosition)

        val gameteIndicesHit = bestHitsForChrom.map { it.first.split("_")[1] }.map { gameteToIdxMap[it]!! }.toSortedSet().toList()
        return Pair(encodedPosition,gameteIndicesHit )
        //Associate the gametes with the average positions
//        return Pair(encodedPositions)


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

    /**
     * Function to convert the count map to a PS4GData class for easy export
     */
    fun convertCountMapToPS4GData(countMap: Map<Pair<Int,List<Int>>, Int>) : List<PS4GData> {
        return countMap.map { (pair, count) ->
            PS4GData(pair.second, pair.first, count)
        }
    }
}