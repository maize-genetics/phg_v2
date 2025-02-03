package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.parseALTHeader
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction
import org.apache.logging.log4j.LogManager
import java.io.File

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


    override fun run() {
        logCommand(this)

        myLogger.info("Convert Ropebwt to PS4G")

        val (splineLookup, chrIndexMap, gameteIndexMap) = buildSplineLookup(hvcfDir)

        //build the output file name
        val outputFile = PS4GUtils.buildOutputFileName(ropebwtBed, outputDir)

        //build and write out the PS4G file
        buildPS4GData(ropebwtBed, splineLookup, chrIndexMap, gameteIndexMap,minMemLength, maxNumHits)


    }

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

    private fun processHvcfFileIntoSplines(
        hvcfFile: File?,
        splineMap: MutableMap<String, PolynomialSplineFunction>,
        chrIndexMap : MutableMap<String,Int>,
        gameteIndexMap: MutableMap<String, Int>
    ) {
        VCFFileReader(hvcfFile, false).use { reader ->
            val header = reader.header
            val headerParsed = parseALTHeader(header)
            val splineBuilder = SplineInterpolator()
            var currentChrom = ""
            val sampleName = reader.fileHeader.sampleNamesInOrder[0]
            checkMapAndAddToIndex(gameteIndexMap, sampleName)

            //This is ASM,REF
            val listOfPoints = mutableListOf<Pair<Double, Double>>()

            val iterator = reader.iterator()
            while (iterator.hasNext()) {
                val currentVariant = iterator.next()
                val chrom = currentVariant.contig
                if (chrom != currentChrom) {
                    if (listOfPoints.isEmpty()) {
                        currentChrom = chrom
                        continue
                    }
                    buildSpline(listOfPoints, splineBuilder, splineMap, currentChrom, sampleName)
                    currentChrom = chrom
                    //Check to see if we need to add a new entry in our chrIndexMap
                    checkMapAndAddToIndex(chrIndexMap, currentChrom)
                }

                val stPosition = currentVariant.start
                val endPosition = currentVariant.end

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
                    val asmEnd = regions.first().second.position
                    listOfPoints.add(Pair(asmStart.toDouble(), stPosition.toDouble()))
                    listOfPoints.add(Pair(asmEnd.toDouble(), endPosition.toDouble()))
                }

            }
            iterator.close()

            if (listOfPoints.isNotEmpty()) {
                buildSpline(listOfPoints, splineBuilder, splineMap, currentChrom, sampleName)
                //Check to see if we need to add a new entry in our chrIndexMap
                checkMapAndAddToIndex(chrIndexMap, currentChrom)
            }
        }
    }

    private fun checkMapAndAddToIndex(
        stringToIndexMap: MutableMap<String, Int>,
        sampleName: String
    ) {
        if (!stringToIndexMap.containsKey(sampleName)) {
            stringToIndexMap[sampleName] = stringToIndexMap.size
        }
    }

    private fun buildSpline(
        listOfPoints: MutableList<Pair<Double, Double>>,
        splineBuilder: SplineInterpolator,
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

    fun buildPS4GData(ropebwtBed: String,  splineLookup: Map<String, PolynomialSplineFunction>, chrIndexMap:Map<String,Int>,
                      gameteToIdxMap: Map<String,Int>,
                      minMEMLength: Int, maxNumHits: Int) : Map<Pair<Int,List<Int>>, Int> {
        val bedFileReader = bufferedReader(ropebwtBed)
        var currentLine = bedFileReader.readLine()
        val tempMems = mutableListOf<MEM>()
        val countMap = mutableMapOf<Pair<Int, List<Int>>, Int>()
        while (currentLine != null) {
            if(currentLine.isEmpty()) {
                currentLine = bedFileReader.readLine()
                continue
            }
            val alignmentParsed = RopeBWTUtils.parseStringIntoMem(currentLine)
            if (tempMems.isNotEmpty() && tempMems[0].readName != alignmentParsed.readName) {
                //write out the tempMems
                val pairPosAndGameteSet = processMemsForRead(tempMems, splineLookup, chrIndexMap,minMEMLength, maxNumHits)
                if(pairPosAndGameteSet.first != -1) {
                    countMap[pairPosAndGameteSet] = countMap.getOrDefault(pairPosAndGameteSet, 0) + 1
                }
                tempMems.clear()
            }
            tempMems.add(alignmentParsed)
            currentLine = bedFileReader.readLine()
        }

        val pairPosAndGameteSet = processMemsForRead(tempMems, splineLookup, chrIndexMap,minMEMLength, maxNumHits)
        if(pairPosAndGameteSet.first != -1) {
            countMap[pairPosAndGameteSet] = countMap.getOrDefault(pairPosAndGameteSet, 0) + 1
        }
        return countMap
    }

    fun processMemsForRead(tempMems: List<MEM>, splineLookup: Map<String, PolynomialSplineFunction>, chrIndexMap: Map<String, Int>, minMEMLength: Int, maxNumHits: Int): Pair<Int, List<Int>> {
        val bestHits = findBestMems(tempMems, minMEMLength, maxNumHits)
        if(bestHits.isEmpty()) {
            return Pair(-1, listOf())
        }

        //convert MEMHits into Encoded Positions
        val encodedPositions = encodeHitsToPosition(bestHits, splineLookup)

        //Create consensus Position
        val (consensusPosition, consensusGametes) = createConsensusPositionAndGametes(encodedPositions, chrIndexMap)

        TODO()

    }

    private fun encodeHitsToPosition(bestHits: List<MEMHit>, splineLookup: Map<String, PolynomialSplineFunction>) : List<Pair<String,Int>> {
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

    private fun createConsensusPositionAndGametes(encodedPositions: List<Pair<String,Int>>, chrIndexMap: Map<String,Int>) : Pair<Int, List<Int>> {
        val decodePositions = encodedPositions.map { Pair(it.first,PS4GUtils.decodePosition(it.second)) }

        //determine chromosome majority
        val bestChromosome = decodePositions.map { it.second }.maxOf { it.contig }
        //remove hits that don't hit our chosen chromosome
        val bestHitsForChrom = decodePositions.filter { it.second.contig == bestChromosome }
        //compute the average position
        val averagePosition = bestHitsForChrom.sumOf { it.second.position } / bestHitsForChrom.size
        //TODO remove hits that are too far from average...

        TODO()
        //Associate the gametes with the average positions
//        return Pair(encodedPositions)


    }

    /**
     * Function to remove MEMs that are not the longest and then return the best hits
     */
    fun findBestMems(tempMems: List<MEM>, minMEMLength: Int, maxNumHits: Int): List<MEMHit> {
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