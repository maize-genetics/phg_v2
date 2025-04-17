package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.choice
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

    val vcfDir by option(help = "Directory containing the hvcf or gvcf files")
        .required()

    val vcfType by option(help = "Type of vcfs to build the splines")
        .choice("hvcf","gvcf")
        .default("hvcf")

    val minMemLength by option(help = "Minimum length of a match to be considered a match.")
        .int()
        .default(148)

    val maxNumHits by option(help = "Number of hits to report.  Note ropebwt can hit more than --max-num-hits but any alignment hitting more haplotypes than this will be ignored.")
        .int()
        .default(50)

    val sortPositions by option(help = "Sort positions in the resulting PS4G file.")
        .flag(default = true)

    val minIndelLength by option(help="Minimum length of an indel to break up the running block for spline creation of gvcfs.  If --vcf-type is hvcf this option is ignored.")
        .int()
        .default(10)

    val maxNumPointsPerChrom by option(help = "Number of points per chrom.  If there are more points for each sample's chromosomes we will downsample randomly..")
        .int()
        .default(250_000)

    /**
     * Function to run the command.  This goes from a RopeBWT3 BED file for reads aligned against the whole chromosomes to a PS4G file.
     */
    override fun run() {
        logCommand(this)

        val command = headerCommand(this)

        myLogger.info("Convert Ropebwt to PS4G")

        myLogger.info("Building Spline Lookup")
        val (splineLookup, chrIndexMap, gameteIndexMap) = buildSplineLookup(vcfDir, vcfType, minIndelLength, maxNumPointsPerChrom)


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
     * Function to build spline lookups from the hvcf files in the directory.
     * This will compute Cubic splines based on the Akima algorithm, as originally formulated by Hiroshi Akima, 1970
     * (http://doi.acm.org/10.1145/321607.321609) implemented via the Apache Commons Math3 library.
     */
    fun buildSplineLookup(hvcfDir: String, vcfType: String, minIndelLength: Int = 10, maxNumPointsPerChrom: Int = 250_000) : Triple<Map<String, PolynomialSplineFunction>, Map<String,Int>, Map<String,Int>> {
        val vcfFiles = buildVCFFileList(hvcfDir, vcfType)
        val splineMap = mutableMapOf<String, PolynomialSplineFunction>()
        val chrIndexMap = mutableMapOf<String,Int>()
        val gameteIndexMap = mutableMapOf<String,Int>()
        for (vcfFile in vcfFiles!!) {
            myLogger.info("Reading ${vcfFile.name}")
            processVCFFileIntoSplines(vcfFile, vcfType, splineMap, chrIndexMap, gameteIndexMap, minIndelLength, maxNumPointsPerChrom)
        }
        return Triple(splineMap, chrIndexMap, gameteIndexMap)
    }

    private fun buildVCFFileList(hvcfDir: String, vcfType: String): List<File> {
        return File(hvcfDir).listFiles()?.filter {
            fileMatchesType(it, vcfType)
        }?: throw Exception("No $vcfType files found in $hvcfDir")
    }

    private fun fileMatchesType(file: File, vcfType: String): Boolean {
        return if(vcfType == "hvcf") {
            file.name.endsWith("h.vcf") || file.name.endsWith("h.vcf.gz") ||
                    file.name.endsWith("hvcf") || file.name.endsWith("hvcf.gz")
        } else {
            file.name.endsWith("g.vcf") || file.name.endsWith("g.vcf.gz") ||
                    file.name.endsWith("gvcf") || file.name.endsWith("gvcf.gz")
        }
    }

    private fun processVCFFileIntoSplines(
        vcfFile: File,
        vcfType: String,
        splineMap: MutableMap<String, PolynomialSplineFunction>,
        chrIndexMap: MutableMap<String, Int>,
        gameteIndexMap: MutableMap<String, Int>,
        minIndelLength: Int=10,
        maxNumPointsPerChrom: Int = 250_000
    ) {
        if(vcfType == "hvcf") {
            processHvcfFileIntoSplines(vcfFile, splineMap, chrIndexMap, gameteIndexMap, maxNumPointsPerChrom)
        }
        else if(vcfType == "gvcf") {
            processGvcfFileIntoSplines(vcfFile, splineMap, chrIndexMap, gameteIndexMap, minIndelLength, maxNumPointsPerChrom)
        }
        else {
            throw Exception("Unknown VCF type $vcfType")
        }
    }

    /**
     * Function to process a single HVCF file into the spline map for the PS4G file.
     */
    fun processHvcfFileIntoSplines(
        hvcfFile: File?,
        splineMap: MutableMap<String, PolynomialSplineFunction>,
        chrIndexMap : MutableMap<String,Int>,
        gameteIndexMap: MutableMap<String, Int>,
        maxNumPointsPerChrom: Int = 250_000
    ) {
        VCFFileReader(hvcfFile, false).use { reader ->
            val header = reader.header
            val headerParsed = parseALTHeader(header)
            val splineBuilder = AkimaSplineInterpolator()
            val sampleName = reader.fileHeader.sampleNamesInOrder[0]
            checkMapAndAddToIndex(gameteIndexMap, sampleName)

            //This is ASM,REF
            val mapOfASMChrToListOfPoints = mutableMapOf<String, MutableList<Pair<Double, Double>>>()

            val iterator = reader.iterator()
            while (iterator.hasNext()) {
                val currentVariant = iterator.next()
                val chrom = currentVariant.contig
                //Check to see if we need to add a new entry in our chrIndexMap
                checkMapAndAddToIndex(chrIndexMap, chrom)

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

                    //We can treat each end point separately in case the assembly chromosome is different
                    val asmStartChr = regions.first().first.contig
                    val asmEndChr = regions.last().second.contig

                    addPointsToMap(mapOfASMChrToListOfPoints,asmStartChr, asmStart, stPositionEncoded)
                    addPointsToMap(mapOfASMChrToListOfPoints,asmEndChr, asmEnd, endPositionEncoded)
                }

            }
            iterator.close()

            //build the splines
            //Downsample the number of points
            downsamplePoints(mapOfASMChrToListOfPoints, maxNumPointsPerChrom)
            //Now we need to build the splines for each of the assembly chromosomes
            //loop through each of the assembly coordinates and make splines for each
            for (entry in mapOfASMChrToListOfPoints.entries) {
                val asmChr = entry.key
                val listOfPoints = entry.value
                checkMapAndAddToIndex(gameteIndexMap, sampleName)
                buildSpline(listOfPoints, splineBuilder, splineMap, asmChr, sampleName)
            }
        }
    }

    fun processGvcfFileIntoSplines(
        gvcfFile: File?,
        splineMap: MutableMap<String, PolynomialSplineFunction>,
        chrIndexMap : MutableMap<String,Int>,
        gameteIndexMap: MutableMap<String, Int>,
        minIndelLength: Int=10,
        maxNumPoints: Int = 250_000
    ) {
        var nextIndex = 0

        //This is ASM,REF
        val mapOfASMChrToListOfPoints = mutableMapOf<String, MutableList<Pair<Double, Double>>>()

        // Block tracking variables for regular (positive stranded) SNVs or END variants.
        var blockRefStart: Int? = null
        var blockRefEnd: Int? = null
        var blockAsmStart: Int? = null
        var blockAsmEnd: Int? = null
        var blockAsmChr: String? = null
        var blockAsmChrIdx: Int? = null
        var currentRefChr: String? = null

        // Flush the current regular block if it exists.
        fun flushBlock() {
            if (currentRefChr != null && blockRefStart != null && blockAsmStart != null &&
                blockAsmChr != null && blockAsmChrIdx != null
            ) {
                val refStPositionEncoded = PS4GUtils.encodePosition(Position(currentRefChr!!, blockRefStart!!), chrIndexMap)
                val refEndPositionEncoded = PS4GUtils.encodePosition(Position(currentRefChr!!, blockRefEnd!!), chrIndexMap)

                if (blockAsmStart == blockAsmEnd || blockRefStart == blockRefEnd) {
                    addPointsToMap(
                        mapOfASMChrToListOfPoints,
                        blockAsmChr!!,
                        blockAsmStart!!,
                        refStPositionEncoded
                    )
                }
                else {
                    //add both sets of points to the list
                    addPointsToMap(mapOfASMChrToListOfPoints, blockAsmChr!!, blockAsmStart!!, refStPositionEncoded)
                    addPointsToMap(mapOfASMChrToListOfPoints, blockAsmChr!!, blockAsmEnd!!, refEndPositionEncoded)
                }
            }
            blockRefStart = null
            blockRefEnd = null
            blockAsmStart = null
            blockAsmEnd = null
            blockAsmChr = null
            blockAsmChrIdx = null
        }

        // Helper function to update the regular block values.
        fun updateBlock(newRefStart: Int, newRefEnd: Int, newAsmStart: Int, newAsmEnd: Int, newAsmChr: String, newAsmChrIdx: Int) {
            if (blockAsmStart == null) {
                blockRefStart = newRefStart
                blockRefEnd = newRefEnd
                blockAsmStart = newAsmStart
                blockAsmEnd = newAsmEnd
                blockAsmChr = newAsmChr
                blockAsmChrIdx = newAsmChrIdx
            } else {
                blockRefEnd = newRefEnd
                blockAsmEnd = newAsmEnd
            }
        }

        // Open the gVCF file using HTSJDK.
        VCFFileReader(gvcfFile, false).use { reader ->
            for (variant in reader) {
                val refChr = variant.contig
                val refPosStart = variant.start
                val refPosEnd = variant.end

                //Check to see if we need to add a new entry in our chrIndexMap
                checkMapAndAddToIndex(chrIndexMap, refChr)

                // Flush regular block if the reference chromosome changes.
                if (currentRefChr != null && currentRefChr != refChr) {
                    //Here is where it deviates from the previous idea
                    //We need to flush the block because they are no longer on the same chromosome
                    flushBlock()
                }
                currentRefChr = refChr

                // Get ASM_Chr; if missing, default to "NA".
                val asmChr = variant.getAttribute("ASM_Chr", null)?.toString() ?: "NA"
                val refChrIndex = chrIndexMap.getOrPut(refChr) { nextIndex++ }

                // Parse ASM_Start and ASM_End safely.
                val asmPosStart = variant.getAttribute("ASM_Start", null)?.toString()?.toIntOrNull()
                val asmPosEnd = variant.getAttribute("ASM_End", null)?.toString()?.toIntOrNull()
                if (asmPosStart == null || asmPosEnd == null) {
                    myLogger.info("WARN - Skipping variant at $refChr:$refPosStart due to invalid ASM_Start/ASM_End.")
                    continue
                }

                // Ensure there is at least one alternate allele.
                if (variant.alternateAlleles.isEmpty()) {
                    myLogger.info("WARN - Skipping variant at $refChr:$refPosStart due to missing alternate allele.")
                    continue
                }

                // Get string attributes for length check (if > 1: indel; else: continue block)
                val altAllele = variant.getAlternateAllele(0).baseString
                val refAllele = variant.reference.baseString
                val hasEnd = variant.hasAttribute("END")

                // Determine the strand, defaulting to "+".
                val strand = variant.getAttribute("ASM_Strand", "+").toString()

                when {
                    // Negative stranded variant with an END attribute is written immediately as an inv_block.
                    (hasEnd && strand == "-") -> {
                        flushBlock()

                        val encodedRefStart = PS4GUtils.encodePosition(Position(refChr, refPosStart), chrIndexMap)
                        val encodedRefEnd = PS4GUtils.encodePosition(Position(refChr, refPosEnd), chrIndexMap)

                        addPointsToMap(mapOfASMChrToListOfPoints, asmChr, asmPosStart, encodedRefStart)
                        addPointsToMap(mapOfASMChrToListOfPoints, asmChr, asmPosEnd, encodedRefEnd)
                    }
                    // Collapse SNV and explicit END cases (for positive stranded variants).
                    (refAllele.length == 1 && (altAllele.length == 1 || hasEnd)) -> {
                        val effectiveRefEnd = if (hasEnd) refPosEnd else refPosStart
                        updateBlock(refPosStart, effectiveRefEnd, asmPosStart, asmPosEnd, asmChr, refChrIndex)
                    }
                    // Insertion: alternate allele is longer than one base.
                    (altAllele.length > 1) -> {
                        // Just update running block if small insertions are encountered
                        if (altAllele.length <= minIndelLength) {
                            val effectiveRefEnd = if (hasEnd) refPosEnd else refPosStart
                            updateBlock(refPosStart, effectiveRefEnd, asmPosStart, asmPosEnd, asmChr, refChrIndex)
                        } else {
                            // Report block and reset tracker
                            flushBlock()

                            // Generate midpoint to ensure monotonicity for spline
                            val asmPosMid = ((asmPosStart + asmPosEnd) * 0.5).toInt()

                            val encodedRefStart = PS4GUtils.encodePosition(Position(refChr, refPosStart), chrIndexMap)
                            addPointsToMap(mapOfASMChrToListOfPoints, asmChr, asmPosMid, encodedRefStart)
                        }
                    }
                    // Deletion: reference allele is longer than one base.
                    (refAllele.length > 1) -> {
                        // Just update running block if small deletions are encountered
                        if (refAllele.length <= minIndelLength) {
                            val effectiveRefEnd = if (hasEnd) refPosEnd else refPosStart
                            updateBlock(refPosStart, effectiveRefEnd, asmPosStart, asmPosEnd, asmChr, refChrIndex)
                        } else {
                            flushBlock()
                            val refPosMid = (((refPosStart + refAllele.length - 1) + refPosStart) * 0.5).toInt()
                            val encodedRefMidPoint = PS4GUtils.encodePosition(Position(refChr, refPosMid), chrIndexMap)
                            addPointsToMap(mapOfASMChrToListOfPoints, asmChr, asmPosStart, encodedRefMidPoint)
                        }
                    }
                    else -> {
                        flushBlock()
                        val encodedRefStart = PS4GUtils.encodePosition(Position(refChr, refPosStart), chrIndexMap)
                        addPointsToMap(mapOfASMChrToListOfPoints, asmChr, asmPosStart, encodedRefStart)
                    }
                }
            }
        }
        flushBlock()

        //Downsample the number of points
        downsamplePoints(mapOfASMChrToListOfPoints, maxNumPoints)
        //Now we need to build the splines for each of the assembly chromosomes
        val splineBuilder = AkimaSplineInterpolator()
        //loop through each of the assembly coordinates and make splines for each
        for (entry in mapOfASMChrToListOfPoints.entries) {
            val asmChr = entry.key
            val listOfPoints = entry.value
            val sampleName = gameteIndexMap.keys.firstOrNull { it.contains(asmChr) } ?: "NA"
            checkMapAndAddToIndex(gameteIndexMap, sampleName)
            buildSpline(listOfPoints, splineBuilder, splineMap, asmChr, sampleName)
        }
    }

    fun addPointsToMap(
        mapOfASMChrToListOfPoints: MutableMap<String, MutableList<Pair<Double, Double>>>,
        asmChr: String,
        asmPos: Int,
        refPos: Int
    ) {
        val listOfPoints = mapOfASMChrToListOfPoints.getOrPut(asmChr) { mutableListOf() }
        listOfPoints.add(Pair(asmPos.toDouble(), refPos.toDouble()))
    }

    fun downsamplePoints(map:MutableMap<String, MutableList<Pair<Double,Double>>>, maxNumPoints: Int = 250_000) {
        for (entry in map.entries) {
            val listOfPoints = entry.value
            val numPointsToRemove = listOfPoints.size - maxNumPoints

            if(numPointsToRemove > 0) {
                (0 until numPointsToRemove).forEach {
                    val randomIndex = (0 until listOfPoints.size).random()
                    listOfPoints.removeAt(randomIndex)
                }
            }
            map[entry.key] = listOfPoints
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
     * Function to build the splines based on a list of Pair<Double, Double>
     */
    fun buildSpline(
        listOfPoints: MutableList<Pair<Double, Double>>,
        splineBuilder: AkimaSplineInterpolator,
        splineMap: MutableMap<String, PolynomialSplineFunction>,
        currentChrom: String,
        sampleName: String?
    ) {
        val sortedPoints = listOfPoints.sortedBy { it.first }.distinctBy { it.first }
        if (sortedPoints.size <= 4) {
            myLogger.warn("Not enough points to build spline for $currentChrom $sampleName")
            listOfPoints.clear()
            return
        }
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

        val ps4gDataList = convertCountMapToPS4GData(countMap, sortPositions)


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

    /**
     * Function to convert the count map to a PS4GData class for easy export
     */
    fun convertCountMapToPS4GData(countMap: Map<Pair<Int,List<Int>>, Int>, sortPositions: Boolean = true) : List<PS4GData> {
        return if(sortPositions) {
            countMap.map { (pair, count) ->
                PS4GData(pair.second.sorted(), pair.first, count)
            }.sortedBy { PS4GUtils.decodePosition(it.pos) } //Need to decode it because chromosome might be in an unexpected order
        }
        else {
            countMap.map { (pair, count) ->
                PS4GData(pair.second, pair.first, count)
            }
        }
    }
}