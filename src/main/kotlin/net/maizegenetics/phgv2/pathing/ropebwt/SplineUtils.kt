package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import com.google.common.collect.TreeRangeMap
import htsjdk.variant.vcf.VCFFileReader
import kotlinx.serialization.Serializable
import kotlinx.serialization.json.Json
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.parseALTHeader
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction
import org.apache.logging.log4j.LogManager
import java.io.*
import kotlin.random.Random

@Serializable
data class SplineKnotLookup(
    val splineKnotMap: Map<String, Pair<DoubleArray, DoubleArray>>,
    val chrIndexMap: Map<String, Int>,
    val gameteIndexMap: Map<String, Int>
)

@Serializable
data class IndexMaps(
    val chrIndexMap: Map<String, Int>,
    val gameteIndexMap: Map<String, Int>
)

/**
 * Class to hold utility functions for building and saving splines built from hvcfs or gvcfs
 */
class SplineUtils{

    companion object {
        val myLogger = LogManager.getLogger(SplineUtils::class.java)

        /**
         * Function to build spline lookups from the hvcf files in the directory.
         * This will compute Cubic splines based on the Akima algorithm, as originally formulated by Hiroshi Akima, 1970
         * (http://doi.acm.org/10.1145/321607.321609) implemented via the Apache Commons Math3 library.
         * This does not build the full splines just yet but rather sets up the spline knots for each chromosome.
         */
        fun buildSplineKnots(vcfDir: String, vcfType: String, outputDir: String, minIndelLength: Int = 10, numBpsPerKnot: Int = 50_000, contigSet : Set<String> = emptySet(), randomSeed: Long = 12345) {
            val vcfFiles = buildVCFFileList(vcfDir, vcfType)

            var chrIndexMap = mutableMapOf<String,Int>()
            var gameteIndexMap = mutableMapOf<String,Int>()


            for (vcfFile in vcfFiles!!) {

                val splineOutputFile = "${outputDir}/${vcfFile.nameWithoutExtension}_spline_knots.json.gz"


                myLogger.info("Reading ${vcfFile.name}")
                val splineKnotLookup = processVCFFileIntoSplineKnots(vcfFile, vcfType, chrIndexMap, gameteIndexMap, minIndelLength, numBpsPerKnot, contigSet, randomSeed)

                myLogger.info("Done processing ${vcfFile.name}")
                myLogger.info("Number of splines: ${splineKnotLookup.splineKnotMap.size}")
                myLogger.info("Number of chromosomes: ${splineKnotLookup.chrIndexMap.size}")
                myLogger.info("Number of gametes: ${splineKnotLookup.gameteIndexMap.size}")

                //update the index maps
                chrIndexMap = splineKnotLookup.chrIndexMap.toMutableMap()
                gameteIndexMap = splineKnotLookup.gameteIndexMap.toMutableMap()

                //need to export the spline knots to a file
                myLogger.info("Writing spline knots to $splineOutputFile")
                writeSplineKnotsToFile(splineKnotLookup.splineKnotMap, splineOutputFile)
            }
            myLogger.info("Done reading VCF files")
            myLogger.info("Total number of chromosomes: ${chrIndexMap.keys.size}")
            myLogger.info("Total number of gametes: ${gameteIndexMap.keys.size}")

            val outputIndexFile = "${outputDir}/index_maps.json.gz"

            myLogger.info("Writing out Chromosome and Gamete Index Maps to: $outputIndexFile")
            writeIndexMapsToFile(IndexMaps(chrIndexMap, gameteIndexMap), outputIndexFile)
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

        private fun processVCFFileIntoSplineKnots(
            vcfFile: File,
            vcfType: String,
            chrIndexMap: MutableMap<String, Int>,
            gameteIndexMap: MutableMap<String, Int>,
            minIndelLength: Int=10,
            numBpsPerKnot: Int = 50_000,
            contigSet: Set<String> = emptySet(),
            randomSeed: Long = 12345
        ) : SplineKnotLookup {
            return if(vcfType == "hvcf") {
                processHvcfFileIntoSplineKnots(vcfFile, chrIndexMap, gameteIndexMap, numBpsPerKnot, contigSet, randomSeed)
            } else if(vcfType == "gvcf") {
                processGvcfFileIntoSplineKnots(vcfFile,chrIndexMap, gameteIndexMap, minIndelLength, numBpsPerKnot, contigSet, randomSeed)
            } else {
                throw IllegalArgumentException("Unknown VCF type $vcfType")
            }
        }


        /**
         * Function to process a single HVCF file into the spline map for the PS4G file.
         */
        fun processHvcfFileIntoSplineKnots(
            hvcfFile: File?,
            chrIndexMap: MutableMap<String, Int>,
            gameteIndexMap: MutableMap<String, Int>,
            numBpsPerKnot: Int = 50_000,
            contigSet: Set<String> = emptySet(),
            randomSeed: Long = 12345
        ) : SplineKnotLookup {
            val splineKnotMap = mutableMapOf<String, Pair<DoubleArray, DoubleArray>>()

            VCFFileReader(hvcfFile, false).use { reader ->
                val header = reader.header
                val headerParsed = parseALTHeader(header)
                val sampleName = reader.fileHeader.sampleNamesInOrder[0]
                checkMapAndAddToIndex(gameteIndexMap, sampleName)

                //This is ASM,REF
                val mapOfASMChrToListOfPoints = mutableMapOf<String, MutableList<Pair<Double, Double>>>()

                val iterator = reader.iterator()
                while (iterator.hasNext()) {
                    val currentVariant = iterator.next()
                    val chrom = currentVariant.contig

                    if(contigSet.isNotEmpty() && !contigSet.contains(chrom)) {
                        continue
                    }

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
                downsamplePointsByChrLength(mapOfASMChrToListOfPoints, numBpsPerKnot, randomSeed)

                //Now we need to build the splines for each of the assembly chromosomes
                //loop through each of the assembly coordinates and make splines for each
                for (entry in mapOfASMChrToListOfPoints.entries) {
                    val asmChr = entry.key
                    val listOfPoints = entry.value
                    checkMapAndAddToIndex(gameteIndexMap, sampleName)
                    //add to the splineList
                    buildSplineKnotsForASMChrom(listOfPoints, splineKnotMap, asmChr, sampleName)
                }
            }
            return SplineKnotLookup(splineKnotMap, chrIndexMap, gameteIndexMap)
        }

        fun processGvcfFileIntoSplineKnots(
            gvcfFile: File?,
            chrIndexMap: MutableMap<String, Int>,
            gameteIndexMap: MutableMap<String, Int>,
            minIndelLength: Int=10,
            numBpsPerKnot: Int = 50_000,
            contigSet: Set<String> = emptySet(),
            randomSeed: Long = 12345
        ) : SplineKnotLookup {

            val splineKnotMap = mutableMapOf<String, Pair<DoubleArray, DoubleArray>>()

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
                val sampleName = reader.fileHeader.sampleNamesInOrder[0]

                for (variant in reader) {
                    val refChr = variant.contig
                    if (contigSet.isNotEmpty() && !contigSet.contains(refChr)) {
                        continue
                    }
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
                    checkMapAndAddToIndex(chrIndexMap, refChr)
                    val refChrIndex = chrIndexMap[refChr]
                        ?: throw IllegalStateException("Reference chromosome $refChr not found in chrIndexMap.")

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
                flushBlock()

                //Downsample the number of points
                downsamplePointsByChrLength(mapOfASMChrToListOfPoints, numBpsPerKnot, randomSeed)

                //loop through each of the assembly coordinates and make splines for each
                for (entry in mapOfASMChrToListOfPoints.entries) {
                    val asmChr = entry.key
                    val listOfPoints = entry.value
                    myLogger.info("Building spline for $asmChr $sampleName")
                    checkMapAndAddToIndex(gameteIndexMap, sampleName)

                    buildSplineKnotsForASMChrom(listOfPoints, splineKnotMap, asmChr, sampleName)
                }
            }

            return SplineKnotLookup(splineKnotMap, chrIndexMap, gameteIndexMap)
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

        /**
         * Function to downsample the points so the splines are not too big for each contig.  This is randomly selected using the randomSeed provided by the user or a default of 12345.
         *
         * We first determine how many points need to be removed then we randomly select that many indices to remove from the list of points for each chromosome.
         *
         * If there are fewer points in the list for a given chromosome than the maxNumPoints, it will not downsample.
         */
        @Deprecated("Use downsamplePointsByChrLength instead for more accurate downsampling based on chromosome length.", replaceWith = ReplaceWith(
            "downsamplePointsByChrLength(splineKnotMap, maxNumPoints, randomSeed)"
        ))
        fun downsamplePoints(splineKnotMap:MutableMap<String, MutableList<Pair<Double,Double>>>, maxNumPoints: Int = 250_000, randomSeed : Long = 12345) {
            for (entry in splineKnotMap.entries) {
                val listOfPoints = entry.value
                val numPointsToRemove = listOfPoints.size - maxNumPoints

                if(numPointsToRemove > 0) {
                    //Build a list of unique indices to remove
                    val indicesToRemove = buildRemoveIndexSet(randomSeed, numPointsToRemove, listOfPoints)

                    //Filter out any of the indices that are in the set and write to the map overwriting the existing contig.
                    listOfPoints.filterIndexed { index, pair -> !indicesToRemove.contains(index) }
                        .let { filteredList ->
                            splineKnotMap[entry.key] = filteredList.toMutableList()
                        }

                }
                else {
                    // If there are fewer points than the maxNumPoints, we do not need to downsample
                    splineKnotMap[entry.key] = listOfPoints
                }
            }
        }

        /**
         * Function to downsample the points in the spline taking into account the number of BPs per knot.
         *
         * First we determine how many knots we need to remove to hit the threshold then we randomly remove that many knots.
         * A lower numBpsPerKnot value will result in more knots being created which in turn will result in more detailed and accurate splines.
         */
        fun downsamplePointsByChrLength(splineKnotMap:MutableMap<String, MutableList<Pair<Double,Double>>>, numBpsPerKnot: Int = 50_000, randomSeed : Long = 12345) {
            //numBpsPerKnot is the number of base pairs per knot
            for (entry in splineKnotMap.entries) {
                val listOfPoints = entry.value

                if(listOfPoints.isEmpty()) {
                    //If there are no points for this chromosome, skip it
                    myLogger.info("Skipping ${entry.key} as there are no points to downsample.")
                    continue
                }

                //get the length of the chromosome.  Should be the max of the first values in the pairs
                val chromLength = listOfPoints.maxOf { it.first }

                //Calculate the number of knots we need to keep.  If numBpsPerKnot is greater than the chromosome length, we do not need to downsample.
                val numKnots = (chromLength / numBpsPerKnot).toInt()
                if( chromLength <= numBpsPerKnot || numKnots <= 1) {
                    //If the chromosome length is less than the number of base pairs per knot, we do not need to downsample
                    continue
                }


                //Calculate the number of points to remove
                val numPointsToRemove = listOfPoints.size - numKnots
                if(numPointsToRemove > 0) {
                    val indicesToRemove = buildRemoveIndexSet(randomSeed, numPointsToRemove, listOfPoints)

                    //Filter out any of the indices that are in the set and write to the map overwriting the existing contig.
                    listOfPoints.filterIndexed { index, pair -> !indicesToRemove.contains(index) }
                        .let { filteredList ->
                            splineKnotMap[entry.key] = filteredList.toMutableList()
                        }

                }
            }
        }

        private fun buildRemoveIndexSet(
            randomSeed: Long,
            numPointsToRemove: Int,
            listOfPoints: MutableList<Pair<Double, Double>>
        ): MutableSet<Int> {
            //Build a list of unique indices to remove
            val indicesToRemove = mutableSetOf<Int>()
            val random = Random(randomSeed)
            // Loop until we have enough unique indices to remove
            while (indicesToRemove.size < numPointsToRemove) {
                val randomIndex = random.nextInt(listOfPoints.size)
                if (!indicesToRemove.contains(randomIndex)) {
                    indicesToRemove.add(randomIndex)
                }
            }
            return indicesToRemove
        }

        /**
         * Function to build the splines based on a list of Pair<Double, Double>
         */
        fun buildSplineKnotsForASMChrom(
            listOfPoints: MutableList<Pair<Double, Double>>,
            splineKnotMap: MutableMap<String, Pair<DoubleArray, DoubleArray>>,
            currentChrom: String,
            sampleName: String?
        ) {

            val sortedPoints = sortAndSplitPoints(listOfPoints)
            val asmArray = sortedPoints.map { it.first }.toDoubleArray()
            val refArray = sortedPoints.map { it.second }.toDoubleArray()
            splineKnotMap["${currentChrom}_${sampleName}"] = Pair(asmArray, refArray)
            listOfPoints.clear()
        }

        fun sortAndSplitPoints(listOfPoints: MutableList<Pair<Double, Double>>) : List<Pair<Double, Double>> {
            val sortedPoints =  listOfPoints.sortedBy { it.first }.distinctBy { it.first }

            if(sortedPoints.size <= 4) {
                val outputPoints = mutableListOf<Pair<Double, Double>>()
                //drop that many points between each set of points if the ref chroms match
                sortedPoints.zipWithNext().map { (firstPair, secondPair) ->
                    val totalPos = secondPair.first - firstPair.first + 1

                    val firstRef = PS4GUtils.decodePosition(firstPair.second.toInt())
                    val secondRef = PS4GUtils.decodePosition(secondPair.second.toInt())

                    //Add the first point
                    outputPoints.add(firstPair)

                    if(firstRef.contig == secondRef.contig) {
                            addIntermediateSplineKnots(firstPair, totalPos, firstRef, secondRef, outputPoints)
                    }
                }
                //Add the last point
                if(sortedPoints.isNotEmpty()) {
                    outputPoints.add(sortedPoints.last())
                }
                return outputPoints
            }
            else {
                return sortedPoints
            }
        }


        /**
         * Function to add the spline knots to the output points.
         * This makes sure there are at least 4 points for each spline so the Akima spline function will work.
         */
        private fun addIntermediateSplineKnots(
            firstPair: Pair<Double, Double>,
            totalPos: Double,
            firstRef: Position,
            secondRef: Position,
            outputPoints: MutableList<Pair<Double, Double>>
        ) {
            //Contigs are the same so we want to split the region into 4 parts.
            // We need 4 as it is the minumum number of knots that the Akima spline function requires.
            // If there are less it will throw an error.
            val numPointsToAdd = minOf(4, totalPos.toInt()-1)
            for (i in 1 until numPointsToAdd) {
                val newPos = firstPair.first + (totalPos * (i.toDouble()) / (numPointsToAdd))

                val newRefPos = if (firstRef.position <= secondRef.position) {
                    //positive strand
                    firstRef.position + (totalPos * ((i.toDouble()) / (numPointsToAdd)))
                } else {
                    //negative strand
                    firstRef.position - (totalPos * ((i.toDouble()) / (numPointsToAdd)))
                }
                val newRef = PS4GUtils.encodePositionNoLookup(Position(secondRef.contig, newRefPos.toInt()))
                outputPoints.add(Pair(newPos, newRef.toDouble()))
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
         * Function to write out the spline knot map to a file.
         */
        fun writeSplineKnotsToFile(
            splineKnotMap: Map<String, Pair<DoubleArray, DoubleArray>>,
            outputFile: String
        ) {
            bufferedWriter(outputFile).use { writer ->
                writer.write(Json.encodeToString(splineKnotMap))
            }
        }

        /**
         * File to write out the global index maps to a file.
         */
        fun writeIndexMapsToFile(
            indexMaps: IndexMaps,
            outputFile: String
        ) {
            bufferedWriter(outputFile).use { writer ->
                writer.write(Json.encodeToString(indexMaps))
            }
        }



        /**
         * Function to load in the spline knot lookups using the directory of a spline file for each assembly.
         * This adds a marginal amount of time to processing but allows for more flexibility and faster writes.
         */
        fun loadSplineKnotLookupFromDirectory(inputDir: String) : SplineKnotLookup {
            //Check to see if there is a directory
            val inputDirectory = File(inputDir)
            if (!inputDirectory.isDirectory) {
                throw IllegalArgumentException("Input directory $inputDir is not a directory")
            }



            //load in the indexMaps
            //index_maps.json.gz
            val indexMapsFile = "${inputDir}/index_maps.json.gz"

            //Check to see if we have index maps
            if (!File(indexMapsFile).exists()) {
                throw IllegalArgumentException("Index maps file $indexMapsFile does not exist")
            }


            var indexMaps: IndexMaps
            bufferedReader(indexMapsFile).use { reader ->
                indexMaps = Json.decodeFromString(IndexMaps.serializer(), reader.readText())
            }

            val allSplineKnots = mutableMapOf<String, Pair<DoubleArray, DoubleArray>>()

            //get list of all spline knot json files
            //"${outputDir}/${vcfFile.nameWithoutExtension}_spline_knots.json.gz"

            val splineKnotFiles = File(inputDir).listFiles { file ->
                file.isFile && file.name.endsWith("_spline_knots.json.gz")
            } ?: throw IllegalArgumentException("No spline knot files found in directory $inputDir")

            //Make sure we have at least one spline knot file to process
            if (splineKnotFiles.isEmpty()) {
                throw IllegalArgumentException("No spline knot files found in directory $inputDir")
            }

            for(splineKnotFile in splineKnotFiles) {
                bufferedReader(splineKnotFile.toString()).use { reader ->
                    val splineKnots = Json.decodeFromString<Map<String, Pair<DoubleArray, DoubleArray>>>(reader.readText())
                    allSplineKnots.putAll(splineKnots)
                }
            }

            return SplineKnotLookup(
                allSplineKnots,
                indexMaps.chrIndexMap,
                indexMaps.gameteIndexMap
            )

        }

        fun insureMonotonicity(
            asmPositions: DoubleArray,
            refPositions: DoubleArray
        ) : Pair<DoubleArray, DoubleArray> {
            val removeSet = mutableSetOf<Int>()
            for(idx in 1 until asmPositions.size) {
                if(asmPositions[idx] <= asmPositions[idx - 1]) {
                    removeSet.add(idx)
                }
            }

            return if(removeSet.isNotEmpty()) {
                val newAsmPositions = asmPositions.filterIndexed { index, _ -> !removeSet.contains(index) }.toDoubleArray()
                val newRefPositions = refPositions.filterIndexed { index, _ -> !removeSet.contains(index) }.toDoubleArray()
                if(newAsmPositions.size != newRefPositions.size) {
                    throw IllegalStateException("ASM and REF positions are not the same size after removing non-monotonic points.")
                }
                Pair(newAsmPositions, newRefPositions)
            }
            else {
                Pair(asmPositions, refPositions)
            }
        }

        /**
         * Function that converts a Spline Knot Map into a Spline Map
         */
        fun convertKnotsToSpline(knots: Map<String, Pair<DoubleArray, DoubleArray>>) : Map<String, PolynomialSplineFunction> {
            val splineMap = mutableMapOf<String, PolynomialSplineFunction>()
            val interpolator = AkimaSplineInterpolator()
            knots.forEach { (key, value) ->
                println("Processing spline for $key")
                val asmPositions = value.first
                val refPositions = value.second
//                val (asmPositions, refPositions) = insureMonotonicity(value.first, value.second)
                val splineFunction = interpolator.interpolate(asmPositions,refPositions)
                splineMap[key] = splineFunction
            }
            return splineMap
        }

        fun convertKnotsToLinearSpline(knots: Map<String,Pair<DoubleArray,DoubleArray>>) : RangeMap<Position,Pair<Position,Position>> {
            //Need to make the splines into linear blocks

            //loop through each assembly chromosome:
            val knotMap = TreeRangeMap.create<Position,Pair<Position,Position>>()
            knots.forEach { (key, value) ->
                val asmPositions = value.first
                val refPositions = value.second

                if(asmPositions.size != refPositions.size) {
                    throw IllegalStateException("ASM and REF positions are not the same size for $key")
                }

                //We need to create a range for each pair of asm and ref positions
                for (i in 0 until asmPositions.size - 1) {
                    val asmStart = Position(key, asmPositions[i].toInt())
                    val asmEnd = Position(key, asmPositions[i + 1].toInt())

                    //Convert the ref positions to Position objects as they are encoded
                    val refStartPos = PS4GUtils.decodePosition(refPositions[i].toInt())
                    val refEndPos = PS4GUtils.decodePosition(refPositions[i + 1].toInt())

                    if(refStartPos.contig != refEndPos.contig) {
                        continue // Skip if the contigs are not the same  We do not want a spline between them
                    }

                    knotMap.put(Range.closed(asmStart, asmEnd), Pair(refStartPos, refEndPos))
                }
            }
            return knotMap
        }

        fun convertKnotMapToLinearLookupFunction(
            knotMap: RangeMap<Position,Pair<Position,Position>>
        ) : LinearLookupFunction {
            return LinearLookupFunction(knotMap)
        }
    }
}
