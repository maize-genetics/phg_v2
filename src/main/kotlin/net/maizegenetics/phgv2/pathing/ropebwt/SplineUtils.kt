package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import htsjdk.variant.vcf.VCFFileReader
import kotlinx.serialization.Serializable
import kotlinx.serialization.json.Json
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.parseALTHeader
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction
import org.apache.logging.log4j.LogManager
import java.io.*

@Serializable
data class SplineKnotLookup(
    val splineKnotMap: Map<String, Pair<DoubleArray, DoubleArray>>,
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
        fun buildSplineKnots(vcfDir: String, vcfType: String, minIndelLength: Int = 10, maxNumPointsPerChrom: Int = 250_000, contigSet : Set<String> = emptySet()) : SplineKnotLookup {
            val vcfFiles = buildVCFFileList(vcfDir, vcfType)
            val splineKnotMap = mutableMapOf<String, Pair<DoubleArray, DoubleArray>>()
            val chrIndexMap = mutableMapOf<String,Int>()
            val gameteIndexMap = mutableMapOf<String,Int>()
            for (vcfFile in vcfFiles!!) {
                myLogger.info("Reading ${vcfFile.name}")
                processVCFFileIntoSplineKnots(vcfFile, vcfType, splineKnotMap, chrIndexMap, gameteIndexMap, minIndelLength, maxNumPointsPerChrom, contigSet)
            }
            myLogger.info("Done reading VCF files")
            myLogger.info("Number of splines: ${splineKnotMap.size}")
            myLogger.info("Number of chromosomes: ${chrIndexMap.size}")
            myLogger.info("Number of gametes: ${gameteIndexMap.size}")
            return SplineKnotLookup(splineKnotMap, chrIndexMap, gameteIndexMap)
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
            splineKnotMap: MutableMap<String, Pair<DoubleArray, DoubleArray>>,
            chrIndexMap: MutableMap<String, Int>,
            gameteIndexMap: MutableMap<String, Int>,
            minIndelLength: Int=10,
            maxNumPointsPerChrom: Int = 250_000,
            contigSet: Set<String> = emptySet()
        ) {
            if(vcfType == "hvcf") {
                processHvcfFileIntoSplineKnots(vcfFile, splineKnotMap, chrIndexMap, gameteIndexMap, maxNumPointsPerChrom, contigSet)
            }
            else if(vcfType == "gvcf") {
                processGvcfFileIntoSplineKnots(vcfFile, splineKnotMap, chrIndexMap, gameteIndexMap, minIndelLength, maxNumPointsPerChrom, contigSet)
            }
            else {
                throw IllegalArgumentException("Unknown VCF type $vcfType")
            }
        }


        /**
         * Function to process a single HVCF file into the spline map for the PS4G file.
         */
        fun processHvcfFileIntoSplineKnots(
            hvcfFile: File?,
            splineMap: MutableMap<String, Pair<DoubleArray, DoubleArray>>,
            chrIndexMap : MutableMap<String,Int>,
            gameteIndexMap: MutableMap<String, Int>,
            maxNumPointsPerChrom: Int = 250_000,
            contigSet: Set<String> = emptySet()
        ) {
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

                myLogger.info ("Done reading HVCF file: ${hvcfFile?.name} with ${mapOfASMChrToListOfPoints.size} assembly chromosomes.")
                //build the splines
                //Downsample the number of points
                downsamplePoints(mapOfASMChrToListOfPoints, maxNumPointsPerChrom)

                //Now we need to build the splines for each of the assembly chromosomes
                //loop through each of the assembly coordinates and make splines for each
                for (entry in mapOfASMChrToListOfPoints.entries) {
                    val asmChr = entry.key
                    val listOfPoints = entry.value
                    checkMapAndAddToIndex(gameteIndexMap, sampleName)
                    //add to the splineList
                    buildSplineKnotsForASMChrom(listOfPoints, splineMap, asmChr, sampleName)
                }
            }
        }

        fun processGvcfFileIntoSplineKnots(
            gvcfFile: File?,
            splineMap: MutableMap<String, Pair<DoubleArray,DoubleArray>>,
            chrIndexMap : MutableMap<String,Int>,
            gameteIndexMap: MutableMap<String, Int>,
            minIndelLength: Int=10,
            maxNumPoints: Int = 250_000,
            contigSet: Set<String> = emptySet()
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
                flushBlock()

                myLogger.info ("Done reading GVCF file: ${gvcfFile?.name} with ${mapOfASMChrToListOfPoints.size} assembly chromosomes.")


                //Downsample the number of points
                downsamplePoints(mapOfASMChrToListOfPoints, maxNumPoints)

                //loop through each of the assembly coordinates and make splines for each
                for (entry in mapOfASMChrToListOfPoints.entries) {
                    val asmChr = entry.key
                    val listOfPoints = entry.value
                    myLogger.info("Building spline for $asmChr $sampleName")
                    checkMapAndAddToIndex(gameteIndexMap, sampleName)
                    buildSplineKnotsForASMChrom(listOfPoints, splineMap, asmChr, sampleName)
                }
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
            val numPointsToAdd = 4
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
         * Function to write the SplineKnotLookup to a file
         * This serializes the object using JSON
         */
        fun writeSplineLookupToFile(splineKnotLookup: SplineKnotLookup, outputFile:String) {
            bufferedWriter(outputFile).use { writer ->
                writer.write(Json.encodeToString(splineKnotLookup))
            }
        }


        /**
         * Function to load the spline knots from a file
         * This assumes that the SplineKnotLookup serialized file was serialized using JSON
         */
        fun loadSplineKnotLookupFromFile(inputFile: String): SplineKnotLookup {
            var splineKnotLookup: SplineKnotLookup

            //Json.decodeFromString<Data>
            bufferedReader(inputFile).use { reader ->
                splineKnotLookup = Json.decodeFromString(SplineKnotLookup.serializer(), reader.readText())
            }
            return splineKnotLookup
        }

        /**
         * Function that converts a Spline Knot Map into a Spline Map
         */
        fun convertKnotsToSpline(knots: Map<String, Pair<DoubleArray, DoubleArray>>) : Map<String, PolynomialSplineFunction> {
            val splineMap = mutableMapOf<String, PolynomialSplineFunction>()
            val interpolator = AkimaSplineInterpolator()
            knots.forEach { (key, value) ->
                val splineFunction = interpolator.interpolate(value.first, value.second)
                splineMap[key] = splineFunction
            }
            return splineMap
        }
    }
}