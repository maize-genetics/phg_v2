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

//@Serializable
//data class SplineKnotLookup(
//    val splineKnotMap: Map<String, Pair<DoubleArray, DoubleArray>>,
//    val chrIndexMap: Map<String, Int>,
//    val gameteIndexMap: Map<String, Int>
//)
@Serializable
data class SplineKnotLookup(
    val splineKnotMap: Map<String, List<Triple<Int,String,Int>>>, //List of (asmPos, refChrom, refPos)
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
            val splineKnotMap = mutableMapOf<String, MutableList<Triple<Int,String,Int>>>()

            VCFFileReader(hvcfFile, false).use { reader ->
                val header = reader.header
                val headerParsed = parseALTHeader(header)
                val sampleName = reader.fileHeader.sampleNamesInOrder[0]
                checkMapAndAddToIndex(gameteIndexMap, sampleName)

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

                        addPointsToMap(splineKnotMap, asmStartChr, asmStart, chrom, stPosition)
                        addPointsToMap(splineKnotMap, asmEndChr, asmEnd, chrom, endPosition)
                    }

                }
                iterator.close()

                //build the splines
                //Downsample the number of points
                downsamplePointsByChrLength(splineKnotMap, numBpsPerKnot, randomSeed)

                checkMapAndAddToIndex(gameteIndexMap, sampleName)

                //Need to sort the points by asm position just in case
                for (entry in splineKnotMap.entries) {
                    val asmChr = entry.key
                    val sortedKnots = entry.value.sortedBy { it.first }.toMutableList()
                    splineKnotMap[asmChr] = sortedKnots
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

            val splineKnotMap = mutableMapOf<String, MutableList<Triple<Int,String,Int>>>()

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
                    if (blockAsmStart == blockAsmEnd || blockRefStart == blockRefEnd) {
                        addPointsToMap(splineKnotMap, blockAsmChr!!, blockAsmStart!!, currentRefChr!!, blockRefStart!!)
                    }
                    else {
                        //add both sets of points to the list
                        addPointsToMap(splineKnotMap, blockAsmChr!!, blockAsmStart!!, currentRefChr!!, blockRefStart!!)
                        addPointsToMap(splineKnotMap, blockAsmChr!!, blockAsmEnd!!, currentRefChr!!, blockRefEnd!!)
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

                            addPointsToMap(splineKnotMap, asmChr, asmPosStart, refChr, refPosStart)
                            addPointsToMap(splineKnotMap, asmChr, asmPosEnd, refChr, refPosEnd)
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
                                addPointsToMap(splineKnotMap, asmChr, asmPosMid, refChr, refPosStart)
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
                                addPointsToMap(splineKnotMap, asmChr, asmPosStart, refChr, refPosMid)
                            }
                        }
                        else -> {
                            flushBlock()
                            addPointsToMap(splineKnotMap, asmChr, asmPosStart, refChr, refPosStart)
                        }
                    }
                }
                flushBlock()

                //Downsample the number of points
                downsamplePointsByChrLength(splineKnotMap, numBpsPerKnot, randomSeed)

                checkMapAndAddToIndex(gameteIndexMap, sampleName)

                //Need to sort the points by asm position just in case
                for (entry in splineKnotMap.entries) {
                    val asmChr = entry.key
                    val sortedKnots = entry.value.sortedBy { it.first }.toMutableList()
                    splineKnotMap[asmChr] = sortedKnots
                }

            }

            return SplineKnotLookup(splineKnotMap, chrIndexMap, gameteIndexMap)
        }


        /**
         * This function will do the ref position binning as well as it adds the point to the map.
         */
        fun addPointsToMap(
            splineKnotMap: MutableMap<String, MutableList<Triple<Int,String, Int>>>,
            asmChr: String,
            asmPos: Int,
            refChrom: String,
            refPos: Int
        ) {
            val listOfPoints = splineKnotMap.getOrPut(asmChr) { mutableListOf() }
            listOfPoints.add(Triple(asmPos, refChrom, refPos/256)) //divide by 256 to bin up the reference positions
        }

        /**
         * Function to downsample the points in the spline taking into account the number of BPs per knot.
         *
         * First we determine how many knots we need to remove to hit the threshold then we randomly remove that many knots.
         * A lower numBpsPerKnot value will result in more knots being created which in turn will result in more detailed and accurate splines.
         */
        fun downsamplePointsByChrLength(splineKnotMap:MutableMap<String, MutableList<Triple<Int,String,Int>>>, numBpsPerKnot: Int = 50_000, randomSeed : Long = 12345) {
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
                val numKnots = (chromLength / numBpsPerKnot)
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
            listOfPoints: MutableList<Triple<Int,String,Int>>
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
            splineKnotMap: Map<String,List<Triple<Int,String, Int>>>,
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

            val allSplineKnots = mutableMapOf<String, List<Triple<Int, String, Int>>>()

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
                    val splineKnots = Json.decodeFromString<Map<String, List<Triple<Int,String,Int>>>>(reader.readText())
                    allSplineKnots.putAll(splineKnots)
                }
            }

            return SplineKnotLookup(
                allSplineKnots,
                indexMaps.chrIndexMap,
                indexMaps.gameteIndexMap
            )

        }
    }
}
