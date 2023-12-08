package net.maizegenetics.phgv2.utils

import it.unimi.dsi.fastutil.longs.Long2LongOpenHashMap
import net.maizegenetics.phgv2.api.HaplotypeGraph
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.*

/**
 * data class to hold the KmerMap information.  From the PHGv1 source code:
 */
data class KmerMapData(val haplotypeListId: Int, val rangeToBitSetMap: Map<Int, BitSet>, val kmerHashToLongMap: Long2LongOpenHashMap)


private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.utils.AlignmentUtils")


fun alignReadsToHaplotypes(hvcfDir: String, kmerIndexFile:String, readFiles:String, paired:Boolean, outputDir:String) {
    //loop through all files in hvcfDir and create a list of hvcf files
    val hvcfFiles = File(hvcfDir).walkTopDown().filter { it.isFile }.filter { it.extension == "h.vcf" }.map { "${it.path}/${it.name}" }.toList()

    //create a HaplotypeGraph from the list of hvcf files
    val graph = HaplotypeGraph(hvcfFiles)

    val kmerIndexMap = importKmerIndex(kmerIndexFile)
}



/**
 * Loads a kmer hash map file to memory for use in mapping reads.  Initially based on PHGv1 source code.
 */
fun importKmerIndex(filename: String): KmerMapData {
    //Load the contents of the file into
    //rangeHapidMap: a map of refRangeId to an IntArray of the haplotype ids in the ReferenceRange
    // and the BitSet of all hapid sets
    //kmerHashmap: a map of kmer hash to reference range and offset into its BitSet, encoded as a long
    //These data structures all the reference range and haplotype set to be looked up for a kmer has

    val rangeToBitSetMap = mutableMapOf<Int, BitSet>()
    val kmerHashMap = Long2LongOpenHashMap()
    var lineCount = 0
    var totalLineCount = 0
    var refrangeId = 0
    var haplotypeListId = -1

    File(filename).bufferedReader().useLines {
        it.forEach { inputStr ->
            totalLineCount++

            lineCount = when(inputStr.first()) {
                '>' -> 1
                '#' -> 4
                else -> lineCount + 1
            }

            when (lineCount) {
                1 -> {
                    refrangeId = inputStr.substring(1).toInt()                  }
                2-> {
                    try {
                        val parsedLine = inputStr.split(",")
                        val myBitset = BitSet.valueOf(parsedLine.map { it.toLong() }.toLongArray())
                        rangeToBitSetMap[refrangeId] = myBitset
                    } catch(e: Exception) {
                        myLogger.info("error at line $totalLineCount for input = $inputStr")
                        throw java.lang.IllegalArgumentException(e)
                    }
                }
                3 -> {
                    val refrangeLong = (refrangeId.toLong()) shl 32
                    val parsedLine = inputStr.split(",")
                    //values are pairs of hash, offset (long,int)
                    //add the pairs of entries to kmerHashMap
                    for (datapair in parsedLine) {
                        val hashOffset = datapair.split("@")

                        if (hashOffset.size < 2) {
                            myLogger.info("improperly formatted datapair at line $totalLineCount")
                        }

                        val hash = hashOffset[0].toLong()
                        val offset = refrangeLong or hashOffset[1].toLong()
                        kmerHashMap.put(hash, offset)
                    }
                }
                4 -> {
                    if (inputStr.startsWith("#haplotypeListId="))
                        haplotypeListId = inputStr.substringAfter('=', "-1").toInt()

                }
            }
        }
    }

    return KmerMapData(haplotypeListId, rangeToBitSetMap, kmerHashMap)
}