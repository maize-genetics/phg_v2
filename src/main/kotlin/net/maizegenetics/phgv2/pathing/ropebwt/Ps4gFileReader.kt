package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.utils.getBufferedReader
import org.apache.logging.log4j.LogManager

/**
 * Reads a PS4G (version 2.0) file, as written by [PS4GUtils.writeOutPS4GFile], into memory. A PS4G file records,
 * for a single sample, the set of gametes hit by each read mapping, along with the reference position of the hit.
 * The whole file is parsed by the constructor and held in memory, indexed by contig and binned reference position.
 * If supplied, [contigSet] filters out any contigs not in the list.
 *
 * The file has three parts: the two header lines "#PS4G" and "#version=2.0", a gamete index block that starts with
 * a "#gamete" line and maps each gamete name to the integer index used in the body, and the body itself, one line
 * per gamete set with the columns gameteSet, refContig, refPosBinned, and count. Any other "#" lines between the
 * version line and the gamete index block are metadata and are ignored.
 *
 * @param filename the path to the PS4G file to read.
 * @throws IllegalArgumentException if the file does not start with the "#PS4G" and "#version=2.0" lines, or if the
 * file ends before a gamete index block is found.
 */
class Ps4gFileReader(val filename: String, contigSet: Set<String> = setOf()) {
    /**
     * The gametes hit by a set of reads that all map to the same binned reference position.
     *
     * @param gameteIndices the indices of the gametes hit, as used in the PS4G file. [gameteIndexMap] maps those
     *                      indices back to gamete names.
     * @param count         the number of reads that produced this hit.
     */
    data class Ps4gGameteSet(val gameteIndices: IntArray, val count: Int) {
        override fun equals(other: Any?): Boolean {
            if (this === other) return true
            if (javaClass != other?.javaClass) return false

            other as Ps4gGameteSet

            if (count != other.count) return false
            if (!gameteIndices.contentEquals(other.gameteIndices)) return false

            return true
        }

        override fun hashCode(): Int {
            var result = count
            result = 31 * result + gameteIndices.contentHashCode()
            return result
        }
    }

    private val gameteIndexMap = mutableMapOf<Int, String>()
    private val contigToDataMap = mutableMapOf<String, MutableMap<Int, MutableList<Ps4gGameteSet>>>()
    private val filterContigs = contigSet.isNotEmpty()
    private val myLogger = LogManager.getLogger(Ps4gFileReader::class.java)

    //contig list
    //method to get read map by contig

    init {
        //is first line #PS4G ?
        //is second line #version=2.0 ?
        //read gamete index
        //read body into list
        //body could also be implemented as a Sequence

        getBufferedReader(filename).use { bufferedReader ->
            val firstRow = bufferedReader.readLine().trim()
            require(firstRow == "#PS4G") {"$filename is not a valid PS4G file. First line is $firstRow"}

            val secondRow = bufferedReader.readLine().trim()
            require(secondRow.equals("#version=2.0", ignoreCase = true)) {"Second row is not a valid PS4G version string."}

            var currentLine = bufferedReader.readLine()
            while(currentLine != null && !currentLine.startsWith("#gamete")) currentLine = bufferedReader.readLine()
            require(currentLine != null) {"End of ps4g file reached without finding gamete indices."}

            //populate the gamete index map
            currentLine = bufferedReader.readLine()
            while (currentLine.startsWith("#")) {
                val gameteInfo = currentLine.split("\t")
                gameteIndexMap[gameteInfo[1].toInt()] = gameteInfo[0].drop(1)
                currentLine = bufferedReader.readLine()
            }

            //read in the data
            currentLine = bufferedReader.readLine() //skip column header
            while (currentLine != null) {
                val parsedLine = currentLine.split("\t")
                val gameteArray = parsedLine[0].split(",").map { it.toInt() }.toIntArray()
                val refContig = parsedLine[1]
                if (!filterContigs || contigSet.contains(refContig)) {
                    val refPosition = parsedLine[2].toInt()
                    val count = parsedLine[3].toInt()
                    val contigMap = contigToDataMap.getOrPut(refContig) { mutableMapOf() }
                    val dataList = contigMap.getOrPut(refPosition) {mutableListOf()}
                    dataList.add(Ps4gGameteSet(gameteArray, count))
                }

                currentLine = bufferedReader.readLine()
            }

        }

        val numberOfContigs = contigToDataMap.size
        val numberOfPositions = contigToDataMap.map { (_, dataList) -> dataList.size }.sum()

        if (contigSet().isNotEmpty()) {
            myLogger.info("Ps4gFileReader read ${filename} using contig set: $contigSet.\n" +
                    "This resulted in a data set containing $numberOfContigs contigs and $numberOfPositions positions.")

        } else {
            myLogger.info("Ps4gFileReader read ${filename}.\nThis resulted in a data set containing " +
                    "$numberOfContigs contigs and $numberOfPositions positions.")
        }

    }

    /**
     * Returns the read mappings for [contig], or null if [contig] is not in the PS4G file.
     *
     * @return a map of binned reference position to the gamete sets recorded at that position. A position can hold
     * more than one gamete set because different reads mapping to the same bin can hit different gametes.
     */
    fun readMapForContig(contig: String) : Map<Int, MutableList<Ps4gGameteSet>>? {
        return contigToDataMap[contig]?.toMap()
    }

    /**
     * Returns the names of the contigs that have read mappings in the PS4G file.
     */
    fun contigSet() : Set<String> {
        return contigToDataMap.keys.toSet()
    }

    /**
     * Returns the gamete index block of the PS4G file as a map of gamete index to gamete name. The indices are the
     * ones used in [Ps4gGameteSet.gameteIndices].
     */
    fun gameteIndexMap() : Map<Int, String> {
        return gameteIndexMap.toMap()
    }
}