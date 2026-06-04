package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.getBufferedReader

class Ps4gFileReader(val filename: String) {
    data class Ps4gGameteSet(val gameteIndices: IntArray, val count: Int)

    private val gameteIndexMap = mutableMapOf<Int, String>()
    private val contigToDataMap = mutableMapOf<String, MutableMap<Int, MutableList<Ps4gGameteSet>>>()

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
            while(!currentLine.startsWith("#gamete") && currentLine != null) currentLine = bufferedReader.readLine()
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
                val refPosition = parsedLine[2].toInt()
                val count = parsedLine[3].toInt()
                val contigMap = contigToDataMap.getOrPut(refContig) { mutableMapOf() }
                val dataList = contigMap.getOrPut(refPosition) {mutableListOf()}
                dataList.add(Ps4gGameteSet(gameteArray, count))

                currentLine = bufferedReader.readLine()
            }

        }

    }

    fun readMapForContig(contig: String) : Map<Int, MutableList<Ps4gGameteSet>>? {
        return contigToDataMap[contig]?.toMap()
    }

    fun contigSet() : Set<String> {
        return contigToDataMap.keys.toSet()
    }

    fun gameteIndexMap() : Map<Int, String> {
        return gameteIndexMap.toMap()
    }
}