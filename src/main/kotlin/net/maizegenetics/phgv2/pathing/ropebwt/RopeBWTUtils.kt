package net.maizegenetics.phgv2.pathing.ropebwt


/**
 * Some classes to hold the data from the ropebwt3 mem output
 * MEMs here are Maximal Exact Matches.  These are the regions of the read that have a full similarity to the reference.
 * MEMHits are specific contig, strand and start positions of the MEM result.
 */
data class MEM(val readName: String, val readStart: Int, val readEnd: Int, val numHits: Int, val listMemHits: List<MEMHit>)
data class MEMHit(val contig: String, val strand: String, val pos: Int)


class RopeBWTUtils {
    companion object {
        /**
         * Function to parse the current alignment line from ropebwt3 mem into a usable object
         */
        fun parseStringIntoMem(string: String) : MEM {
            val split = string.split("\t")
            val readName = split[0]
            val readStart = split[1].toInt()
            val readEnd = split[2].toInt()
            val numHits = split[3].toInt()
            val listMemHits = split.subList(5, split.size).filter { it.isNotEmpty() }.map {
                val hitSplit = it.split(":")
                MEMHit(hitSplit[0], hitSplit[1], hitSplit[2].toInt())
            }
            return MEM(readName, readStart, readEnd, numHits, listMemHits)
        }

    }
}