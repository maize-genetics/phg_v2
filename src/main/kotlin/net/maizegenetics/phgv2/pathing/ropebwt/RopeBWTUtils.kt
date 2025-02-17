package net.maizegenetics.phgv2.pathing.ropebwt

import org.apache.logging.log4j.LogManager
import java.io.File


/**
 * Some classes to hold the data from the ropebwt3 mem output
 * MEMs here are Maximal Exact Matches.  These are the regions of the read that have a full similarity to the reference.
 * MEMHits are specific contig, strand and start positions of the MEM result.
 */
data class MEM(val readName: String, val readStart: Int, val readEnd: Int, val numHits: Int, val listMemHits: List<MEMHit>)
data class MEMHit(val contig: String, val strand: String, val pos: Int)


class RopeBWTUtils {
    companion object {
        val myLogger = LogManager.getLogger(RopeBWTUtils::class.java)
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

        /**
         * Function to run the ropebwt3 build command
         * This builds the initial BWT index file in an update-able format
         */
        fun runBuildStep(inputFasta:String, indexFilePrefix:String, numThreads: Int, condaEnvPrefix:String) {
            //"conda","run","-p","phgv2-conda"
            val prefixArg = if(condaEnvPrefix.isNotBlank()) {
                Pair("-p",condaEnvPrefix)
            }
            else {
                Pair("-n", "phgv2-ropebwt-conda") }
            val buildCommand = listOf("conda","run",prefixArg.first,prefixArg.second,"ropebwt3", "build", "-t$numThreads", "-bo", "$indexFilePrefix.fmr", "$inputFasta")
            myLogger.info("Running ropebwt3 build command: ${buildCommand.joinToString(" ")}")
            try {
                val process = ProcessBuilder(buildCommand)
                    .inheritIO()
                    .start()
                process.waitFor()
            } catch (e: Exception) {
                myLogger.error("Error running ropebwt3 build command: ${buildCommand.joinToString(" ")}")
                throw e
            }
        }

        /**
         * Function to run the ropebwt3 build Command but writing to a temp file and then deleting it.
         *
         */
        fun runBuildUpdateStep(inputFasta:String, indexFilePrefix:String, numThreads: Int, condaEnvPrefix:String) {
            val tempIndex = "${indexFilePrefix}_temp.fmr"
            val prefixArg = if(condaEnvPrefix.isNotBlank()) {
                Pair("-p",condaEnvPrefix)
            }
            else {
                Pair("-n", "phgv2-ropebwt-conda")
            }
            //time ../ropebwt3/ropebwt3 build -t24 -i phg_ASMs_input.fmr -bo phg_ASMs.fmr /workdir/zrm22/phgv2/ropeBWT/fullASMTests/ASMs/A188
            val buildCommand = listOf("conda","run",prefixArg.first,prefixArg.second,"ropebwt3", "build", "-t$numThreads", "-i", "$indexFilePrefix.fmr", "-bo", tempIndex, "$inputFasta")
            myLogger.info("Running ropebwt3 build command: ${buildCommand.joinToString(" ")}")
            try {
                val process = ProcessBuilder(buildCommand)
                    .inheritIO()
                    .start()
                process.waitFor()
            } catch (e: Exception) {
                myLogger.error("Error running ropebwt3 build command: ${buildCommand.joinToString(" ")}")
                throw e
            }
            File(tempIndex).renameTo(File("$indexFilePrefix.fmr"))
        }

        /**
         * Function to convert the BWT index file to a format that is static but more efficient
         * This is recommended by the ropebwt3 documentation
         */
        fun convertBWTIndex(indexFilePrefix: String, condaEnvPrefix: String) {
            //Convert the fmr to fmt
            //time ../ropebwt3/ropebwt3 build -i /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmr -do /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmd
            val prefixArg = if(condaEnvPrefix.isNotBlank()) {
                Pair("-p",condaEnvPrefix)
            }
            else {
                Pair("-n", "phgv2-ropebwt-conda") }
            val convertCommand = listOf("conda","run",prefixArg.first,prefixArg.second,"ropebwt3", "build", "-i", "$indexFilePrefix.fmr", "-do", "$indexFilePrefix.fmd")
            myLogger.info("Running ropebwt3 convert command: ${convertCommand.joinToString(" ")}")
            try {
                val process = ProcessBuilder(convertCommand)
                    .inheritIO()
                    .start()
                process.waitFor()
            } catch (e: Exception) {
                myLogger.error("Error running ropebwt3 convert command: ${convertCommand.joinToString(" ")}")
                throw e
            }
        }


        fun deleteFMRIndex(indexFilePrefix: String) {
            try {
                File("$indexFilePrefix.fmr").delete()
            }
            catch (e: Exception) {
                myLogger.error("Error deleting $indexFilePrefix.fmr")
                throw e
            }
        }

        /**
         * Function to build the suffix array for the BWT index
         */
        fun buildSuffixArray(indexFilePrefix: String, numThreads: Int, condaEnvPrefix: String) {
            //Build suffix array
            //time ../ropebwt3/ropebwt3 ssa -o /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmd.ssa -s8 -t32 /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmd
            val prefixArg = if(condaEnvPrefix.isNotBlank()) {
                Pair("-p",condaEnvPrefix)
            }
            else {
                Pair("-n", "phgv2-ropebwt-conda") }
            val ssaCommand = listOf("conda","run",prefixArg.first,prefixArg.second,"ropebwt3", "ssa", "-o", "$indexFilePrefix.fmd.ssa", "-s8", "-t${numThreads}", "$indexFilePrefix.fmd")
            myLogger.info("Running ropebwt3 ssa command: ${ssaCommand.joinToString(" ")}")
            try {
                val process = ProcessBuilder(ssaCommand)
                    .inheritIO()
                    .start()
                process.waitFor()
            } catch (e: Exception) {
                myLogger.error("Error running ropebwt3 ssa command: ${ssaCommand.joinToString(" ")}")
                throw e
            }
        }

    }
}