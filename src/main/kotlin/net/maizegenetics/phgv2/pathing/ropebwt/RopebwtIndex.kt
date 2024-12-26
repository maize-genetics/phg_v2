package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.seqIO.NucSeqIO
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.cli.logCommand
import org.apache.logging.log4j.LogManager
import java.io.File

class RopebwtIndex : CliktCommand(help="Create a ropeBWT3 index") {

    private val myLogger = LogManager.getLogger(RopebwtIndex::class.java)

    val inputFasta by option(help = "Pangenome Fasta File")
        .required() //Needs to be required now due to the agc archive

    val indexFilePrefix by option(help = "The full path of the ropebwt3 index file.")
        .required()

    val numThreads by option(help = "Number of threads to use for the index creation.")
        .int()
        .default(3)

    val deleteFmrIndex by option(help = "Delete the fmr index file after conversion to fmd.")
        .flag(default = true)

    val condaEnvPrefix by option (help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
        .default("")

    override fun run() {
        logCommand(this)

        myLogger.info("Creating ropeBWT3 index for $inputFasta")

        createInitialIndex(inputFasta, indexFilePrefix, numThreads, deleteFmrIndex, condaEnvPrefix)

    }

    fun createInitialIndex(inputFasta:String, indexFilePrefix:String, numThreads: Int, deleteFMRIndex:Boolean, condaEnvPrefix:String) {
        //Build fmt file
        //time ../ropebwt3/ropebwt3 build -t24 -bo phg_ASMs.fmr /workdir/zrm22/phgv2/ropeBWT/fullASMTests/ASMs/B73.fa
        runBuildStep(inputFasta, indexFilePrefix, numThreads, condaEnvPrefix)

        //Convert the fmr to fmt
        //time ../ropebwt3/ropebwt3 build -i /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmr -do /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmd
        convertBWTIndex(indexFilePrefix, condaEnvPrefix)

        //Delete the FMR if requested
        if (deleteFMRIndex) {
            deleteFMRIndex(indexFilePrefix)
        }

        //Build suffix array
        //time ../ropebwt3/ropebwt3 ssa -o /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmd.ssa -s8 -t32 /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmd
        buildSuffixArray(indexFilePrefix, numThreads, condaEnvPrefix)

        //Build chromosome seq lengths
        //cat /workdir/zrm22/phgv2/maize2_1/minimap2Tests2/pangenome.fa | /programs/seqtk/seqtk comp | cut -f1,2 | gzip > /workdir/zrm22/phgv2/ropeBWT/phg_generalSingleFile_out.fmd.len.gz
        buildChrLengthFile(inputFasta, indexFilePrefix)
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
            Pair("-n", "phgv2-conda") }
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
     * Function to convert the BWT index file to a format that is static but more efficient
     */
    fun convertBWTIndex(indexFilePrefix: String, condaEnvPrefix: String) {
        //Convert the fmr to fmt
        //time ../ropebwt3/ropebwt3 build -i /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmr -do /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmd
        val prefixArg = if(condaEnvPrefix.isNotBlank()) {
            Pair("-p",condaEnvPrefix)
        }
        else {
            Pair("-n", "phgv2-conda") }
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
            Pair("-n", "phgv2-conda") }
        val ssaCommand = listOf("conda","run",prefixArg.first,prefixArg.second,"ropebwt3", "ssa", "-o", "$indexFilePrefix.bwt.ssa", "-s8", "-t${numThreads}", "$indexFilePrefix.fmd")
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

    /**
     * Function to build the chromosome length file for the BWT index
     * This is needed to get the chromosome lengths for the pangenome
     * We may need to remove this with a biokotlin based function in the future as cat does some odd stuff if the EOF is missing
     */
//    fun buildChrLengthFile(inputFasta: String, indexFilePrefix: String, condaEnvPrefix: String) {
//        //Build chromosome seq lengths
//        //cat /workdir/zrm22/phgv2/maize2_1/minimap2Tests2/pangenome.fa | /programs/seqtk/seqtk comp | cut -f1,2 | gzip > /workdir/zrm22/phgv2/ropeBWT/phg_generalSingleFile_out.fmd.len.gz
//        val prefixArg = if(condaEnvPrefix.isNotBlank()) { "-p" } else { "-n" }
//        val chrLengthCommand = listOf("conda","run",prefixArg,"phgv2-conda","cat", inputFasta, "|", "seqtk", "comp", "|", "cut", "-f1,2", "|", "gzip", ">", "$indexFilePrefix.fmd.len.gz")
//        myLogger.info("Running chromosome length command: ${chrLengthCommand.joinToString(" ")}")
//        try {
//            val process = ProcessBuilder(chrLengthCommand)
//                .inheritIO()
//                .start()
//            process.waitFor()
//        } catch (e: Exception) {
//            myLogger.error("Error running chromosome length command: ${chrLengthCommand.joinToString(" ")}")
//            throw e
//        }
//    }

    /**
     * Function to build the chromosome length file for the BWT index
     * This is needed to get the chromosome lengths for the pangenome
     */
    fun buildChrLengthFile(inputFasta: String, indexFilePrefix: String) {
        //load the input fasta into nucSeq
        //iterate over the nucSeq and find the lengths of each chromosome
        //write the lengths to a file
        bufferedWriter("$indexFilePrefix.fmd.len.gz").use { writer ->
            //write the lengths to the file
            val neqSeqIO = NucSeqIO(inputFasta)
            neqSeqIO.forEach { seq ->
                writer.write("${seq.name}\t${seq.sequence.size()}\n")
            }
        }
    }

}