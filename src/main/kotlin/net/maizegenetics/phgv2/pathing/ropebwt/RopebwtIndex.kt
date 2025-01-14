package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.seqIO.NucSeqIO
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.cli.CreateFastaFromHvcf
import net.maizegenetics.phgv2.cli.logCommand
import org.apache.logging.log4j.LogManager
import java.io.File


/**
 * Class to hold the input for the ropeBWT3 index creation
 * This can give you a pangenome file or a hvcf dir
 */
sealed class RopeBWTIndexInput {
    data class PangenomeFile(val pangenomeFile: String) : RopeBWTIndexInput()
    data class HvcfDir(val hvcfDir: String) : RopeBWTIndexInput()
}


/**
 * Class to call the appropriate ropebwt3 commands to create a ropeBWT3 index
 * First it will run the build step to make the index
 * Then it converts the fmr index to fmd for speed and efficiency
 * Then it builds the suffix array for the index
 * Finally it builds the chromosome length file for the index
 *
 * With the output files we can then run ropebwt3 mem and get read mappings.
 */
class RopebwtIndex : CliktCommand(help="BETA: Create a ropeBWT3 index") {

    private val myLogger = LogManager.getLogger(RopebwtIndex::class.java)

    // Pre-compile the Regex pattern - used when creating the output fasta file names
    val HVCF_PATTERN = Regex("""(\.hvcf|\.h\.vcf|\.hvcf\.gz|\.h\.vcf\.gz)$""")

    val ropeBWTIndexInput : RopeBWTIndexInput? by mutuallyExclusiveOptions<RopeBWTIndexInput>(
        option(
            "--pangenome-file",
            help = "Path to pangenome fasta file created by create-fasta-from-hvcf."
        ).convert { RopeBWTIndexInput.PangenomeFile(it) },
        option(
            "--hvcf-dir",
            help = "Path to directory holding hVCF files. This will build a pangenome fasta file from the hvcf files."
        ).convert { RopeBWTIndexInput.HvcfDir(it) }
    ).single().required()

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  " +
            "If not provided, the current working directory is used.  Note that if --pangenome-file is used this is not used as the sequence is already extracted into a fasta file")
        .default("")

    val outputDir by option(help = "Output Directory")
        .required()

    val indexFilePrefix by option(help = "Prefix of the ropebwt3 index file.  This will be added to the output directory and ropeBWT3 will make a number of output files.")
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

        val pangenomeFastaFile = if(ropeBWTIndexInput is RopeBWTIndexInput.HvcfDir) {
            myLogger.info("Creating pangenome fasta file from hvcf files in ${(ropeBWTIndexInput as RopeBWTIndexInput.HvcfDir).hvcfDir}")
            val hvcfDir = (ropeBWTIndexInput as RopeBWTIndexInput.HvcfDir).hvcfDir
            createPangenomeFasta(hvcfDir, outputDir, dbPath, condaEnvPrefix)
            "$outputDir/pangenome.fa"
        }
        else {
            (ropeBWTIndexInput as RopeBWTIndexInput.PangenomeFile).pangenomeFile
        }

        myLogger.info("Creating ropeBWT3 index for $pangenomeFastaFile")

        createInitialIndex(pangenomeFastaFile, "${outputDir}/${indexFilePrefix}", numThreads, deleteFmrIndex, condaEnvPrefix)

    }

    fun createPangenomeFasta(hvcfDir: String, outputDir: String, dbPath: String, condaEnvPrefix: String) {
        val hvcfFiles : Array<File> = File(hvcfDir).listFiles { file ->
            HVCF_PATTERN.containsMatchIn(file.name)
        } ?: throw Exception("No hvcf files found in $hvcfDir")
        CreateFastaFromHvcf().createPangenomeHaplotypeFile(outputDir, hvcfFiles, dbPath, condaEnvPrefix)
    }

    fun createInitialIndex(inputFasta:String, indexFilePrefix:String, numThreads: Int, deleteFMRIndex:Boolean, condaEnvPrefix:String) {
        //Build fmt file
        //time ../ropebwt3/ropebwt3 build -t24 -bo phg_ASMs.fmr /workdir/zrm22/phgv2/ropeBWT/fullASMTests/ASMs/B73.fa
        runBuildStep(inputFasta, indexFilePrefix, numThreads, condaEnvPrefix)

        //Convert the fmr to fmd to make it static and efficient
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
                writer.write("${seq.id}\t${seq.sequence.size()}\n")
            }
        }
    }

}