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
class RopeBwtIndex : CliktCommand(help="BETA: Create a ropeBWT3 index") {

    private val myLogger = LogManager.getLogger(RopeBwtIndex::class.java)

    // Pre-compile the Regex pattern - used when creating the output fasta file names
    val HVCF_PATTERN = Regex("""(\.hvcf|\.h\.vcf|\.hvcf\.gz|\.h\.vcf\.gz)$""")

    val ropeBWTIndexInput : RopeBWTIndexInput by mutuallyExclusiveOptions(
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

    val threads by option(help = "Number of threads to use for the index creation.")
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

        createInitialIndex(pangenomeFastaFile, "${outputDir}/${indexFilePrefix}", threads, deleteFmrIndex, condaEnvPrefix)

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
        RopeBWTUtils.runBuildStep(inputFasta, indexFilePrefix, numThreads, condaEnvPrefix)

        //Convert the fmr to fmd to make it static and efficient
        //time ../ropebwt3/ropebwt3 build -i /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmr -do /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmd
        RopeBWTUtils.convertBWTIndex(indexFilePrefix, condaEnvPrefix)

        //Delete the FMR if requested
        if (deleteFMRIndex) {
            RopeBWTUtils.deleteFMRIndex(indexFilePrefix)
        }

        //Build suffix array
        //time ../ropebwt3/ropebwt3 ssa -o /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmd.ssa -s8 -t32 /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmd
        RopeBWTUtils.buildSuffixArray(indexFilePrefix, numThreads, condaEnvPrefix)

        //Build chromosome seq lengths
        //cat /workdir/zrm22/phgv2/maize2_1/minimap2Tests2/pangenome.fa | /programs/seqtk/seqtk comp | cut -f1,2 | gzip > /workdir/zrm22/phgv2/ropeBWT/phg_generalSingleFile_out.fmd.len.gz
        buildChrLengthFile(inputFasta, indexFilePrefix)
    }





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