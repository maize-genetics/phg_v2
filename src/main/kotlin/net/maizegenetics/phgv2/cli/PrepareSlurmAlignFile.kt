package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.int
import java.io.File

/**
 * This is a convenience class for creating a file that can be submitted to slurm as a data array job.
 * It is not required to be run - users can create the file manually if they prefer.
 *
 * This class provides a means of preparing data for aligning assemblies in a slurm data array job.
 * It takes as input the values for parameters needed to run the align-assemblies command.  This
 * could be done manually, but if the list of assemblies to align is long, this class provides
 * a quick way to create the file needed for a slurm data array job.
 *
 * Output:  A file that contains a line for each assembly in the list of assemblies to align.  Each line
 * is a call to the phg align-assemblies command.  The parameters for the call are the parameters supplied
 * to this class.
 *
 * The output file can be submitted to slurm as a data array job.
 */

class PrepareSlurmAlignFile: CliktCommand(help = "create files for aligning assemblies in a slurm data array job") {

    val phgLocation by option(help = " full path to the phg executable to be used in the slurm command file. " +
            "Default assumes command is run from the folder containing the executable.  This will be the first" +
            " part of the command in the slurm command file. ")
        .default("./phg")

    val gff by option(help = "Full path to the reference gff file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--gff must not be blank"
            }
        }

    val referenceFile by option(help = "Full path to reference file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--reference-file must not be blank"
            }
        }

    val referenceSam by option(help = "Full path to reference SAM file created by align-assemblies command when the just-ref-prep option is used." +
    "This is the SAM file created when aligning the reference fasta to the reference CDS regions.")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--reference-sam must not be blank"
            }
        }

    val referenceCdsFasta by option(help = "Full path to reference fasta file created by AlignAssemblies class when aligning to CDS regions")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--reference-cds-fasta must not be blank"
            }
        }


    // THis is a list of assemblies vs the option of a single file because the output
    // is a text file that contains a line for each assembly to be aligned.
    val assemblies by option(
        "-a",
        "--assemblies",
        help = "File containing list of assemblies to align, 1 per line, full path to updated file created via the phg prepare-assemblies command."
    )
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--assemblies must not be blank"
            }
        }

    val slurmCommandFile by option( help = "Full path and name for a file where the individual align-assembly commands will be printed")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--slurm-command-file must not be blank"
            }
        }

    val outputDir by option("-o", "--output-dir", help = "Directory where temporary and final files from anchorwave alignment will be written")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--output-dir must not be blank"
            }
        }

    val totalThreads by option(help = "Number of threads available.  These will be split among the alignments that are run in parallel")
        .int()
        .default(1)

    val refMaxAlignCov by option(help = "Anchorwave proali parameter R, indicating reference genome maximum alignment coverage.")
        .int()
        .default(1)

    val queryMaxAlignCov by option(help = "Anchorwave proali parameter Q, indicating query genome maximum alignment coverage.")
        .int()
        .default(1)

    val condaEnvPrefix by option (help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment." +
            " This is necessary if the conda env is not named phgv2-conda and/or is not in the conda default location")
        .default("")

    override fun run() {
        // This method creates the files needed for aligning assemblies in a slurm data array job.
        // For each assembly in the list, add an entry to the slurm data array job file that calls align-assemblies clikt command
        // THis is implemented in the createSlrumDatArrayJobFile() function
        println("PrepareSlurmALignFile: in run, calling createSlurmDataArrayJobFile")
        createSlurmDataArrayJobFile(
            phgLocation,
            gff,
            referenceFile,
            referenceSam,
            referenceCdsFasta,
            assemblies,
            slurmCommandFile,
            outputDir,
            totalThreads,
            1, // default for in-paralel as there is only one assembly per line
            refMaxAlignCov,
            queryMaxAlignCov,
            condaEnvPrefix
        )
    }

    // This method creates the files needed for aligning assemblies in a slurm data array job.
    // it creates a single file, where there is a line for each assembly to be aligned.
    // In each of these lines, a call is made to the phg align-assemblies command using the parameters supplied here
    // The output of this method is a file that can be submitted to slurm as a data array job.
    //The input to the function is all the parameters that are defined above.
    //The output is a file that can be submitted to slurm as a data array job.
    //The file will contain a line for each assembly to be aligned.  Each line will be a call to the phg align-assemblies command.
    //The parameters for the call will be the parameters supplied to this function.
    fun createSlurmDataArrayJobFile(
        phgLocation: String,
        gff: String,
        referenceFile: String,
        referenceSam: String,
        referenceCdsFasta: String,
        assemblies: String,
        slurmCommandFile: String,
        outputDir:String,
        totalThreads: Int,
        inParallel: Int,
        refMaxAlignCov: Int,
        queryMaxAlignCov: Int,
        condaEnvPrefix:String
    ) {
        // This method creates the file needed for aligning assemblies in a slurm data array job.
        // For each assembly in the list, add an entry to the slurm data array job file that calls align-assemblies clikt command

        // Create a list of assemblies
        val assembliesList = File(assemblies).readLines().filter { it.isNotBlank() }
        val writer = File(slurmCommandFile).bufferedWriter()
        // for each assembly in the list, write a line to the file that calls the align-assemblies command
        // We use the parameters supplied to this function as the parameters for the align-assemblies command,

        // set the conda-env-prefix if one was supplied
        val condaEnv = if (condaEnvPrefix.isNotBlank()) {
            "--conda-env-prefix $condaEnvPrefix"
        } else {
            ""
        }
        assembliesList.forEach {
            // THis might not be correct - check copilots' parameters.
            // TODO: Should I even have in-parallel and total-threads as parameters?
            writer.write("$phgLocation align-assemblies --gff $gff --output-dir $outputDir --reference-file $referenceFile --reference-sam $referenceSam --reference-cds-fasta $referenceCdsFasta --assembly-file $it --total-threads $totalThreads --in-parallel $inParallel --ref-max-align-cov $refMaxAlignCov --query-max-align-cov $queryMaxAlignCov ${condaEnv}\n")
        }
        writer.close()
    }
}
