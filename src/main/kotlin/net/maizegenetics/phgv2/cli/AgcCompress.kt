package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import org.apache.logging.log4j.LogManager
import java.io.File
import java.nio.file.Files

/**
 * This class uses the Assembled Genome Compressor (AGC) program to create
 * a compressed representation of the input assembly genome files.  The AGC
 * documentation is here:
 *    https://github.com/refresh-bio/agc
 *
 * This program will create a single compressed file from the input fasta files.
 * It is expected to live in the same parent folder as the TileDB datasets.
 *
 * If a compressed file named "assemblies.agc" already exists, the fastas in the
 * folder will be added to the existing file.
 *
 * What checks should this contain to verify the list of fastas don't contain duplicates
 * to what already exists in the folder?
 */
class AgcCompress : CliktCommand(){

    private val myLogger = LogManager.getLogger(AgcCompress::class.java)

    val dbPath by option(help = "Folder name where AGC compressed files and tiledb datasets will be created")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--dbpath must not be blank"
            }
        }

    val fastaDir by option(help = "Folder containg fasta file to compress into a single file. Reference fasta should be uncompressed, other genome fastas may be compressed or uncompressed files.")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--fasta-dir must not be blank"
            }
        }

    override fun run() {

        // process the input
        createAppendToAGC(dbPath,fastaDir)

    }

    // This function takes a list of fasta file names and returns a list of sample names
    // Sample names are the file names without extension
    fun getCurrentSampleNames(fastaFiles:List<String>):List<String> {
        //create a list of names from the fasta files, where the name is
        // just the file name, no path, and without the extension
        val sampleNames = mutableListOf<String>()
        fastaFiles.forEach { sampleNames.add(it.split("/").last().split(".").first()) }
        return sampleNames
    }

    fun createAppendToAGC(dbPath:String,fastaDir:String) {
        // create a list of fasta files in the fastDir.  The list should be of only regular files
        // that end with .fa or .fasta, and the list should be the full path to the file, but each as a String
        val fastaFiles = File(fastaDir).walk().filter { it.isFile }.filter{it.name.endsWith(".fa") || it.name.endsWith(".fasta")}.map({it.toString()}).toList()

        //TODO: below is created by co-pilot - fix it up to make the agc process Builder calls!
        // check if ${dbPath}/assemblies.agc file exist
        if (Files.exists(File("${dbPath}/assemblies.agc").toPath())) {
            // if it exists, check if the fasta files in the fastaDir are already in the agc file
            // if they are, throw an error
            // if they are not, append the fasta files to the agc file
            val duplicateList =
                compareNewExistingSampleNames("${dbPath}/assemblies.agc", getCurrentSampleNames(fastaFiles))
            if (duplicateList.isNotEmpty()) {
                println("The following fasta files are already in the AGC compressed file: ${duplicateList}")
                println("Please remove these fasta files from the input fasta directory and try again.")
                System.exit(1)
            } else {
                // append the fasta files to the agc file
                // this is done by running the agc program with the -a option
                // the -a option is used to append to an existing agc file
                // the -o option is used to specify the output file name
                // the -i option is used to specify the input fasta files
                // the -t option is used to specify the number of threads to use
                // the -v option is used to specify the verbosity
            }
        } else {
            // if it does not exist, create the agc file
            // this is done by running the agc program with the -c option
            // the -c option is used to create a new agc file
            // the -o option is used to specify the output file name
            // the -i option is used to specify the input fasta files
            // the -t option is used to specify the number of threads to use
            // the -v option is used to specify the verbosity
        }

    }

    fun getSampleListFromAGC(agcFile:String): List<String> {
        // This function will return a list of samples from the AGC compressed file.
        // This will be used to verify that the new list of fastas has nothing overlapping
        // the exsiting fastas in the AGC compressed file.

        val sampleList = mutableListOf<String>()
        // Query the agc compressed file to get list of sample names

        // TODO:  query agc file to get list of sample names
        return sampleList
    }

    // THis function will take a list of fasta names compiled from the input fasta
    // files and compare them to those name already in the agc compressed file.
    // It returns a list of any duplicates
    fun compareNewExistingSampleNames(agcFile:String, newFastas:List<String>): List<String> {
        val duplicateList = getSampleListFromAGC(agcFile).intersect(newFastas).toList()

        return duplicateList
    }
}