package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import org.apache.logging.log4j.LogManager
import java.io.File
import java.lang.Exception
import java.nio.file.Files
import java.util.stream.Collectors

/**
 * This class uses the Assembled Genome Compressor (AGC) program to create
 * a compressed representation of the input assembly genome files.  The AGC
 * documentation is here:
 *    https://github.com/refresh-bio/agc
 *
 * This program will create a single compressed file from the input fasta files.
 * It is expected to live in the same parent folder as the TileDB datasets.
 *
 * The sample name for each fasta is determined by the file name, minus extension.  Users
 * should consider this when naming their fasta files.
 *
 * If a compressed file named "assemblies.agc" already exists in the parent folder,
 * any fasta in the fasta folder not currently in the agc compressed file will be added.
 * A list of duplicate fasta files, those not added, will be printed to the log.
 * Duplicates are determined based on the fastas name
 *
 * When appending to an agc file:  the append command takes both an input and an output agc file.
 * It does not actually append, it makes a copy of the old file and adds everything to it.
 */
class AgcCompress : CliktCommand(){

    private val myLogger = LogManager.getLogger(AgcCompress::class.java)

    val dbPath by option(help = "Folder name where AGC compressed files will be created.  Should be the same parent folder where tiledb datasets are stored.")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--db-path must not be blank"
            }
        }

    val fastaList by option(help = "File containing full path name fasta files, one per line, to compress into a single agc file.  Fastas may be compressed or uncompressed files. Reference fasta should NOT be included.")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--fasta-list must not be blank"
            }
        }

    val refFasta by option(help = "Full path to the reference fasta file to be added with other fastas to the agc compressed file. Reference fasta should NOT be compressed.")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--ref-fasta must not be blank"
            }
        }


    override fun run() {
        myLogger.info("Starting AGC compression")
        // process the input
        processAGCFiles(dbPath,fastaList,refFasta)

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



    fun processAGCFiles(dbPath:String, fastaList:String, refFasta:String) {

        val tempDir = "${dbPath}/temp"
        File(tempDir).mkdirs()

        // read the fastaList file into a list of strings
        val fastaFiles = File(fastaList).readLines()

        // check if ${dbPath}/assemblies.agc file exist
        if (Files.exists(File("${dbPath}/assemblies.agc").toPath())) {
            // if it exists, check if the fasta files in the fastaList/fastaFiles are already in the agc file
            // Print a list of the duplicates, add the non-duplicates to the compressed file.
            val duplicateList =
                compareNewExistingSampleNames("${dbPath}/assemblies.agc", getCurrentSampleNames(fastaFiles),tempDir)
            if (duplicateList.isNotEmpty()) {
                println("The following fasta files are already represented in the AGC compressed file and will not be loaded: ${duplicateList}")
            }
            val listToLoad = fastaFiles.filter { !duplicateList.contains(it) }
            // Write the listToLoad to a file named tempDir/nonDuplicateFastaFiles.txt
            val fileToLoad = File("${tempDir}/nonDuplicateFastaFiles.txt")
            fileToLoad.writeText(listToLoad.joinToString("\n"))

            if (listToLoad.isNotEmpty()) {
                // Call method to load AGC files with the list of fasta files and load option
                val success = loadAGCFiles(fileToLoad.toString(), "append", dbPath, refFasta,tempDir)
            } else {
                myLogger.info("No new fasta files to load -returning")
            }
        } else {
            // call method to load AGC files with the list of fasta files and load option
            val success = loadAGCFiles(fastaList, "create",dbPath,refFasta, tempDir)
        }

    }

    fun loadAGCFiles(fastaFiles: String, loadOption: String, dbPath:String, refFasta:String, tempDir:String): Boolean {

        val agcFile = "${dbPath}/assemblies.agc"
        val agcFileOut = "${dbPath}/assemblies_tmp.agc"
        val builder = ProcessBuilder()
        var redirectOutput = tempDir + "/agc_create_output.log"
        var redirectError = tempDir + "/agc_create_error.log"
        when (loadOption) {
            "create" -> {
                builder.command("conda","run","-n","phgv2-conda","agc","create","-i",fastaFiles,"-o",agcFile,refFasta)
            }
            "append" -> {
                builder.command("conda","run","-n","phgv2-conda","agc","append","-i",fastaFiles,agcFile,"-o",agcFileOut)
                redirectOutput = tempDir + "/agc_append_output.log"
                redirectError = tempDir + "/agc_append_error.log"
            }
            else -> {
                println("Invalid load option: ${loadOption}")
                return false
            }
        }

        builder.redirectOutput( File(redirectOutput))
        builder.redirectError( File(redirectError))

        println("begin Command to create/append:" + builder.command().stream().collect(Collectors.joining(" ")))
        try {
            var process = builder.start()
            var error = process.waitFor()
            if (error != 0) {
                println("agc create run via ProcessBuilder returned error code $error")
                throw IllegalArgumentException("Error running ProcessBuilder to compress agc files: ${error}")
            }

        } catch (exc: Exception) {
            println("Error: could not create agc compressed file.")
            throw IllegalArgumentException("Error running ProcessBuilder for agc create or append: ${exc.message}")
        }

        // TODO: if command was append, copy the agcFile to agcFile_backup_<date>.agc and then copy agcFileOut to agcFile


        return true
    }
    fun getSampleListFromAGC(agcFile:String,tempDir:String): List<String> {
        // This function will return a list of samples from the AGC compressed file.
        // This will be used to verify that the new list of fastas has nothing overlapping
        // the exsiting fastas in the AGC compressed file.

        var sampleList = listOf<String>()
        // Query the agc compressed file to get list of sample names
        var builder = ProcessBuilder("conda","run","-n","phgv2-conda","agc","listset",agcFile)
        var redirectOutput = tempDir + "/agc_create_output.log"
        var redirectError = tempDir + "/agc_create_error.log"
        builder.redirectOutput( File(redirectOutput))
        builder.redirectError( File(redirectError))

        myLogger.info("begin Command to list agc sample names:" + builder.command().stream().collect(Collectors.joining(" ")))
        try {
            var process = builder.start()
            var error = process.waitFor()
            if (error != 0) {
                println("agc listset run via ProcessBuilder returned error code $error")
                throw IllegalArgumentException("Error running ProcessBuilder to list agc samples: ${error}")
            }
            // read the output file to get the list of samples
            sampleList = File(redirectOutput).readLines()

        } catch (exc: Exception) {
            println("Error: could not list agc samples.")
            throw IllegalArgumentException("Error running ProcessBuilder for agc listset: ${exc.message}")
        }

        return sampleList
    }

    // This function will take a list of fasta names compiled from the input fasta
    // files and compare them to those name already in the agc compressed file.
    // It returns a list of any duplicates
    fun compareNewExistingSampleNames(agcFile:String, newFastas:List<String>, tempDir:String): List<String> {
        // Need to process the newFastas list to just the name without extension or path
        // that is what is stored as the sample name in the AGC compressed file.
        val newSampleNames = mutableListOf<String>()
        newFastas.forEach { newSampleNames.add(it.split("/").last().split(".").first()) }

        // Query the agc compressed file to get list of sample names
        val duplicateList = getSampleListFromAGC(agcFile,tempDir).intersect(newSampleNames).toList()

        // Match the duplicateList to the newFastas list to get the full path and extension
        // of the fasta files that are duplicates
        val duplicateFastas = newFastas.filter { duplicateList.contains(it.split("/").last().split(".").first()) }

        return duplicateFastas
    }
}