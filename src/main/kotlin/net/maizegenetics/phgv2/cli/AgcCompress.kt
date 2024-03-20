package net.maizegenetics.phgv2.cli

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import net.maizegenetics.phgv2.utils.verifyURI
import org.apache.logging.log4j.LogManager
import java.io.File
import java.lang.Exception
import java.nio.file.Files
import java.nio.file.Paths
import java.time.LocalDate

/**
 * This class uses the Assembled Genome Compressor (AGC) program to create
 * a compressed representation of the input assembly genome files.  The AGC
 * documentation is here:
 *    https://github.com/refresh-bio/agc
 *
 * This program will create a single compressed file from the input fasta files.
 * It is expected to live in the same parent folder as the TileDB datasets.
 *
 * The input fasta files should be the files created from the AnnotateFasta command.
 * AGC does not keep track of the sample names when it pulls sequence from multiple fasta files.
 * For the phg software to identify the sample from which a sequence came, we annotated the fasta
 * id line with the sample name.  This is done in the AnnotateFasta command.  This AgcCompress
 * command requires the fasta files to be annotated.
 *
 * The sample name for each fasta is determined by the file name, minus extension.  Users
 * should consider this when naming their fasta files.
 *
 * If a compressed file named "assemblies.agc" already exists in the parent folder,
 * any fasta in the fasta list not currently in the agc compressed file will be added.
 * This is verified based on the fasta name (minus extension), not the contents of the fasta.
 * A list of duplicate sample names (fasta names) will be printed to the log.  These
 * duplicates will not be added to the compressed file.
 *
 * When appending to an agc file:  the append command takes both an input and an output agc file.
 * It does not actually append, it makes a copy of the old file and adds everything to it.
 * This code will copy the old file to a backup file, then create a new file with the name
 * assemblies.agc.  The backup file will be named assemblies_backup_YYYY-MM-DD.agc, where
 * YYYY-MM-DD is the date the backup was created.
 *
 */
class AgcCompress : CliktCommand(help = "Create a single AGC compressed file from an input of FASTA files created with the phg annotate-fasta command") {

    private val myLogger = LogManager.getLogger(AgcCompress::class.java)

    val dbPath by option(help = "Folder name where AGC compressed files will be created.  Should be the same parent folder where tiledb datasets are stored." )
        .default("")


    val fastaList by option(help = "File containing full path name for the annotated fasta files, one per line, to compress into a single agc file.  Fastas may be compressed or uncompressed files. Reference fasta should NOT be included.\nAll fastas must be fastas created via the phg annotate-fastas command")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--fasta-list must not be blank"
            }
        }

    val referenceFile by option(help = "Full path to the reference fasta file to be added with other fastas to the agc compressed file. Reference fasta should NOT be compressed and must be the annotated fasta from the annotate-fastas command.")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--reference-file must not be blank"
            }
        }


    override fun run() {
        myLogger.info("Starting AGC compression: validate the URI")
        // If the dbPath is not provided, use the current working directory
        val tiledbFolder = if (dbPath.isBlank()) {
            System.getProperty("user.dir")
        } else {
            dbPath
        }
        // Verify the dbPath contains valid tiledb created datasets
        // If it doesn't an exception will be thrown
        val validDB = verifyURI(tiledbFolder,"hvcf_dataset")
        // process the input
        processAGCFiles(tiledbFolder,fastaList,referenceFile)

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
                myLogger.info("The following fasta files are already represented in the AGC compressed file and will not be loaded: ${duplicateList}")
            }

            val listToLoad = fastaFiles.filter { !duplicateList.contains(it.substringAfterLast("/").substringBefore(".")) }
            // Write the listToLoad to a file named tempDir/nonDuplicateFastaFiles.txt
            val fileToLoad = File("${tempDir}/nonDuplicateFastaFiles.txt")
            fileToLoad.writeText(listToLoad.joinToString("\n"))

            if (listToLoad.isNotEmpty()) {
                // verify that the fasta files are annotated
                val startTime = System.nanoTime()
                verifyFileAnnotation(listToLoad) // this with throw an exception if the files are not annotated
                // print out time it took to verify the fasta files in seconds
                myLogger.info("VerifyFileAnnotation: time: " + (System.nanoTime() - startTime).toDouble() / 1_000_000_000.0 + " secs.")
                // Call method to load AGC files with the list of fasta files and load option
                myLogger.info("calling loadAGCFiles")
                val success = loadAGCFiles(fileToLoad.toString(), "append", dbPath, refFasta,tempDir)
            } else {
                myLogger.info("No new fasta files to load -returning")
            }
        } else {
            // call method to load AGC files with the list of fasta files and load option
            // verify that the fasta files are annotated
            val startTime = System.nanoTime()
            verifyFileAnnotation(fastaFiles) // this with throw an exception if the files are not annotated
            // print out time it took to verify the fasta files in seconds
            myLogger.info("VerifyFileAnnotation: time: " + (System.nanoTime() - startTime).toDouble() / 1_000_000_000.0 + " secs.")
            myLogger.info("calling loadAGCFiles")
            val success = loadAGCFiles(fastaList, "create",dbPath,refFasta, tempDir)
        }

    }

    // This function will verify that the fasta files are annotated
    // It looks for the first line of each fasta file that begins with ">"
    // and checks if that line contains "sampleName="
    // If it does, we assume the file is annotated
    // if it doesn't, we assume the file is not annotated
    // we are only checking the first idline of the fasta file, assuming all or none are annotated.
    // The code throws an exception on the first fasta file that is not annotated.
    fun verifyFileAnnotation(fastaFiles:List<String>): Boolean {
        // This function will verify that the fasta files are annotated
        // This is done by finding the first line of each fasta file that begins with ">"
        // and checking if that line contains "sampleName="
        // If it does not, the file is not annotated and an exception is thrown
        myLogger.info("Verifying fasta files are annotated")

        fastaFiles.forEach {
            // find the first line that begins with ">"
            bufferedReader(it).use { reader ->
                val line = reader.readLine()
                while (line != null) {
                    if (line.startsWith(">")) {
                        if (!line.contains("sampleName=")) {
                            myLogger.error("Fasta file ${it} is not annotated.  Please use the phg annotate-fastas command to annotate the fasta files.")
                            throw IllegalStateException("Fasta file ${it} is not annotated.  Please use the phg annotate-fastas command to annotate the fasta files.")
                        }
                        else {
                            // go to next fasta file
                            break
                        }
                    } else {
                        // read next line
                        reader.readLine()
                    }
                }
            }
        }

        return true
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
                myLogger.error("Invalid load option: ${loadOption}")
                return false
            }
        }

        builder.redirectOutput( File(redirectOutput))
        builder.redirectError( File(redirectError))

        myLogger.info("begin Command to create/append:" + builder.command().joinToString(" "))
        try {
            var process = builder.start()
            var error = process.waitFor()
            if (error != 0) {
                myLogger.error("agc create run via ProcessBuilder returned error code $error")
                throw IllegalArgumentException("Error running ProcessBuilder to compress agc files: ${error}")
            }
            // If this was an "append" command, the file was written to assemblies_tmp.agc
            // In this case, we need to move the original assemblies.agc file to assemblies_backup_<date>.agc
            // can then move the assemblies_tmp.agc to assemblies.agc

            if (loadOption == "append") {
                // Use Files.move instead of File.renameTo because File.renameTo does not work across file systems
                // and Files.move is platform independent
                val agcFileBackupFile = Paths.get("${dbPath}/assemblies_backup_${LocalDate.now()}.agc")
                val agcFileOutFile = Paths.get(agcFileOut)

                val agcFileOrig = Paths.get(agcFile)
                Files.move(agcFileOrig,agcFileBackupFile)
                Files.move(agcFileOutFile,agcFileOrig)
            }

        } catch (exc: Exception) {
            myLogger.error("Error: could not create agc compressed file.")
            throw IllegalArgumentException("Error running ProcessBuilder for agc create or append: ${exc.message}")
        }

        return true
    }

    // This function takes a list of fasta file names and returns a list of sample names
    // Sample names are the file names without extension
    fun getCurrentSampleNames(fastaFiles:List<String>):List<String> {
        //create a list of names from the fasta files, where the name is
        // just the file name, no path, and without the extension
        val sampleNames = mutableListOf<String>()
        fastaFiles.forEach{ sampleNames.add(it.substringAfterLast("/").substringBefore(".")) }

        return sampleNames
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

        myLogger.info("begin Command to list agc sample names:" + builder.command().joinToString(" "))
        try {
            var process = builder.start()
            var error = process.waitFor()
            if (error != 0) {
                myLogger.error("agc listset run via ProcessBuilder returned error code $error")
                throw IllegalArgumentException("Error running ProcessBuilder to list agc samples: ${error}")
            }
            // read the output file to get the list of samples
            sampleList = File(redirectOutput).readLines()

        } catch (exc: Exception) {
            myLogger.error("Error: could not list agc samples.")
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
        val newSampleNames = getCurrentSampleNames(newFastas)

        // Query the agc compressed file to get list of sample names
        val duplicateList = getSampleListFromAGC(agcFile,tempDir).intersect(newSampleNames).toList()

        // Match the duplicateList to the newFastas list to get the full path and extension
        // of the fasta files that are duplicates
        val duplicateFastas = newFastas.filter { duplicateList.contains(it.substringAfterLast("/").substringBefore(".")) }

        return duplicateFastas
    }
}