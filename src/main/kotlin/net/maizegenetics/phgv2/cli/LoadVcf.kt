package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import org.apache.logging.log4j.LogManager
import java.io.File
import java.lang.Exception
import java.nio.file.Files
import java.nio.file.Paths
import java.util.stream.Collectors

/**
 * This class will load gvcf and hvcf files into TileDB datasets.
 * It has these steps:
 *  1. create list of gvcf vs hvcf files for loading
 *      a. if no hvcf.gz, h.vcf.gz, gvcf.gz, g.vcf.gz files are found, return error
 *  2. verify the dbpath exists and TileDB datasets exists in that folder
 *      a. if no datasets exist, create them
 *  3. load the gvcf and hvcf files
 */
class LoadVcf : CliktCommand(help="load g.vcf and h.vcf files into tiledb datasets") {

    private val myLogger = LogManager.getLogger(LoadVcf::class.java)

    val vcfDir by option(help = "VCF file directory")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--vcf-dir must not be blank"
            }
        }

    val dbPath by option(help = "Folder holding TileDB datasets")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--db-path must not be blank"
            }
        }

    // threads is not required.  Defaults to 1
    // While technically, this should be an "int", it is passed to ProcessBuilder as a string
    // so we're leaving it defined as a string here.
    val threads by option(help = "Number of threads for use by tiledb")
        .default("1")

    override fun run() {

        // Check the type of files in the vcfDir
        // antying with h.vcf.gz or hvcf.gz is a hvcf file
        // anything with g.vcf.gz or gvcf.gz is a gvcf file
        val fileLists = getFileLists(vcfDir)

        if (fileLists.first.isEmpty() && fileLists.second.isEmpty()) {
            // no gvcf or hvcf files found in the user supplied vcfDir
            myLogger.warn("No files ending in g.vcf.gz or h.vcf.gz found in $vcfDir.  Note that both the bgzipped and indexed files must exist in the specified folder. \nPlease check the folder and try again.")
            throw IllegalArgumentException("LoadVCF: No files ending in g.vcf.gz or h.vcf.gz found in $vcfDir.  Note that both the bgzipped and indexed files must exist in the specified folder. \nPlease check the folder and try again.")
        }

        // verify the URI is correct,
        if (fileLists.first.isNotEmpty()) {
            val gvcfExists = verifyURI(dbPath,"gvcf_dataset")
            if (gvcfExists){
                // Load the gvcf files!
                val datasetURI = "${dbPath}/gvcf_dataset"
                loadVCFToTiledb(fileLists.first, datasetURI, threads)
            }
        }
        if (fileLists.second.isNotEmpty()) {
            val hvcfExists = verifyURI(dbPath,"hvcf_dataset")
            if (hvcfExists){
                // Load the hvcf files!
                val datasetURI = "${dbPath}/hvcf_dataset"
                loadVCFToTiledb(fileLists.second, datasetURI, threads)
            }
        }

    }

    // This function walks the files in the vcf folder and returns a pair of lists,
    // one for gvcf files and one for hvcf files.
    fun getFileLists(vcfDir:String):Pair<List<String>,List<String>> {
        // get a list of files in the vcfDir
        // return a pair of lists, one for gvcf files and one for hvcf files
        val vcfFiles = File(vcfDir).listFiles()
        val gvcfFiles = mutableListOf<String>()
        val hvcfFiles = mutableListOf<String>()
        for (file in vcfFiles) {
            if (file.name.endsWith("g.vcf.gz") || file.name.endsWith("gvcf.gz")) {
                gvcfFiles.add(file.toString())
            } else if (file.name.endsWith("h.vcf.gz") || file.name.endsWith("hvcf.gz")) {
                hvcfFiles.add(file.toString())
            }
        }
        return Pair(gvcfFiles,hvcfFiles)
    }

    // The uri should either be gvcf_dataset or hvcf_dataset
    // The user determines the parent folder name where these datasets live
    // The actual tiledb dataset names are constant and are either gvcf_dataset or hvcf_dataset

    // We are only verifying, not creating the datasets
    fun verifyURI(dbPath:String,uri:String): Boolean {
        // Check that the user supplied db folder exists
        check(File(dbPath).exists()) { "Folder $dbPath does not exist - please send a valid path that indicates the parent folder for your tiledb datasets." }

        // Check if the dataset exists
        val dataset = "${dbPath}/${uri}"
        val datasetPath = Paths.get(dataset)

        if (File(dataset).exists() && Files.isRegularFile(datasetPath)) {
            throw IllegalArgumentException("URI ${dataset}is a file, not a tiledb dataset folder.  The parent folder must not contain any files/folders named gvcf_dataset or hvcf_dataset that is not a tiledb created URI")
        }

        // Create tne temp folder if it doesn't exist
        // This will be used to write output files from ProcessBuilder commands
        // called elsewhere in this class
        val tempDir = "${dbPath}/temp"
        File(tempDir).mkdirs()

        if (File(dataset).exists()  && Files.isDirectory(Paths.get(dataset))){
            // check if is a tiledb dataset
            var builder = ProcessBuilder("conda","run","-n","phgv2-conda","tiledbvcf","stat","--uri",dataset)
            var redirectOutput = tempDir + "/tiledb_statURI_output.log"
            var redirectError = tempDir + "/tiledb_statURI_error.log"
            builder.redirectOutput( File(redirectOutput))
            builder.redirectError( File(redirectError))

            // verify if the output.log contains "Version"
            // if not, then the URI is not a tiledbvcf URI
            myLogger.info("begin Command:" + builder.command().stream().collect(Collectors.joining(" ")))

            try {
                var process = builder.start()
                var error = process.waitFor()
                if (error != 0) {
                    myLogger.error("LoadTiledbH tiledb stat returned error code $error")
                    throw IllegalArgumentException("Error: URI is not a tiledb URI folder created via the tiledb create command: ${error}")
                }

            } catch (exc: Exception) {
                myLogger.error("Error: could not run tiledb stat on ${uri}.")
                throw IllegalArgumentException("Error running ProcessBuilder to stat tiledb URI: ${exc}")
            }

            myLogger.info("Using existing TileDB datasets previously created in folder $dbPath.")
            return true
        } else {
            myLogger.info("TileDB datasets not found in folder $dbPath. Please run InitDB to create the datasets.")
            return false
        }
    }

    fun loadVCFToTiledb(vcfList:List<String>, uri:String, threads:String) {

        // declare tne temp folder
        // This will be used to write output files from ProcessBuilder commands
        // called elsewhere in this class.  It should have been created in the
        // verifyURI function,
        val tempDir = "${dbPath}/temp"

        // get just last part of uri string, ie just the last name in this folder
        val uriName = uri.split("/").last()
        val vcfListFile = "${tempDir}/${uriName}_vcfList.txt"
        File(vcfListFile).writeText(vcfList.joinToString("\n"))
        // print the vcfListFile contents
        myLogger.info("vcfListFile contents: ${File(vcfListFile).readText()}")

        // Store the files to tiledb
        var builder = ProcessBuilder("conda","run","-n","phgv2-conda","tiledbvcf","store","--uri",uri,"-t",threads,"-f",vcfListFile,"--remove-sample-file")
        var redirectOutput = tempDir + "/tiledb_store_output.log"
        var redirectError = tempDir + "/tiledb_store_error.log"
        builder.redirectOutput( File(redirectOutput))
        builder.redirectError( File(redirectError))

        myLogger.info("begin Command:" + builder.command().stream().collect(Collectors.joining(" ")))
        println("begin Command:" + builder.command().stream().collect(Collectors.joining(" ")))
        try {
            var process = builder.start()
            var error = process.waitFor()
            if (error != 0) {
                myLogger.error("loadVCFToTiledb run via ProcessBuilder returned error code ${error}. See ${redirectError} and ${redirectOutput} for more details.")
                throw IllegalArgumentException("Error running ProcessBuilder to store vcfs to tiledb ${uri} array: ${error}")
            }

        } catch (exc: Exception) {
            myLogger.error("Error: loadVCFToTiledb: could not load tiledb array ${uri}.  See ${redirectError} and ${redirectOutput} for more details.")
            throw IllegalArgumentException("Error running ProcessBuilder to store vcfs to tiledb array: ${exc.message}")
        }
        myLogger.info("Finished! vcfs stored to tiledb dataset ${uri}")
    }

}