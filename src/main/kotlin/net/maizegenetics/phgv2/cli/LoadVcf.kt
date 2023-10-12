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
class LoadVcf : CliktCommand() {

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

    val tempDir by option(help = "Folder where temporary files will be written")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--temp-dir must not be blank"
            }
        }

    override fun run() {

        // Check the type of files in the vcfDir
        // antying with h.vcf.gz or hvcf.gz is a hvcf file
        // anything with g.vcf.gz or gvcf.gz is a gvcf file
        val fileLists = getFileLists(vcfDir)

        if (fileLists.first.isEmpty() && fileLists.second.isEmpty()) {
            // no gvcf or hvcf files found in the user supplied vcfDir
            myLogger.warn("No gvcf or hvcf files found in $vcfDir.  Please check the folder and try again.")
            return
        }

        // verify the URI is correct, create the datasets if they don't exist
        // verifyURI() will create both datasets if they don't exist
        if (fileLists.first.isNotEmpty()) {

        }
       // val goodDbPath = verifyURI(dbPath) - add uri, determined by the list you have
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
                gvcfFiles.add(file.name)
            } else if (file.name.endsWith("h.vcf.gz") || file.name.endsWith("hvcf.gz")) {
                hvcfFiles.add(file.name)
            }
        }
        return Pair(gvcfFiles,hvcfFiles)
    }

    // The uri should either be gvcf_dataset or hvcf_dataset
    // The user determines the parent folder name where these datasets live
    // The actual tiledb dataset names are constant and are either gvcf_dataset or hvcf_dataset

    // Are we just verifying the uri here, or creating it as well?  I think it would be ok to
    // create them here.
    fun verifyURI(dbPath:String,uri:String): Boolean {
        // Check that the user supplied db folder exists
        check(File(dbPath).exists()) { "Folder $dbPath does not exist - please send a valid path that indicates the parent folder for your tiledb datasets." }

        // Check if the dataset exists
        val dbfolder = File(dbPath)
        val dataset = "${dbfolder.name}/${uri}"

        val datasetPath = Paths.get(dataset)
        if (Files.isRegularFile(datasetPath)) {
            throw Exception("URI ${dataset}is a file, not a tiledb dataset folder.  The parent folder must not contain any files/folders named gvcf_dataset or hvcf_dataset that is not a tiledb created URI")
        }

        if (File(dataset).exists()  && Files.isDirectory(Paths.get(dataset))){
            var builder = ProcessBuilder("conda","run","-n","phgv2-conda","tiledbvcf","stat","--uri",uri)
            var redirectOutput = tempDir + "/tiledb_statURI_output.log"
            var redirectError = tempDir + "/tiledb_statURI_error.log"
            builder.redirectOutput( File(redirectOutput))
            builder.redirectError( File(redirectError))

            // verify if the output.log contains "Version"
            // if not, then the URI is not a tiledbvcf URI
            println("begin Command:" + builder.command().stream().collect(Collectors.joining(" ")))

            try {
                var process = builder.start()
                var error = process.waitFor()
                if (error != 0) {
                    println("LoadTiledbH tiledb stat returned error code $error")
                    throw IllegalArgumentException("Error: URI is not a tiledb URI folder created via the tiledb create command: ${error}")
                }

            } catch (exc: Exception) {
                println("Error: could not run tiledb stat on ${uri}.")
                throw IllegalArgumentException("Error running ProcessBuilder to stat tiledb URI: ${exc.message}")
            }
            // check if is a tiledb dataset

            myLogger.info("Using existing TileDB datasets previously created in folder $dbPath.")
            return true
        } else {
            myLogger.info("TileDB datasets not found in folder $dbPath. Please run InitDB to create the datasets.")
            return false
        }

    }

}