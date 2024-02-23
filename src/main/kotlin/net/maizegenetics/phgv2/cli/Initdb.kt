package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.int
import org.apache.logging.log4j.LogManager
import java.io.File
import java.lang.Exception
import java.nio.file.Files
import java.nio.file.Paths
import java.util.stream.Collectors

/**
 * This will create Tiledb datasets for the g.vcf and h.vcf files.
 * Input is a folder name where the datasets will be created.  If this
 * folder does not exist, it will be created.
 *
 * Result:  2 datasets created, which appear as sub-directories of the input folder.
 *   gvcf_dataset
 *   hvcf_dataset
 *
 * The code will check for the existence of the datasets and will not overwrite them.
 */
class Initdb : CliktCommand(help = "Create TileDB datasets for g.vcf and h.vcf files") {

    private val myLogger = LogManager.getLogger(Initdb::class.java)

    val dbPath by option(help = "Folder name under which TileDB datasets will be created. If this folder does not exist, it will be created.")
        .default("")

    val gvcfAnchorGap by option("-g", "--gvcf-anchor-gap", help = "tiledbvcf --anchor-gap for gvcf database (default: 1000000)")
        .int()
        .default(1000000)

    val hvcfAnchorGap by option("-h", "--hvcf-anchor-gap", help = "tiledbvcf --anchor-gap for hvcf database (default: 1000)")
        .int()
        .default(1000)

    fun createDataSets(dbpath:String, gvcfAnchorGap: Int = 1000000, hvcfAnchorGap: Int = 1000) {
        // Check that the user supplied folder exists
        val dbfolder = File(dbpath)

        if (!dbfolder.exists()) {
            myLogger.warn("Folder $dbpath does not exist - creating.")
            dbfolder.mkdirs()
        }

        // check if tiledb datasets already exist.  If so, do not overwrite them.
        val gvcf_dataset = dbpath + "/gvcf_dataset"
        val hvcf_dataset = dbpath + "/hvcf_dataset"

        if (File(gvcf_dataset).exists() || File(hvcf_dataset).exists()) {
            myLogger.warn("TileDB datasets already exist in folder $dbpath.  \nIf $gvcf_dataset or $hvcf_dataset are not tiledb datasets, then delete and run again or chose a different base folder to house your tiledb data.")
            return
        }

        val tempDir = "${dbpath}/temp"
        Files.createDirectories(Paths.get(tempDir))

        // create the hvcf dataset.  The processBuilder command must first set the conda
        // environment to phgv2-conda.  It is presumed that the user has already created
        // this environment.
        var logFile = "${tempDir}/tiledbvcf_createHvcf.log"

        var builder = ProcessBuilder("conda","run","-n","phgv2-conda","tiledbvcf","create","--uri",hvcf_dataset,"-n","--log-level","debug","--log-file",logFile,"--anchor-gap",hvcfAnchorGap.toString())
        var redirectOutput = tempDir + "/tiledb_hvcf_createURI_output.log"
        var redirectError = tempDir + "/tiledb_hvcf_createURI_error.log"
        builder.redirectOutput( File(redirectOutput))
        builder.redirectError( File(redirectError))

        println("begin Command to create hvcf dataset:" + builder.command().stream().collect(Collectors.joining(" ")))
        myLogger.info("begin Command to create hvcf dataset:" + builder.command().joinToString(" "))
        try {
            var process = builder.start()
            var error = process.waitFor()
            if (error != 0) {
                myLogger.error("Initdb run via ProcessBuilder returned error code $error")
                throw IllegalArgumentException("Error running ProcessBuilder to create tiledb hvcf array: ${error}")
            }

        } catch (exc: Exception) {
            myLogger.info("Error: could not create tiledb array $hvcf_dataset.")
            throw IllegalArgumentException("Error running ProcessBuilder to create tiledb array: ${exc.message}")
        }

        // Now create the gvcf dataset
        logFile = "${tempDir}/tiledbvcf_createHvcf.log"
        builder = ProcessBuilder("conda","run","-n","phgv2-conda","tiledbvcf","create","--uri",gvcf_dataset,"-n","--log-level","debug","--log-file",logFile,"--anchor-gap",gvcfAnchorGap.toString())
        redirectOutput = tempDir + "/tiledb_gvcf_createURI_output.log"
        redirectError = tempDir + "/tiledb_gvcf_createURI_error.log"
        builder.redirectOutput( File(redirectOutput))
        builder.redirectError( File(redirectError))

        println("begin Command to create gvcf dataset:" + builder.command().stream().collect(Collectors.joining(" ")))
        myLogger.info("begin Command to create gvcf:" + builder.command().joinToString(" "))
        try {
            var process = builder.start()
            var error = process.waitFor()
            if (error != 0) {
                myLogger.error("Initdb run via ProcessBuilder returned error code $error")
                throw IllegalArgumentException("Error running ProcessBuilder to create tiledb gvcf array: ${error}")
            }

        } catch (exc: Exception) {
            myLogger.info("Error: could not create tiledb array $gvcf_dataset.")
            throw IllegalArgumentException("Error running ProcessBuilder to create tiledb gvcf array: ${exc.message}")
        }
    }
    override fun run() {

        // Set the dbPath to the current working directory if it is not provided
        // no validation is done here as the createDataSets method will check for the existence of the folder
        // and create it if doesn't exist.
        val dbPath = if (dbPath.isBlank()) {
            System.getProperty("user.dir")
        } else {
            dbPath
        }

        // call method to create the environment
        createDataSets(dbPath, gvcfAnchorGap, hvcfAnchorGap)
    }

}