package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.arguments.argument
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
class Initdb : CliktCommand() {

    private val myLogger = LogManager.getLogger(Initdb::class.java)

    val dbpath by option(help = "Folder name where TileDB datasets will be created")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--dbpath must not be blank"
            }
        }

    fun createDataSets(dbpath:String) {
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


        // create a temp folder to hold the output and error logs
        val tempDir = dbpath + "/temp"
        Files.createDirectories(Paths.get(tempDir))

        // create the hvcf dataset.  The processBuilder command must first set the conda
        // environment to phgv2-conda.  It is presumed that the user has already created
        // this environment.

        var builder = ProcessBuilder("conda","run","-n","phgv2-conda","tiledbvcf","create","--uri",hvcf_dataset)
        var redirectOutput = tempDir + "/tiledb_hvcf_createURI_output.log"
        var redirectError = tempDir + "/tiledb_hvcf_createURI_error.log"
        builder.redirectOutput( File(redirectOutput))
        builder.redirectError( File(redirectError))

        println("begin Command to create hvcf dataset:" + builder.command().stream().collect(Collectors.joining(" ")))
        myLogger.info("begin Command to create gvcf dataset:" + builder.command().joinToString(" "))
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
        builder = ProcessBuilder("conda","run","-n","phgv2-conda","tiledbvcf","create","--uri",gvcf_dataset)
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

        // call method to create the environment
        createDataSets(dbpath)
    }

}