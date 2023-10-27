package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import java.io.File
import java.util.stream.Collectors

/**
 * Class to create the user conda environment for PHG.
 * The environment will not be activated until a command is run
 * that requires specific packages.
 * The purpose of this script is to ensure the environment is created
 * and available.
 *
 * Check if this environment already exists.  If an error occurs when
 * creating it, do we assume it exists?
 */
class SetupEnvironment : CliktCommand(help="create conda environment for PHG") {
    // The output directory will store the log and error files from the ProcessBUilder conda create command
    val outputDir by option("-o", "--outputDir", help = "Directory where ProcessBuilder, conda create env command log files will be written")
        .required()

    // This function uses ProcssBuilder to setup the phgv2 conda environment
    fun createEnvironment( envFile: String, outputDir:String) {
        // call ProcessBuilder to execute the conda create env command
        val builder = ProcessBuilder("conda", "env", "create","--solver=libmamba", "--file", envFile)
        var redirectOutput = outputDir + "/condaCreate_output.log"
        var redirectError = outputDir + "/condaCreate_error.log"
        builder.redirectOutput( File(redirectOutput))
        builder.redirectError( File(redirectError))
        println(" begin conda create Command:" + builder.command().joinToString(" "));
        var process = builder.start()
        var error = process.waitFor()

        if (error != 0) {
            val errorLines = File(redirectError).readLines()
            var success = errorLines
                .filter { it.contains("prefix already exists") }
                .toList()
            if (success == null || success.size == 0) {
                println("conda env create command for file ${envFile} run via ProcessBuilder returned error code $error")
                println("Verify you have conda installed and the environment phgv2-conda does not already exist")
                throw IllegalStateException("SetupEnvironment: create conda envirionment failed: $error")
            } else {
                println("\nThe phgv2-conda conda environment already exists.  You can activate it with the command: conda activate phgv2-conda\n")
            }
        } else {
            println("Successfully created your environment.  You can activate it with the command: conda activate phgv2-conda")
        }
    }

        override fun run() {
            // This is the official phgv2 conda environment file
            val envFile = "src/main/resources/environment.yml"

            // call method to create the environment
            createEnvironment(envFile,outputDir)
        }
}