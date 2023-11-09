package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import java.io.File
import java.util.logging.Logger


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

    private val myLogger = Logger.getLogger("net.maizegenetics.phgv2.cli.SetupEnvironment")

    val envFile by option("-e", "--envFile", help = "File containing the conda environment definition")
        .default("")

    // This function uses ProcssBuilder to setup the phgv2 conda environment
    fun createEnvironment( envFile: String, outputDir:String) {

        var selectedEnvFile = envFile
        myLogger.info("begin run: selectedEnvFile = $selectedEnvFile")

        // if no user defined environment file, use the default
        if (selectedEnvFile == "") {

            selectedEnvFile = "${outputDir}/environment.yml"
            val inputUrl = SetupEnvironment::class.java.getResource("environment.yml")
            println("inputUrl = ${inputUrl}, selectedEnvFile = $selectedEnvFile")

            val fileContent = this::class.java.classLoader.getResource("environment.yml").readText()
            // This will overwrite an existing file
            File(selectedEnvFile).writeText(fileContent)
        }

        // call ProcessBuilder to execute the conda create env command
        myLogger.info("Creating conda environment from file: $envFile")
        val builder = ProcessBuilder("conda", "env", "create","--solver=libmamba", "--file", selectedEnvFile)
        var redirectOutput = outputDir + "/condaCreate_output.log"
        var redirectError = outputDir + "/condaCreate_error.log"
        builder.redirectOutput( File(redirectOutput))
        builder.redirectError( File(redirectError))
        myLogger.info(" begin conda create Command:" + builder.command().joinToString(" "));
        var process = builder.start()
        var error = process.waitFor()

        if (error != 0) {
            val errorLines = File(redirectError).readLines()
            var success = errorLines
                .filter { it.contains("prefix already exists") }
                .toList()
            if (success == null || success.size == 0) {
                myLogger.severe("conda env create command for file ${envFile} run via ProcessBuilder returned error code $error")
                println("Verify you have conda installed and the environment phgv2-conda does not already exist")
                throw IllegalStateException("SetupEnvironment: create conda envirionment failed: $error")
            } else {
                myLogger.info("\nThe phgv2-conda conda environment already exists.  You can activate it with the command: conda activate phgv2-conda\n")
            }
        } else {
            myLogger.info("Successfully created your environment.  You can activate it with the command: conda activate phgv2-conda")
        }
    }

        override fun run() {
            val workingDir = System.getProperty("user.dir")

            // call method to create the environment
            createEnvironment(envFile,workingDir)
        }
}