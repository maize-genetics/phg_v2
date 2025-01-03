package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import org.apache.logging.log4j.LogManager
import java.io.File


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
class SetupEnvironment : CliktCommand(help = "Create a conda environment for PHG") {

    private val myLogger = LogManager.getLogger(SetupEnvironment::class.java)

    enum class CONDAENVTYPE {
        PHG, ROPEBWT
    }

    val envFile by option("-e", "--env-file", help = "File containing the conda environment definition")
        .default("")

    val ropebwtEnvFile by option(help = "File containing the conda environment definition for ropebwt")
        .default("")

    // This function uses ProcessBuilder to setup the phgv2 conda environment
    fun createEnvironment(envFile: String, outputDir: String, ropeBWTFile:String) {

        // if no user defined environment file, use the default
        var resultEnvFile = if (envFile == "") {

            val selectedEnvFile = "${outputDir}/phg_environment.yml"

            val fileContent = this::class.java.classLoader.getResource("phg_environment.yml").readText()
            myLogger.info("writing default environment file to $selectedEnvFile")
            // This will overwrite an existing file
            File(selectedEnvFile).writeText(fileContent)

            selectedEnvFile

        } else {
            if (!File(envFile).exists()) {
                myLogger.error("The specified file $envFile does not exist.")
                throw IllegalStateException("SetupEnvironment: create conda environment failed: file $envFile does not exist")
            }
            envFile
        }

        var ropeBwtResutlEnvFile = if(ropeBWTFile == "") {
            val selectedEnvFile = "${outputDir}/phg_ropebwt_environment.yml"
            val fileContent = this::class.java.classLoader.getResource("phg_ropebwt_environment.yml").readText()
            myLogger.info("writing default environment file to $selectedEnvFile")
            // This will overwrite an existing file
            File(selectedEnvFile).writeText(fileContent)
            selectedEnvFile
        } else {
            if (!File(ropeBWTFile).exists()) {
                myLogger.error("The specified file $ropeBWTFile does not exist.")
                throw IllegalStateException("SetupEnvironment: create conda environment failed: file $ropeBWTFile does not exist")
            }
            ropeBWTFile

        }

        runCondaCreate(resultEnvFile, outputDir, envFile, CONDAENVTYPE.PHG)
        runCondaCreate(ropeBwtResutlEnvFile, outputDir, ropeBWTFile, CONDAENVTYPE.ROPEBWT)
    }

    private fun runCondaCreate(resultEnvFile: String, outputDir: String, envFile: String, envType: CONDAENVTYPE) {
        // call ProcessBuilder to execute the conda create env command
        myLogger.info("Creating conda environment from file: $resultEnvFile")
        val builder = ProcessBuilder("conda", "env", "create", "--solver=libmamba", "--file", resultEnvFile)
        var redirectOutput = when(envType) {
            CONDAENVTYPE.PHG -> outputDir + "/condaCreate_output.log"
            CONDAENVTYPE.ROPEBWT -> outputDir + "/ropebwtCondaCreate_output.log"
        }

        var redirectError = when(envType) {
            CONDAENVTYPE.PHG -> outputDir + "/condaCreate_error.log"
            CONDAENVTYPE.ROPEBWT -> outputDir + "/ropebwtCondaCreate_error.log"
        }

        builder.redirectOutput(File(redirectOutput))
        builder.redirectError(File(redirectError))
        myLogger.info(" begin conda create Command:" + builder.command().joinToString(" "));
        var process = builder.start()
        var error = process.waitFor()

        if (error != 0) {
            val errorLines = File(redirectError).readLines()
            var success = errorLines
                .filter { it.contains("prefix already exists") }
                .toList()
            if (success == null || success.size == 0) {
                myLogger.error("conda env create command for file ${envFile} run via ProcessBuilder returned error code $error")
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
        logCommand(this)

        val workingDir = System.getProperty("user.dir")

        // call method to create the environment
        createEnvironment(envFile, workingDir, ropebwtEnvFile)
    }
}