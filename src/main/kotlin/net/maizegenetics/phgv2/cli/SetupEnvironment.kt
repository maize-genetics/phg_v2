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
        PHG, TILEDB
    }

    companion object {
        internal fun condaEnvironmentName(envType: CONDAENVTYPE): String = when (envType) {
            CONDAENVTYPE.PHG -> "phgv2-conda"
            CONDAENVTYPE.TILEDB -> "phgv2-tiledb"
        }

        internal fun isExistingEnvironmentFailure(
            errorLines: List<String>,
            outputLines: List<String>
        ): Boolean {
            val combinedLines = errorLines + outputLines
            return combinedLines.any {
                it.contains("prefix already exists", ignoreCase = true) ||
                    it.contains("environment already exists", ignoreCase = true)
            }
        }

        internal fun summarizeCondaFailure(
            errorLines: List<String>,
            outputLines: List<String>,
            errorLogPath: String
        ): String {
            return errorLines.firstOrNull { it.isNotBlank() }
                ?: outputLines.firstOrNull { it.isNotBlank() }
                ?: "No conda output captured. See $errorLogPath for details."
        }

        internal fun isTermsOfServiceFailure(
            errorLines: List<String>,
            outputLines: List<String>
        ): Boolean {
            val combinedLines = errorLines + outputLines
            return combinedLines.any {
                it.contains("CondaToSNonInteractiveError", ignoreCase = true) ||
                    it.contains("Terms of Service have not been accepted", ignoreCase = true)
            }
        }

        internal fun termsOfServiceHelpMessage(): String {
            return buildString {
                appendLine("Conda is blocking access because channel Terms of Service have not been accepted.")
                appendLine("Accept the ToS for the required channels, then rerun setup-environment.")
                appendLine("Typical commands are:")
                appendLine("  conda tos accept")
                appendLine("  conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main")
                appendLine("  conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r")
                append("For CI or other non-interactive environments, you can also set CONDA_PLUGINS_AUTO_ACCEPT_TOS=yes.")
            }
        }
    }

    val envFile by option("-e", "--env-file", help = "File containing the conda environment definition")
        .default("")

    val tiledbEnvFile by option(help = "File containing the conda environment definition for tiledb")
        .default("")

    // This function uses ProcessBuilder to setup the phgv2 conda environment
    fun createEnvironment(envFile: String, outputDir: String, tiledbFile: String) {

        // if no user defined environment file, use the default
        val resultEnvFile = if (envFile == "") {

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

        val tiledbEnvFile = if (tiledbFile == "") {
            val selectedEnvFile = "${outputDir}/phg_tiledb_environment.yml"
            val fileContent = this::class.java.classLoader.getResource("phg_tiledb_environment.yml").readText()
            myLogger.info("writing default environment file to $selectedEnvFile")
            // This will overwrite an existing file
            File(selectedEnvFile).writeText(fileContent)
            selectedEnvFile
        } else {
            if (!File(tiledbFile).exists()) {
                myLogger.error("The specified file $tiledbFile does not exist.")
                throw IllegalStateException("SetupEnvironment: create conda environment failed: file $tiledbFile does not exist")
            }
            tiledbFile
        }

        runCondaCreate(resultEnvFile, outputDir, CONDAENVTYPE.PHG)
        runCondaCreate(tiledbEnvFile, outputDir, CONDAENVTYPE.TILEDB)
    }

    private fun runCondaCreate(resultEnvFile: String, outputDir: String, envType: CONDAENVTYPE) {
        // call ProcessBuilder to execute the conda create env command
        myLogger.info("Creating conda environment from file: $resultEnvFile")
        val builder = ProcessBuilder("conda", "env", "create", "--solver=libmamba", "--file", resultEnvFile)
        val redirectOutput = when(envType) {
            CONDAENVTYPE.PHG -> "$outputDir/condaCreate_output.log"
            CONDAENVTYPE.TILEDB -> "$outputDir/tiledbCondaCreate_output.log"
        }

        val redirectError = when(envType) {
            CONDAENVTYPE.PHG -> "$outputDir/condaCreate_error.log"
            CONDAENVTYPE.TILEDB -> "$outputDir/tiledbCondaCreate_error.log"
        }

        builder.redirectOutput(File(redirectOutput))
        builder.redirectError(File(redirectError))
        myLogger.info(" begin conda create Command:" + builder.command().joinToString(" "));
        val process = builder.start()
        val error = process.waitFor()

        if (error != 0) {
            val errorLines = File(redirectError).takeIf { it.exists() }?.readLines() ?: emptyList()
            val outputLines = File(redirectOutput).takeIf { it.exists() }?.readLines() ?: emptyList()
            val envName = condaEnvironmentName(envType)

            if (!isExistingEnvironmentFailure(errorLines, outputLines)) {
                val failureSummary = summarizeCondaFailure(errorLines, outputLines, redirectError)
                val tosHelp = if (isTermsOfServiceFailure(errorLines, outputLines)) {
                    "\n${termsOfServiceHelpMessage()}"
                } else {
                    ""
                }
                myLogger.error(
                    "conda env create command for file $resultEnvFile returned error code $error. " +
                        "See $redirectError for details. First message: $failureSummary$tosHelp"
                )
                println("Conda environment creation failed for $resultEnvFile.")
                println("See $redirectError for details.")
                if (isTermsOfServiceFailure(errorLines, outputLines)) {
                    println(termsOfServiceHelpMessage())
                }
                throw IllegalStateException(
                    "SetupEnvironment: create conda environment failed for $resultEnvFile " +
                        "with exit code $error. $failureSummary$tosHelp"
                )
            } else {
                myLogger.info(
                    "\nThe $envName conda environment already exists.  " +
                        "You can activate it with the command: conda activate $envName\n"
                )
            }
        } else {
            val envName = condaEnvironmentName(envType)
            myLogger.info(
                "Successfully created your environment.  " +
                    "You can activate it with the command: conda activate $envName"
            )
        }
    }

    override fun run() {
        logCommand(this)

        val workingDir = System.getProperty("user.dir")

        // call method to create the environment
        createEnvironment(envFile, workingDir, tiledbEnvFile)
    }
}