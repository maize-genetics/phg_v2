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
class SetupEnvironment : CliktCommand() {
    // The output directory will store the log and error files from the ProcessBUilder conda create command
    val outputDir by option("-o", "--outputDir", help = "Directory where ProcessBuilder, conda create env command log files will be written")
        .required()


        override fun run() {
            // initial test
            println("This is the SetupEnvironment script, which does nothing yet!")
            val envFile = "src/main/resources/environment.yaml"
            // call ProcessBuilder to execute the conda create env command
            val builder = ProcessBuilder("conda", "env", "create", "--file", envFile)
            var redirectOutput = outputDir + "/condaCreate_output.log"
            var redirectError = outputDir + "/condaCreate_error.log"
            builder.redirectOutput( File(redirectOutput))
            builder.redirectError( File(redirectError))
            println(" begin conda create Command:" + builder.command().stream().collect(Collectors.joining(" ")));
            var process = builder.start()
            var error = process.waitFor()

            if (error != 0) {
                println("conda env create command for file ${envFile} run via ProcessBuilder returned error code $error")
                println("Verify you have conda installed and the environment phgv2-conda does not already exist")
                println("If the environment already exists, you can activate it with the command: conda activate phgv2-conda")
                throw IllegalStateException("SetupEnvironment: create conda envirionment failed: $error")
            }
            println("Successfully created your environment.  You can activate it with the command: conda activate phgv2-conda")

        }
}