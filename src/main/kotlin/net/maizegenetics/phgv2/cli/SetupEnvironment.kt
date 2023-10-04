package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand

/**
 * Class to create the user conda environment for PHG.
 * The environment will not be activated until a command is run
 * that requires specific packages.
 * The purpose of this script is to ensure the environment is created
 * and available.
 */
class SetupEnvironment : CliktCommand() {

        override fun run() {
            // initial test
            println("This is the SetupEnvironment script, which does nothing yet!")
        }
}