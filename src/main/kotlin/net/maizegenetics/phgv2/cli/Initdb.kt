package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.arguments.argument
import org.apache.logging.log4j.LogManager

/**
 * This will create a TileDB folder for the given name.
 */
class Initdb : CliktCommand() {

    private val myLogger = LogManager.getLogger(Initdb::class.java)

    private val name by argument(name = "name")

    override fun run() {

    }

}