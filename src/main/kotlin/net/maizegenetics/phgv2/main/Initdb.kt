package net.maizegenetics.phgv2.main

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.arguments.argument
import org.apache.logging.log4j.LogManager
import java.nio.file.Files
import java.nio.file.Paths
import java.util.stream.Collectors

/**
 * This will create a TileDB folder for the given name.
 */
class Initdb : CliktCommand() {

    private val myLogger = LogManager.getLogger(Initdb::class.java)

    private val name by argument(name = "name")

    override fun run() {

    }

}