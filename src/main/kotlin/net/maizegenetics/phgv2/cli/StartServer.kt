package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import io.ktor.server.engine.*
import io.ktor.server.netty.*
import net.maizegenetics.phgv2.utils.setupDebugLogging
import net.maizegenetics.phgv2.utils.verifyURI
import org.apache.logging.log4j.LogManager
import java.io.File
import java.nio.file.Files
import java.nio.file.Paths

class StartServer : CliktCommand(help = "Starts PHGv2 BrAPI Server") {

    private val myLogger = LogManager.getLogger(StartServer::class.java)

    val dbPath by option(help = "Folder name where TileDB datasets are stored.   ")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--db-path must not be blank"
            }
        }

    override fun run() {
        setupDebugLogging()

        // Verify the uri is valid.  We only care about the hvcf dataset,
        // so check that one explicitly
        val hvcfExists = verifyURI(dbPath,"hvcf_dataset")

        // commandLineEnvironment reads the application.config file
        // https://ktor.io/docs/configuration.html#hocon-file
        embeddedServer(Netty, commandLineEnvironment(emptyArray())).start(wait = true)
    }



}