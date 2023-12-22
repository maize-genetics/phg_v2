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
import java.nio.file.Paths
import java.util.stream.Collectors

class StartServer : CliktCommand(help = "Starts PHGv2 BrAPI Server") {

    private val myLogger = LogManager.getLogger(StartServer::class.java)

    val dbPath by option(help = "Full path to folder where TileDB datasets are stored.   ")
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
        if (!hvcfExists) {
            myLogger.error("hvcf_dataset does not exist in $dbPath.  Exiting.")
            return
        }

        // Create and Args list to pass to the server
        // This tells the endpoint code where the datasets are located.
        val args = listOf("TILEDB_URI=${dbPath}").toTypedArray()

        // commandLineEnvironment reads the application.config file
        // https://ktor.io/docs/configuration.html#hocon-file
        embeddedServer(Netty, commandLineEnvironment(args)).start(wait = true)
    }



}