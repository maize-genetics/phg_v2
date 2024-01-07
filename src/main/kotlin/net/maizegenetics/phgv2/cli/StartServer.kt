package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import io.ktor.serialization.kotlinx.json.*
import io.ktor.server.application.*
import io.ktor.server.engine.*
import io.ktor.server.netty.*
import io.ktor.server.plugins.callloging.*
import io.ktor.server.plugins.contentnegotiation.*
import io.ktor.server.plugins.defaultheaders.*
import io.ktor.server.routing.*
import kotlinx.serialization.json.Json
import net.maizegenetics.phgv2.brapi.api.apiRoute
import net.maizegenetics.phgv2.brapi.module

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

        // Create an Args list to pass to the server
        // This tells the endpoint code where the datasets are located.
        val dbUri = "-P:TILEDB_URI=${dbPath}"
        //val dbUri = "TILEDB_URI=${dbPath}"
        val args = arrayOf(dbUri)

        // commandLineEnvironment reads the application.config file
        // https://ktor.io/docs/configuration.html#hocon-file
        //embeddedServer(Netty, commandLineEnvironment()).start(wait = true)

//        val server = embeddedServer(Netty, commandLineEnvironment(args)) {
//            Application.module(args)
//        }
        val server = embeddedServer(Netty,port=8080) {
            module(args)
        }
        server.start(true)

    }

}