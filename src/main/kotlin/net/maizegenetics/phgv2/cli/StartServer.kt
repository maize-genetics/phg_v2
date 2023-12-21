package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import io.ktor.server.engine.*
import io.ktor.server.netty.*
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.apache.logging.log4j.LogManager

class StartServer : CliktCommand(help = "Starts PHGv2 BrAPI Server") {

    private val myLogger = LogManager.getLogger(StartServer::class.java)

    override fun run() {
        setupDebugLogging()

        // commandLineEnvironment reads the application.config file
        // https://ktor.io/docs/configuration.html#hocon-file
        embeddedServer(Netty, commandLineEnvironment(emptyArray())).start(wait = true)
    }

}