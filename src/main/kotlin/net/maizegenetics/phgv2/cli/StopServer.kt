package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import io.ktor.server.engine.*
import io.ktor.server.netty.*
import org.apache.logging.log4j.LogManager

/**
 * This will stop the PHGv2 BrAPI Server
 * It should stop the same server setup by STartServer as it uses the same
 * commandLineEnvironment.  We do not need to pass any args to the server
 * as it will not make use of the dbUri.
 *
 * https://dev.to/viniciusccarvalho/graceful-shutdown-of-ktor-applications-1h53
 */
class StopServer: CliktCommand(help = "Stops PHGv2 BrAPI Server") {
    private val myLogger = LogManager.getLogger(StopServer::class.java)
    override fun run() {
        myLogger.info("Stopping PHGv2 BrAPI Server")
        embeddedServer(Netty, commandLineEnvironment(emptyArray())).stop()
    }
}