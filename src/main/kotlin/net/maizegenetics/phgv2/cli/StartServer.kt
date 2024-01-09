package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import io.ktor.server.engine.*
import io.ktor.server.netty.*
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.apache.logging.log4j.LogManager
import java.nio.file.Files
import java.nio.file.Paths

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

        val appHome = getClassPath()

        // write dbpath to application.conf
        writeConfigFile(dbPath,appHome)

        println("start-server:  config file successfully written with dbPath = ${dbPath}")
        // commandLineEnvironment reads the application.config file
        // https://ktor.io/docs/configuration.html#hocon-file
        embeddedServer(Netty, commandLineEnvironment(emptyArray())).start(wait = true)
    }

    // write the dbPath to the application.conf file using
    // the appHome variable created above
    fun writeConfigFile(dbPath:String,appHome:String) {

        // the application file is in ${appHome}/resources
        // NOTE for the CLASSPATH, APP_HOME does not include the "phg" resources,
        // but it does when calculated here.

        // Create a configPath using the appHome variable created above
        val configPath = Paths.get("${appHome}/resources/application.conf")


        myLogger.info("writeConfigFile: configPath = ${configPath}")
        var config = Files.readString(Paths.get(configPath.toString()))

        // append to the beginning of config the line "TILEDB_URI=${dbPath}"
        config  = "TILEDB_URI=${dbPath}\n" + config
        //write the new config file - do I have permission ???
        configPath.toFile().writeText(config)
    }


    // Find the path to this class from the jar file
    // Use that to determine the path to the application.conf file
    fun getClassPath():String {
        // little bit of help from ecerer answer in https://stackoverflow.com/questions/227486/find-where-java-class-is-loaded-from
        val c: Class<*> = StartServer::class.java
        val path = c.getResource(c.getSimpleName() + ".class").path.replace(c.getSimpleName() + ".class", "")
        val appHome = path.substringAfter("file:").substringBefore("lib/phg_v2.jar")

        myLogger.info("getClassPath: \n  path = $path \n  appHome = $appHome")
        return appHome
    }
}