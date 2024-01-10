package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import io.ktor.server.engine.*
import io.ktor.server.netty.*
import net.maizegenetics.phgv2.utils.setupDebugLogging
import net.maizegenetics.phgv2.utils.verifyURI
import org.apache.logging.log4j.LogManager
import java.nio.file.Files
import java.nio.file.Paths

class StartServer : CliktCommand(help = "Starts PHGv2 BrAPI Server") {

    private val myLogger = LogManager.getLogger(StartServer::class.java)
    val dbPath by option(help = "Full path to folder where TileDB datasets are stored.  \nThis must be run at least once before starting the server. \nIf you have already run it for this server instance you need not supply this again.")
        .default("")


    override fun run() {
        setupDebugLogging()
        var tiledbPath = dbPath

        val appHome = getClassPath()
        // Check if the application.conf file has a TILEDB_URI line that is not commented out.
        val configPath = getDbPathFromConfigFile(appHome)

        // dbPath must either be passed, or it must already exist in the application.conf file.
        // If application.conf file has a value, use that.
        // Otherwise, use the dbPath value.
        // If dbPath is blank, and there is no value in the application.conf file, error.
        if (dbPath.isBlank()) {
            myLogger.info("start-server:  dbPath is blank.  Reading application.conf file to get TILEDB_URI")
            if (configPath == null) {
                myLogger.error("start-server:  \nTILEDB_URI is not set in application.conf file.  \nPlease re-run start-server with a valid value for dbPath parameter.")
                throw IllegalArgumentException("start-server:  \nTILEDB_URI is not set in application.conf file.  \nPlease re-run start-server with a valid value for dbPath parameter.")
            }
        } else { // dbPath has a value
            if (configPath!=null && configPath.isNotBlank()) {
                myLogger.info("\nstart-server:  TILEDB_URI is already set in application.conf file.  \nUsing the application.conf TILEDB_URI value of ${configPath}.")
            } else {
                // verify the URI is good
                val hvcfExists = verifyURI(tiledbPath,"hvcf_dataset")
                if (!hvcfExists) {
                    myLogger.error("hvcf_dataset does not exist in $dbPath.  Please check your path, and/or run Initdb to create the datasets.")
                    throw IllegalArgumentException("A valid hvcf_dataset does not exist in $dbPath.")
                }
                myLogger.info("start-server:  dbPath = ${dbPath}, writing to application.conf file")
                writeConfigFile(dbPath,appHome)
            }
        }

        // Checks have passed - Ready to start the server!
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

    // Find the TILEDB_URI from the application.conf file
    fun getDbPathFromConfigFile(appHome:String):String? {
        val configPath = Paths.get("${appHome}/resources/application.conf")

        var config = Files.readString(Paths.get(configPath.toString()))
        // Find the line in the config file that starts with "TILEDB_URI"
        // and extract the path
        // This will return null if the line is not found
        val tiledbPath = config.lines().find { it.startsWith("TILEDB_URI") }?.substringAfter("=")?.substringBefore("\n")
        return tiledbPath
    }
}