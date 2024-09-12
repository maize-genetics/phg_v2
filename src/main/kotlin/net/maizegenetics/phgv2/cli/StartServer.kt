package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.types.int
import io.ktor.server.engine.*
import io.ktor.server.netty.*
import net.maizegenetics.phgv2.brapi.service.VariantSetsService
import net.maizegenetics.phgv2.utils.setupDebugLogging
import net.maizegenetics.phgv2.utils.verifyURI
import org.apache.logging.log4j.LogManager
import java.nio.file.Files
import java.nio.file.Paths

/**
 * This class starts the PHGv2 BrAPI server.
 * It is called from the main function in src/main/kotlin/net/maizegenetics/phgv2/cli/PHGv2.kt
 * It takes a single optional parameter, dbPath, which defines the location of the folder containing the TileDB datasets.
 * If dbPath is not supplied, the application.conf file is checked for a TILEDB_URI line.
 * If that is not found, an error is thrown.
 * If dbPath is supplied, we verify that it contains a valid hvcf_dataset.
 *
 */
object StartServer : CliktCommand(help = "Starts PHGv2 BrAPI Server") {

    private val myLogger = LogManager.getLogger(StartServer::class.java)

    val dbPath by option(help = "Full path to folder where TileDB datasets are stored.  \nThis must be run at least once before starting the server. \nIf you have already run it for this server instance you do not need to supply this again.")
        .default("")

    val port by option(help = "The port on which the server will listen.")
        .int()
        .default(8080)

    val condaEnvPrefix by option (help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
        .default("")

    lateinit var server: NettyApplicationEngine

    fun serverURL(): String {
        val actualPort = server.environment.connectors.map { it.port }.firstOrNull()
        val host = server.environment.connectors.map { it.host }.firstOrNull() ?: "localhost"

        return "http://$host:$actualPort/brapi/v2/"
    }

    override fun run() {

        setupDebugLogging()
        // Do not default dbPath to the current working folder.  User may have previously started
        // the server, which will set the TILEDB_URI in the application.conf file.
        // We want to use that if no dbPath is specified.
        // If one is specified, it may be that the user wants to start the server with a different tiledb folder of datasets.
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
            myLogger.info("start-server: using existing TILEDB_URI value from application.conf file: ${configPath}.")
            tiledbPath = configPath // use the configPath already stored in the application file
        } else { // dbPath has a value
            if (configPath != null && configPath.isNotBlank()) {
                if (configPath != dbPath) {
                    if (!verifyURI(dbPath, "hvcf_dataset",condaEnvPrefix)) {
                        myLogger.error("start-server:  \ndp-path does not contain a valid tiledb created hvcf_dataset.  \nPlease re-run start-server with a valid value for db-path parameter.")
                        throw IllegalArgumentException("start-server:  \nTILEDB_URI is not valid.  \nPlease re-run start-server with a valid value for dbPath parameter.")
                    }
                    myLogger.info("start-server:  Updating TILEDB_URI in the application.conf with value from db-path: ${dbPath}")
                } else {
                    myLogger.info("\nstart-server:  Running server with db-path/TILEDB_URI of ${configPath}.")
                }
            } else {
                if (!verifyURI(tiledbPath, "hvcf_dataset",condaEnvPrefix)) {
                    myLogger.error("hvcf_dataset does not exist in $dbPath.  Please check your path, and/or run Initdb to create the datasets.")
                    throw IllegalArgumentException("A valid hvcf_dataset does not exist in $dbPath.")
                }
            }
        }
        // Update the application.conf file with the dbPath value and port number
        myLogger.info("start-server: updating application.conf file with dbPath = $tiledbPath and port = $port")
        updateConfigFile(tiledbPath, port.toString(), appHome)

        VariantSetsService.createAllSamplesHVCF()

        // Checks have passed - Ready to start the server!
        // commandLineEnvironment reads the application.config file
        // https://ktor.io/docs/configuration.html#hocon-file
        server = embeddedServer(Netty, commandLineEnvironment(emptyArray()))
        server.start(wait = true)

    }

    fun updateConfigFile(dbPath: String, port: String, appHome: String) {
        val configPath = Paths.get("${appHome}/resources/main/application.conf")
        // write the existing file minus the TILEDB_URI line
        var config = ""
        Files.lines(Paths.get(configPath.toString())).forEach { line ->
            if (!line.startsWith("TILEDB_URI") && !line.startsWith("PORT")) {
                config += line + "\n"
            }
        }

        // append to the beginning of config the line "TILEDB_URI=${dbPath}"
        config = "TILEDB_URI=${dbPath}\nPORT=${port}\n" + config
        // write the new config file - do I have permission ???
        configPath.toFile().writeText(config)
    }

    // Find the path to this class from the jar file
    // Use that to determine the path to the application.conf file
    fun getClassPath(): String {
        // a little bit of help from ecerer answer in https://stackoverflow.com/questions/227486/find-where-java-class-is-loaded-from
        val c: Class<*> = StartServer::class.java
        val path = c.getResource(c.getSimpleName() + ".class").path.replace(c.getSimpleName() + ".class", "")
        val appHome = path.substringAfter("file:").substringBefore("lib/phg_v2.jar")

        myLogger.info("getClassPath: \n  path = $path \n  appHome = $appHome")
        return appHome
    }

    // Find the TILEDB_URI from the application.conf file
    fun getDbPathFromConfigFile(appHome: String): String? {
        val configPath = Paths.get("${appHome}/resources/main/application.conf")

        var config = Files.readString(Paths.get(configPath.toString()))
        // Find the line in the config file that starts with "TILEDB_URI"
        // and extract the path
        // This will return null if the line is not found
        return config.lines().find { it.startsWith("TILEDB_URI") }?.substringAfter("=")?.substringBefore("\n")
    }
}