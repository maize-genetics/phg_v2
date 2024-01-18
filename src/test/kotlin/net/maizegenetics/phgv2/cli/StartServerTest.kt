package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.brapi.createSmallSeqTiledb
import net.maizegenetics.phgv2.brapi.resetDirs
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.AfterEach
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.BeforeEach
import java.io.File
import java.nio.file.Files
import java.nio.file.Paths
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * It is difficult to test the StartServer command. When we start the server, it doesn't return
 * until the server is stopped.
 * Testing of the dbPath variable depends on checking values in the application.conf file.
 * We identify the location of the application.conf file based on the location of the StartServer class.
 * For non-junits, this will be in <>/phg/resources/application.conf.  But for junit it is in the
 * repository at <>phg/cli/resources/application.conf, so this code fails.
 *
 * That leaves us with the option of testing a couple of the functions in the StartServer class,
 * so that is what is included here.
 *
 */
class StartServerTest {
    companion object {

        @JvmStatic
        @BeforeAll
        fun copyApplicationConfFile() {
            // copy the application.conf file to a saved location
            // in the same folder.
            val appHome = StartServer().getClassPath().substringBefore("classes")
            val configPath = Paths.get("${appHome}/resources/main/application.conf")

            var config = Files.readString(Paths.get(configPath.toString()))

            val configSaved = Paths.get("${appHome}/resources/main/application.conf.saved")
            //write the new config file - do I have permission ???
            configSaved.toFile().writeText(config.toString())

            // Now, for the real file, we need to re-write it with the TILEDB_URI line removed.
            // When CI tests are run, a script writes the TILEDB_URI line to the application.conf file.
            // Copy the original config file with the TILEDB_URI line removed back to the original location.
            config = ""
            Files.lines(Paths.get(configPath.toString())).forEach { line ->
                if (!line.startsWith("TILEDB_URI") && !line.startsWith("PORT")) {
                    config += line + "\n"
                }
            }
            //write the new config file
            configPath.toFile().writeText(config.toString())
        }

        @JvmStatic
        @AfterAll
        fun returnApplicationConfFile() {
            val appHome = StartServer().getClassPath().substringBefore("classes")
            val savedConfigPath = Paths.get("${appHome}/resources/main/application.conf.saved")

            var config = Files.readString(Paths.get(savedConfigPath.toString()))

            val configOrig = Paths.get("${appHome}/resources/main/application.conf")
            //return the original config file
            configOrig.toFile().writeText(config.toString())

            // delete the saved config file
            savedConfigPath.toFile().delete()

        }
    }

    // We do not include a test that calls StartServer as it results in the software hanging while the
    // server is active.  JUnit tests for brAPI endpoints are executed using the Ktor testApplication
    // which doesn't require the server to be running.

    @Test
    fun testGetDbPathFromConfigFile() {
        // Running in junit gives a slightly different path then running from the command line
        // Because of this, we need to substring appHome to remove parts of the path that
        // include the "classes" folder and those below it.

        val appHome = StartServer().getClassPath().substringBefore("classes")
        assertTrue(appHome != null) // will be different based on which user is running the test
        println("appHome: $appHome")
        val configPath = StartServer().getDbPathFromConfigFile(appHome) // this should be consistent
        assertTrue(configPath == null) // null because we have not yet written "TILEDB_URI" line to the application.conf file
        //println("configPath: $configPath")


        // cleanup - should be done in @AfterAll, but failed when we ran the full suite of tests
        // vs running them individually
        val savedConfigPath = Paths.get("${appHome}/resources/main/application.conf.saved")
        var config = Files.readString(Paths.get(savedConfigPath.toString()))

        val configOrig = Paths.get("${appHome}/resources/main/application.conf")
        //return the original config file
        configOrig.toFile().writeText(config.toString())
    }

    @Test
    fun testUpdateConfigFile() {
        // The "beforeAll" saves the original config file, then returns it after the test
        val dbPath = "/Users/lcj34/temp/phgv2Tests/tempDir/testTileDBURI/"
        val port = 8090
        val appHome = StartServer().getClassPath().substringBefore("classes")
        assertTrue(appHome != null)
        println("appHome: $appHome")

        // update the config file with dbPath and port values above
        StartServer().updateConfigFile(dbPath,port.toString(),appHome)

        // Verify the TILEDB_URI and PORT values in the config file
        val tiledbURI = StartServer().getDbPathFromConfigFile(appHome)
        println("tiledbURI: $tiledbURI")
        assertTrue(tiledbURI == dbPath)

        // check the PORT value
        val configPath = Paths.get("${appHome}/resources/main/application.conf")

        var configLines = Files.readString(Paths.get(configPath.toString()))
        // Find the line in the config file that starts with "TILEDB_URI"
        // and extract the path
        // This will return null if the line is not found
        val portFromFile = configLines.lines().find { it.startsWith("PORT") }?.substringAfter("=")?.substringBefore("\n")

        assertTrue(portFromFile != null)
        assertEquals(port,portFromFile.toInt())

        // cleanup
        val savedConfigPath = Paths.get("${appHome}/resources/main/application.conf.saved")
        var config = Files.readString(Paths.get(savedConfigPath.toString()))

        val configOrig = Paths.get("${appHome}/resources/main/application.conf")
        //return the original config file
        configOrig.toFile().writeText(config.toString())

    }
}
