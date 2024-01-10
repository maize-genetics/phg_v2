package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.brapi.createSmallSeqTiledb
import net.maizegenetics.phgv2.brapi.resetDirs
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import java.io.File
import java.nio.file.Files
import java.nio.file.Paths
import kotlin.test.Test
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
                if (!line.startsWith("TILEDB_URI")) {
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

    // This test should not be run as it results in the software hanging while the
    // server is active.
    @Test
    fun testStartServer() {

//        val dbUri = "/Users/lcj34/temp/phgv2Tests/tempDir/testTileDBURI/"
//        val command = StartServer()
//        val result = command.test("--db-path $dbUri")
//        println("result: $result")
//
//        // This doesn't get hit because we don't come back from the server call above
//        // can the server call be put into a ProcessBuilder command?
//        println("Now calling stopServer - will this work?")
//        //StartServer().stopServer()

    }

    @Test
    fun testGetDbPathFromConfigFile() {
        // Running in junit gives a slightly different path then running from the command line
        // Because of this, we need to substring appHome to remove parts of the path that
        // include the "classes" folder and those below it.

        val appHome1 = StartServer().getClassPath()
        println("appHome1 - before substring: $appHome1") // used for debugging
        val appHome = StartServer().getClassPath().substringBefore("classes")
        assertTrue(appHome != null) // will be different based on which user is running the test
        println("appHome: $appHome")
        val configPath = StartServer().getDbPathFromConfigFile(appHome) // this should be consistent
        assertTrue(configPath == null) // null because we have not yet written "TILEDB_URI" line to the application.conf file
        //println("configPath: $configPath")
    }

    @Test
    fun testWriteConfigFile() {
        val dbUri = "/Users/lcj34/temp/phgv2Tests/tempDir/testTileDBURI/"
        val appHome = StartServer().getClassPath().substringBefore("classes")
        assertTrue(appHome != null)
        println("appHome: $appHome")
        StartServer().writeConfigFile(dbUri,appHome)
        val configPath = StartServer().getDbPathFromConfigFile(appHome)
        assertTrue(configPath == dbUri)
        println("configPath: $configPath")
    }
}
