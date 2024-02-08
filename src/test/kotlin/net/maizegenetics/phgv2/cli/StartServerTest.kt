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
    companion object {  // static variables

        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"
        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.tempDir).mkdirs()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(TestExtension.tempDir).deleteRecursively()
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

        val appHome = StartServer.getClassPath().substringBefore("classes")
        assertTrue(appHome != null) // will be different based on which user is running the test
        println("appHome: $appHome")

        // copy the file to a temporary location
        val testConfigPath = Paths.get("${TestExtension.tempDir}/resources/main/")
        val origConfigFile = Paths.get("${appHome}/resources/main/application.conf")
        // make the testConfigPath directory if it doesn't exist
        testConfigPath.toFile().mkdirs()
        // copy the origConfigFile to the testConfigPath
        // But omit any TILEDB_URI lines.
        // junit from Intellij will not have the TILEDB_URI line in the config file
        // But CI writes that line before it starts running tests.
        var config = ""
        Files.lines(Paths.get(origConfigFile.toString())).forEach { line ->
            if (!line.startsWith("TILEDB_URI") ) {
                config += line + "\n"
            }
        }
        val testConfigFile = Paths.get("${testConfigPath}/application.conf")
        testConfigFile.toFile().writeText(config)

        val tileDbPath = StartServer.getDbPathFromConfigFile(TestExtension.tempDir) // this should be consistent

        assertTrue(tileDbPath == null) // null because we have not yet written "TILEDB_URI" line to the application.conf file
        //println("configPath: $configPath")

        // cleanup
        // Delete the temp config file
        val configFile = Paths.get("${TestExtension.tempDir}/resources/main/application.conf")
        configFile.toFile().delete()

        // Delete the testConfigPath directory
        testConfigPath.toFile().deleteRecursively()
    }

    @Test
    fun testUpdateConfigFile() {
        // Get the original config file
        val dbPath = "/Users/lcj34/temp/phgv2Tests/tempDir/testTileDBURI/"
        val port = 8090
        val appHome = StartServer.getClassPath().substringBefore("classes")
        assertTrue(appHome != null)
        println("appHome: $appHome")

        // copy the file to a temporary location
        val testConfigPath = Paths.get("${TestExtension.tempDir}/resources/main/")
        val origConfigFile = Paths.get("${appHome}/resources/main/application.conf")
        // make the testConfigPath directory if it doesn't exist
        testConfigPath.toFile().mkdirs()
        // copy the origConfigFile to the testConfigPath
        origConfigFile.toFile().copyTo(Paths.get("${testConfigPath}/application.conf").toFile())

        // update the config file with dbPath and port values above
        StartServer.updateConfigFile(dbPath,port.toString(),TestExtension.tempDir)

        // Verify the TILEDB_URI and PORT values in the config file
        val tiledbURI = StartServer.getDbPathFromConfigFile(TestExtension.tempDir)
        println("tiledbURI: $tiledbURI")
        assertTrue(tiledbURI == dbPath)

        // check the PORT value
        val configPath = Paths.get("${TestExtension.tempDir}/resources/main/application.conf")

        var configLines = Files.readString(Paths.get(configPath.toString()))
        // Find the line in the config file that starts with "TILEDB_URI"
        // and extract the path
        // This will return null if the line is not found
        val portFromFile = configLines.lines().find { it.startsWith("PORT") }?.substringAfter("=")?.substringBefore("\n")

        assertTrue(portFromFile != null)
        assertEquals(port,portFromFile.toInt())

        // Delete the configPath file
        configPath.toFile().delete()

        // Delete the testConfigPath directory
        testConfigPath.toFile().deleteRecursively()

    }
}
