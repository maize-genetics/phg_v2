package net.maizegenetics.phgv2.utils

import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Assertions
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import org.junit.jupiter.api.extension.ExtendWith
import org.junit.jupiter.api.io.TempDir
import java.io.File
import java.io.FileNotFoundException
import java.nio.file.Path
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class GeneralUtilitiesTest {

    /**
     * Helper function to create a temporary version.properties file for testing
     */
    private fun createVersionPropertiesFile(directory: Path, content: String) {
        val versionFile = File(directory.toFile(), "version.properties")
        versionFile.writeText(content)
    }

    @Test
    fun testWorkingVersionProperties(@TempDir tempDir: Path) {
        // Simulate the file in the current working directory
        createVersionPropertiesFile(tempDir, """
            majorVersion=2
            minorVersion=3
            patchVersion=4
            buildNumber=234
        """.trimIndent())

        System.setProperty("user.dir", tempDir.toFile().absolutePath)

        val version = phgVersion()

        assertEquals("2.3.4.234", version)
    }

    @Test
    fun testMalformedVersionProperties(@TempDir tempDir: Path) {
        // Create a malformed version.properties file
        createVersionPropertiesFile(tempDir, """
            majorVersion=2
            minorVersion=wrongFormat
            patchVersion=4
            buildNumber=234
        """.trimIndent())

        System.setProperty("user.dir", tempDir.toFile().absolutePath)

        val exception = assertThrows<NumberFormatException> {
            phgVersion()
        }

        Assertions.assertNotNull(exception)
    }

    @Test
    fun testMissingVersionPropertiesFile() {
        // Simulate the absence of the version.properties file in both locations
        System.setProperty("user.dir", "/non/existent/directory")

        val exception = assertThrows<FileNotFoundException> {
            phgVersion()
        }

        assertEquals("/non/existent/directory/version.properties (No such file or directory)", exception.message)
    }

}