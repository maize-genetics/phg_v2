package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File

class SetupEnvironmentTest {

    companion object {
        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"

        @JvmStatic
        @BeforeAll
        fun setup() {
            File(tempDir).mkdirs()
        }

        // Comment these out if you need to look at the logs files created by ProcessBuilder
        // comda commands.
        @JvmStatic
        @AfterAll
        fun teardown() {
            File(tempDir).deleteRecursively()
        }
    }

    @Test
    fun testSetupEnvironment() {

        // When running createEnvironment, debug is printed to the screen.
        // IF an error occurs, the user is instructed on checking the log files,
        // and to verify conda is installed.
        // If the environment already exists, this test will pass with a message
        // indicating the environment already exists.
        // If the environment doesn't exist and the test passes, the user will
        // receive a message indicating the environment was created successfully.

        val envFile = "src/test/resources/net/maizegenetics/phgv2/cli/test.yml"

        val setupEnv = SetupEnvironment()
        val result = setupEnv.createEnvironment( envFile,tempDir)
        println("result = $result")
         // To verify, we  check that the outputDir contains the expected files
        val expectedFiles = listOf("condaCreate_error.log", "condaCreate_output.log")
        expectedFiles.forEach { assert(File(tempDir + it).exists()) }
    }

    @Test
    fun testSetupEnvironmentNoEnvFile() {
        // This test verifies the default environment file is used when no user defined
        // environment file is provided.

        val envFile = "" // no user defined file - should use the default file

        val setupEnv = SetupEnvironment()

        val result = setupEnv.createEnvironment(envFile, tempDir)
        println("result = $result")
        // To verify, we  check that the outputDir contains the expected files
        val expectedFiles = listOf("condaCreate_error.log", "condaCreate_output.log")
        expectedFiles.forEach { assert(File(tempDir + it).exists()) }
    }
}