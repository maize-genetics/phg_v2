package net.maizegenetics.phgv2.cli

import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

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
        val result = setupEnv.createEnvironment(envFile, tempDir,"")
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

        val result = setupEnv.createEnvironment(envFile, tempDir,"")
        println("result = $result")
        // To verify, we  check that the outputDir contains the expected files
        val expectedFiles = listOf("condaCreate_error.log", "condaCreate_output.log")
        expectedFiles.forEach { assert(File(tempDir + it).exists()) }
    }

    @Test
    fun testEnvFileDoesNotExist() {

        val envFile = "src/test/resources/net/maizegenetics/phgv2/cli/test_does_not_exist.yml"

        val setupEnv = SetupEnvironment()
        val exceptionThrow = try {
            setupEnv.createEnvironment(envFile, tempDir,"")
            false
        } catch (e: IllegalStateException) {
            assert(e.message!!.contains("file $envFile does not exist"))
            true
        }
        assert(exceptionThrow)

    }

    @Test
    fun testExistingEnvironmentFailureDetectedFromStderr() {
        val isExistingEnvFailure = SetupEnvironment.isExistingEnvironmentFailure(
            listOf("CondaValueError: prefix already exists: /tmp/phgv2-conda"),
            emptyList()
        )

        assertTrue(isExistingEnvFailure)
    }

    @Test
    fun testExistingEnvironmentFailureIsFalseForOtherCondaErrors() {
        val isExistingEnvFailure = SetupEnvironment.isExistingEnvironmentFailure(
            listOf("PackagesNotFoundError: The following packages are not available from current channels:"),
            emptyList()
        )

        assertFalse(isExistingEnvFailure)
    }

    @Test
    fun testSummarizeCondaFailurePrefersStderr() {
        val summary = SetupEnvironment.summarizeCondaFailure(
            listOf("PackagesNotFoundError: anchorwave=1.2.5"),
            listOf("Collecting package metadata"),
            "/tmp/condaCreate_error.log"
        )

        assertEquals("PackagesNotFoundError: anchorwave=1.2.5", summary)
    }

    @Test
    fun testTermsOfServiceFailureDetectedFromStderr() {
        val isTermsOfServiceFailure = SetupEnvironment.isTermsOfServiceFailure(
            listOf("CondaToSNonInteractiveError: Terms of Service have not been accepted"),
            emptyList()
        )

        assertTrue(isTermsOfServiceFailure)
    }

    @Test
    fun testTermsOfServiceHelpMessageContainsExpectedCommands() {
        val helpMessage = SetupEnvironment.termsOfServiceHelpMessage()

        assertTrue(helpMessage.contains("conda tos accept"))
        assertTrue(helpMessage.contains("https://repo.anaconda.com/pkgs/main"))
        assertTrue(helpMessage.contains("CONDA_PLUGINS_AUTO_ACCEPT_TOS=yes"))
    }

    @Test
    fun testCondaEnvironmentNamesMatchExpectedValues() {
        assertEquals(
            "phgv2-conda",
            SetupEnvironment.condaEnvironmentName(SetupEnvironment.CONDAENVTYPE.PHG)
        )
        assertEquals(
            "phgv2-tiledb",
            SetupEnvironment.condaEnvironmentName(SetupEnvironment.CONDAENVTYPE.TILEDB)
        )
    }

}