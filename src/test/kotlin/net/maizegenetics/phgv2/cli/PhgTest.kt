package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class PhgTest {

    private val myLogger = LogManager.getLogger(PhgTest::class.java)

    @Test
    fun testPhgCommand() {
        val result = Phg().test(
            "--help"
        )

        myLogger.info("testPhgCommand: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")
    }

    @Test
    fun testPhgVersionCommand() {
        val result = Phg().test(
            "--version"
        )

        myLogger.info("testPhgVersionCommand: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")
        assert(result.output.startsWith("phg version 2")) { "version not found in output" }
    }

}