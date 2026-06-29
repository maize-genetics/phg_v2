package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File

class ImputePathFromPs4gTest {
    companion object {
        val tempTestDir = "${TestExtension.tempDir}ropebwtTest/"


        //Setup/download  files
        //Resetting on both setup and teardown just to be safe.
        @JvmStatic
        @BeforeAll
        fun setup() {
            resetDirs()
            setupDebugLogging()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            //resetDirs()
        }

        private fun resetDirs() {

            File(TestExtension.tempDir).deleteRecursively()
            File(TestExtension.testOutputFastaDir).deleteRecursively()
            File(TestExtension.testOutputDir).deleteRecursively()
            File(tempTestDir).deleteRecursively()

            File(TestExtension.tempDir).mkdirs()
            File(TestExtension.testOutputFastaDir).mkdirs()
            File(TestExtension.testOutputDir).mkdirs()
            File(tempTestDir).mkdirs()
        }
    }

    @Test
    fun testCliktParams() {

    }

    @Test
    fun testImputeHaploidPath() {
        val testArgs = "--read-file /Users/pjbra/temp/CML247XOh43-sr-simulated5000.ps4g --out-path-dir ${TestExtension.testOutputDir} " +
                "--prob-same 0.99 --prob-correct 0.99 --bin-size 5000"
        val imputeResult = ImputePathFromPs4g().test(testArgs)

    }

    @Test
    fun testDiploidPath() {
        val testArgs = "--read-file /Users/pjbra/temp/CML247XOh43-sr-simulated5000.ps4g --out-path-dir ${TestExtension.testOutputDir} " +
                "--prob-same 0.9999 --path-type diploid --prob-correct 0.99 --inbreed-coef 0.5 --n-parents 4 --bin-size 5000"
        val imputeResult = ImputePathFromPs4g().test(testArgs)

    }
}

