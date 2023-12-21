package net.maizegenetics.phgv2.brapi.api

import io.kotest.common.runBlocking
import io.kotest.core.spec.style.AnnotationSpec
import io.ktor.client.call.*
import io.ktor.client.plugins.contentnegotiation.*
import io.ktor.client.statement.*
import io.ktor.server.testing.*
import net.maizegenetics.phgv2.brapi.createSmallSeqTiledb
import net.maizegenetics.phgv2.brapi.resetDirs
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Assertions
import java.io.File

import io.ktor.client.request.get
import io.ktor.http.HttpStatusCode
import io.ktor.serialization.kotlinx.json.*
import net.maizegenetics.phgv2.brapi.model.SampleListResponse
import org.junit.jupiter.api.AfterAll
import kotlin.test.*
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.BeforeAll

/**
 * Test the brapi samples endpoint
 * This test uses the ktor test harness
 *
 * This works when the application.conf file has been modified to include
 * a value for TILEDB_URI.  Need to determine how to do this programmatically
 * for the test harness (but not when running the server)
 */
class SamplesTest {
    companion object {

        @JvmStatic
        @BeforeAll
        fun setup() {
            // delete, reset the directories
            resetDirs()

            // create the tiledb datasets, load them with from the vcf files
            // This will also create the AGC compressed file
            createSmallSeqTiledb()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testSamples() = testApplication {

        // This is needed or you get "NoTransformationFoundException" from ktor HttpClient
        val client = createClient {
            install(ContentNegotiation) {
                json()
            }
        }

        val response = client.get("/brapi/v2/samples")
        assertEquals(HttpStatusCode.OK, response.status)
        val samples = response.body<SampleListResponse>().result
        println("samples: $samples")
        assertEquals(3, samples.data.size)
        val sampleEntries = samples.data
        assertEquals("LineA", sampleEntries[0].sampleName)
        assertEquals("LineB", sampleEntries[1].sampleName)
        assertEquals("Ref", sampleEntries[2].sampleName)

    }

}