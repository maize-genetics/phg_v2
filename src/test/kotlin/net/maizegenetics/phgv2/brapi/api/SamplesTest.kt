package net.maizegenetics.phgv2.brapi.api

import io.ktor.client.call.*
import io.ktor.client.request.*
import io.ktor.http.*
import io.ktor.serialization.kotlinx.json.*
import io.ktor.server.routing.*
import io.ktor.server.testing.*
import net.maizegenetics.phgv2.brapi.createSmallSeqTiledb
import net.maizegenetics.phgv2.brapi.model.SampleListResponse
import net.maizegenetics.phgv2.brapi.resetDirs
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import io.ktor.client.plugins.contentnegotiation.ContentNegotiation as ClientContentNegotiation
import io.ktor.server.plugins.contentnegotiation.ContentNegotiation as ServerContentNegotiation

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

        application {
            routing {
                apiRoute()
            }
        }

        // This is needed, or you get "NoTransformationFoundException" from ktor HttpClient
        val client = createClient {
            install(ClientContentNegotiation) {
                json()
            }
        }

        val response = client.get("/brapi/v2/samples") {
            contentType(ContentType.Application.Json)
        }
        assertEquals(HttpStatusCode.OK, response.status)
        val samples = response.body<SampleListResponse>().result
        println("samples: $samples")
        assertEquals(3, samples.data.size)
        val sampleEntries = samples.data
        assertEquals("LineA", sampleEntries[0].sampleName)
        assertEquals("LineB", sampleEntries[1].sampleName)
        assertEquals("Ref", sampleEntries[2].sampleName)

    }

    @Test
    fun testSampleID() = testApplication {

        application {
            routing {
                apiRoute()
            }
        }

        // This is needed, or you get "NoTransformationFoundException" from ktor HttpClient
        val client = createClient {
            install(ClientContentNegotiation) {
                json()
            }
        }

        val response = client.get("/brapi/v2/samples/1")
        assertEquals(HttpStatusCode.OK, response.status)
        val samples = response.body<SampleListResponse>().result
        println("samples: $samples")
        assertEquals(1, samples.data.size)
        val sampleEntries = samples.data
        assertEquals("LineB", sampleEntries[0].sampleName)

    }

}