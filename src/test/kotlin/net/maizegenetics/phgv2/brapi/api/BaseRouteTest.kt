package net.maizegenetics.phgv2.brapi.api

import io.ktor.client.request.*
import io.ktor.client.statement.*
import io.ktor.http.*
import io.ktor.server.response.*
import io.ktor.server.routing.*
import io.ktor.server.testing.*
import org.junit.jupiter.api.Assertions
import org.junit.jupiter.api.Test

/**
 * This class is the initial test for the Brapi server.
 * It makes use of the ktor testApplication method to test the server.
 */
class BaseRouteTest {
    @Test
    fun testBaseRoute() = testApplication {

        application {
            routing {
                get("/brapi/v2/") {
                    call.respondText(
                        "Hello Brapi User! Please send your requests to /brapi/v2/",
                        ContentType.Text.Plain
                    )
                }
            }
        }

        // Example from:
        //  https://github.com/ktorio/ktor-documentation/blob/2.3.7/codeSnippets/snippets/engine-main/src/test/kotlin/EngineMainTest.kt
        val response = client.get("/brapi/v2/")
        Assertions.assertEquals(HttpStatusCode.OK, response.status)
        Assertions.assertEquals("Hello Brapi User! Please send your requests to /brapi/v2/", response.bodyAsText())
    }
}