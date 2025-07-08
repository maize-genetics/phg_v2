package net.maizegenetics.phgv2.brapi.api

import io.ktor.client.call.*
import io.ktor.client.plugins.contentnegotiation.*
import io.ktor.client.request.*
import io.ktor.http.*
import io.ktor.serialization.kotlinx.json.*
import io.ktor.server.application.*
import io.ktor.server.plugins.calllogging.*
import io.ktor.server.routing.*
import io.ktor.server.testing.*
import net.maizegenetics.phgv2.brapi.model.ServerInfoResponse
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test
import org.slf4j.event.Level
import io.ktor.server.plugins.contentnegotiation.ContentNegotiation as ServerContentNegotiation

/**
 * This class tests the Brapi ServerInfo endpoint.
 * It makes use of the ktor testApplication method to test the server.
 * No need to have an external server running - this bypasses the server and
 * executes the tests.
 */

class ServerInfoTest {
    companion object {
        // Don't need anything here for now
    }


    @Test
    fun testServerInfo() = testApplication {

        application {
            this@application.install(CallLogging) {
                level = Level.DEBUG
            }
            this@application.install(ServerContentNegotiation) {
                json()
            }
            routing {
                apiRoute()
            }
        }

        // This is needed or you get "NoTransformationFoundException" from ktor HttpClient
        val client = createClient {
            install(ContentNegotiation) {
                json()
            }
        }

        val response = client.get("/brapi/v2/serverinfo")
        assertEquals(HttpStatusCode.OK, response.status)
        val serverInfo = response.body<ServerInfoResponse>().result
        println("serverInfo: $serverInfo")

        // This assert will need to change when more endpoints are supported
        assertEquals(4, serverInfo.calls.size)
        val calls = serverInfo.calls
        assertEquals("/samples", calls[0].service)
    }

}