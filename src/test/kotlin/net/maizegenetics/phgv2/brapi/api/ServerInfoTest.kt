package net.maizegenetics.phgv2.brapi.api

import net.maizegenetics.phgv2.brapi.model.ServerInfoResponse
import io.ktor.client.call.*
import io.ktor.client.plugins.contentnegotiation.*
import io.ktor.client.request.get
import io.ktor.http.HttpStatusCode
import io.ktor.serialization.kotlinx.json.*
import io.ktor.server.testing.*

import kotlin.test.*
import org.junit.jupiter.api.Assertions.assertEquals

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
        assertEquals(3, serverInfo.calls.size)
        val calls = serverInfo.calls
        assertEquals("/samples", calls[0].service)
    }

}