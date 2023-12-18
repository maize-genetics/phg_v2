package net.maizegenetics.phgv2.brapi.api

import com.fasterxml.jackson.module.kotlin.jacksonObjectMapper
import com.fasterxml.jackson.module.kotlin.readValue
import net.maizegenetics.phgv2.brapi.createSmallSeqTiledb
import net.maizegenetics.phgv2.brapi.model.ServerInfoResponse
import net.maizegenetics.phgv2.brapi.resetDirs
import net.maizegenetics.phgv2.cli.StartServer
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import java.net.http.HttpResponse

import io.kotest.common.runBlocking
import io.kotest.matchers.shouldBe
import io.ktor.client.HttpClient
import io.ktor.client.call.*
import io.ktor.client.engine.cio.*
import io.ktor.client.plugins.contentnegotiation.*
import org.junit.jupiter.api.DisplayName
import io.ktor.client.request.get
import io.ktor.http.HttpStatusCode
import io.ktor.serialization.kotlinx.json.*
import org.junit.jupiter.api.Assertions.assertEquals


class ServerInfoTest {
    companion object {

        @JvmStatic
        @BeforeAll
        fun setup() {
            // This is not needed for ServerInfoTest
            // It is here as an example of what shoudl be put in the @BeforeAll for other tests
            // that must access a tiledb dataset
            resetDirs()

            // create the tiledb datasets, load them with from the vcf files
            createSmallSeqTiledb()

            // now start the server
            // This is currently not working.  You must start the server in a separate
            // terminal window before executing these tests
            //StartServer()

        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testServerInfo() {
        val query = "http://localhost:8080/brapi/v2/serverinfo"

        runBlocking {
            // This is needed or you get "NoTransformationFOundException" from ktor HttpClient
            val client = HttpClient(CIO) {
                install(ContentNegotiation) {
                    json()
                }
            }
            client.use { myClient ->
                val response = myClient.get("http://localhost:8080/brapi/v2/serverinfo")
                assertEquals(HttpStatusCode.OK, response.status)
                val serverInfo = response.body<ServerInfoResponse>().result
                println("serverInfo: $serverInfo")
                assertEquals(4, serverInfo.calls.size)
            }
        }
    }
}