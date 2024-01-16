package net.maizegenetics.phgv2.brapi.api

import com.fasterxml.jackson.module.kotlin.jacksonObjectMapper
import io.kotest.matchers.shouldBe
import io.ktor.client.call.*
import io.ktor.client.plugins.contentnegotiation.*
import io.ktor.client.request.*
import io.ktor.http.*
import io.ktor.serialization.kotlinx.json.*
import io.ktor.server.testing.*
import model.VariantsListResponse
import net.maizegenetics.phgv2.brapi.createSmallSeqTiledb
import net.maizegenetics.phgv2.brapi.model.SampleListResponse
import net.maizegenetics.phgv2.brapi.resetDirs
import net.maizegenetics.phgv2.brapi.service.VariantsService
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File

class VariantsTest {
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
    fun testCreateRefRangesForCache() {
        // THis test will verify the VarinetService.createRefRangesForCache() method
        // the createSmallSeqTiledb() should have created the tileDB folder, and copied
        // the bed file to it.  Use that for testing
        val bedFile = File("${TestExtension.testTileDBURI}/reference/").walk().filter { it.name.endsWith(".bed") }.toList()[0]
        val referenceRanges = VariantsService().createRefRangesForCache(bedFile)
        // Verify the number of lines in the file bedFile matches the number of reference ranges
        // in the referenceRanges list
        val bedFileLines = bedFile.readLines()
        assert(bedFileLines.size == referenceRanges.size)
    }

    @Test
    fun testVariantsNoBedFile() = testApplication {

        // This is needed, or you get "NoTransformationFoundException" from ktor HttpClient
        val client = createClient {
            install(ContentNegotiation) {
                json()
            }
        }

        // copy the bedFile created by createSMallSeqTiledb() to a file named the same but with ".save" appended
        // This will allow us to delete the bedFile and then restore it after the test
        val bedFile = File("${TestExtension.testTileDBURI}/reference/").walk().filter { it.name.endsWith(".bed") }.toList()[0]
        val bedFileSave = File("${TestExtension.testTileDBURI}/reference/${bedFile.name}.save")
        bedFile.copyTo(bedFileSave)

        // Delete the bedFile created by createSmallSeqTiledb()
        bedFile.delete()

        // Run test to verify that the server returns no data
        val response = client.get("/brapi/v2/variants")
        Assertions.assertEquals(HttpStatusCode.OK, response.status)
        val variants = response.body<VariantsListResponse>().result
        println("variants: $variants")
        Assertions.assertEquals(0, variants.data.size)

        //Restore the bedFile
        bedFileSave.copyTo(bedFile)

    }
    @Test
    fun testVariantsDefaultPageSize() = testApplication {

        // This is needed, or you get "NoTransformationFoundException" from ktor HttpClient
        val client = createClient {
            install(ContentNegotiation) {
                json()
            }
        }

        val response = client.get("/brapi/v2/variants")
        Assertions.assertEquals(HttpStatusCode.OK, response.status)
        val variants = response.body<VariantsListResponse>().result
        println("variants: $variants")
        Assertions.assertEquals(40, variants.data.size)
        val variantEntries = variants.data
        Assertions.assertEquals(1, variantEntries[0].start)
        Assertions.assertEquals(1001, variantEntries[1].start)
        Assertions.assertEquals(5501, variantEntries[2].start)

    }

    @Test
    fun testVariantsPageSize3() = testApplication {

        // This is needed, or you get "NoTransformationFoundException" from ktor HttpClient
        val client = createClient {
            install(ContentNegotiation) {
                json()
            }
        }

        val response = client.get("/brapi/v2/variants?pageSize=3")
        Assertions.assertEquals(HttpStatusCode.OK, response.status)

        val variants = response.body<VariantsListResponse>().result
        println("variants: $variants")
        Assertions.assertEquals(3, variants.data.size)
        val variantEntries = variants.data

        Assertions.assertEquals(1, variantEntries[0].start)
        Assertions.assertEquals(1001, variantEntries[1].start)
        Assertions.assertEquals(5501, variantEntries[2].start)

        val pagination = response.body<VariantsListResponse>().metadata.pagination
        val nextPageToken = pagination!!.nextPageToken
        nextPageToken shouldBe "4"

    }
}