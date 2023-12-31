package net.maizegenetics.phgv2.brapi.api

import io.ktor.http.*
import io.ktor.server.application.*
import io.ktor.server.response.*
import io.ktor.server.routing.*
import model.Metadata
import net.maizegenetics.phgv2.brapi.model.SampleListResponse
import net.maizegenetics.phgv2.brapi.model.SampleListResponseResult
import net.maizegenetics.phgv2.brapi.service.SamplesService

fun Route.samples() {

    val samplesService = SamplesService

    route("/samples") {

        get("") {
            call.respond(
                SampleListResponse(
                    Metadata(),
                    SampleListResponseResult(samplesService.allTaxaNames().toTypedArray())
                )
            )
        }

        get("/{sampleDbId}") {
            val id = call.parameters["sampleDbId"] ?: throw IllegalStateException("Must provide id")
            val sample = samplesService.taxa(id)

            if (sample == null) {
                call.respond(" ${HttpStatusCode.NotFound}: The requested object $id was not found in the database")
            } else {
                call.respond(SampleListResponse(Metadata(), SampleListResponseResult(arrayOf(sample))))
            }
        }

    }
}