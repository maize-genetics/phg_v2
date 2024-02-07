package net.maizegenetics.phgv2.brapi.api

import io.ktor.http.*
import io.ktor.server.application.*
import io.ktor.server.response.*
import io.ktor.server.routing.*
import net.maizegenetics.phgv2.brapi.model.Metadata
import net.maizegenetics.phgv2.brapi.model.VariantSetResponse
import net.maizegenetics.phgv2.brapi.service.VariantSetsService
import net.maizegenetics.phgv2.brapi.utilities.BrAPIConfig

/**
 * This method handles the "variantsets" endpoint.  It is used to get a list of all the variant sets in the database.
 * We are using this to return a URL that points to the multi-sample HVCF file.
 */

fun Route.variantSets() {

    route("/variantsets") {

        val variantSetService = VariantSetsService()
        get("") {
            // No pagination needed - we are only returning a URI that points to a file,
            // not the actual file contents (and not an alleleMatrix)
            val variantSet = variantSetService.generateVariantSets(BrAPIConfig.tiledbURI)
            var metadata = Metadata()
            call.respond(VariantSetResponse(metadata, variantSet))
        }

        get("/{variantSetId}") {
            call.response.status(HttpStatusCode.NotImplemented)
            // val id = call.parameters["variantSetId"] ?: throw IllegalStateException("Must provide id")
            // No pagination needed - we are only returning a URI that points to a file,
            // not the actual file contents (and not an alleleMatrix)
            // val variantSet = variantSetService.generateVariantSets(BrAPIConfig.tiledbURI, id)
            // var metadata = Metadata()
            // call.respond(VariantSetResponse(metadata, variantSet))
        }

    }

}