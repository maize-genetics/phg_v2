package net.maizegenetics.phgv2.brapi.api

import com.typesafe.config.ConfigFactory
import io.ktor.http.*
import io.ktor.server.application.*
import io.ktor.server.config.*
import io.ktor.server.response.*
import io.ktor.server.routing.*
import net.maizegenetics.phgv2.brapi.model.*
import net.maizegenetics.phgv2.brapi.service.VariantsService

/**
 * This method handles the "variants" endpoint.  It is used to get a list of all the variants in the database.
 * "variants" in the PHG tiledb database are reference ranges.
 */
private val config = HoconApplicationConfig(ConfigFactory.load())

fun Route.variants() {
    route("/variants") {
        val variantService = VariantsService()

        //Endpoint to get all the variants in the DB.
        //This can be used to get the list of all the ids to further get more information using the other endpoints.
        get("") {
            val pageToken =
                call.parameters["page"]?.toInt() ?: 1 // token is index into the list of reference ranges, which are stored in sorted order
            val pageSize = call.parameters["pageSize"]?.toInt() ?: defaultVariantsPageSize

            // These are array indices, but user will think of them as 1 based, not 0-based.
            // the software will handle this.
            if (pageToken < 1) {
                call.respond(" ${HttpStatusCode.NotFound}: Invalid page token: ${pageToken}. The page token must be greater than 0")
            }

            val variantsAndPagination = variantService.generateVariantsListFromCache(pageToken, pageSize,"all")

            val pagination = variantsAndPagination.first
            val variants = variantsAndPagination.second
            var metadata = MetadataTokenPagination(pagination = pagination)
            call.respond(VariantsListResponse(metadata, VariantsListResponseResult(variants.toTypedArray())))
        }

        //This end point will return data for a specific variant corresponding to a variantDbId.
        //If the id is not found a 404 will be thrown.
        get("/{variantDbId}") {
            val pageToken = call.parameters["page"]?.toInt() ?: 1 // page token is reference range id, which starts at 1
            val pageSize = call.parameters["pageSize"]?.toInt() ?: defaultVariantsPageSize

            if(pageToken < 1) {
                call.respond(" ${HttpStatusCode.NotFound}: The page ${pageToken} is invalid.  It must be 1 or greater.")
            }

            val variantDbId = call.parameters["variantDbId"] ?: throw IllegalStateException("Must provide variantDbId")
            val variant = variantService.generateVariantFromID(variantDbId, pageToken, pageSize, "all")

            if(variant == null) {
                call.respond("${HttpStatusCode.NotFound}: The requested Variant object ${variantDbId} was not found in the database")
            }
            else {
                // There will only be 1 page, so nextPageToken will be defaulted to null
                call.respond(VariantSingleResponse(MetadataTokenPagination(pagination= TokenPagination(currentPageToken="1")), variant))
            }
        }

    }
}