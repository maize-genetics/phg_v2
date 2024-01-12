package net.maizegenetics.phgv2.brapi.api

import com.typesafe.config.ConfigFactory
import io.ktor.http.*
import io.ktor.server.application.*
import io.ktor.server.config.*
import io.ktor.server.response.*
import io.ktor.server.routing.*
import model.Metadata
import model.MetadataTokenPagination
import model.VariantsListResponse
import model.VariantsListResponseResult
import net.maizegenetics.phgv2.brapi.model.SampleListResponse
import net.maizegenetics.phgv2.brapi.model.SampleListResponseResult
import net.maizegenetics.phgv2.brapi.service.VariantsService

private val config = HoconApplicationConfig(ConfigFactory.load())

fun Route.variants() {
    route("/variants") {
        val variantService = VariantsService()

        //Endpoint to get all the variants in the DB.
        //This can be used to get the list of all the ids to further get more information using the other endpoints.
        get("") {
            val pageToken =
                call.parameters["page"]?.toInt() ?: 1 // token is reference range ID, they are stored in order
            val pageSize = call.parameters["pageSize"]?.toInt() ?: defaultVariantsPageSize

            // These are tokens, not array indices.  The ranges are pulled from the database based on
            // reference_range_id BETWEEN start and end.  0 as a start won't cause an error, but it also
            // isn't truly valid.
            if (pageToken < 1) {
                call.respond(" ${HttpStatusCode.NotFound}: Invalid page token: ${pageToken}. The page token represents a reference range id, which is 1 or greater.")
            }

            val variantsAndPagination = variantService.generateVariantsListFromCache(pageToken, pageSize,"all")

            val pagination = variantsAndPagination.first
            val variants = variantsAndPagination.second
            var metadata = MetadataTokenPagination(pagination = pagination)
            call.respond(VariantsListResponse(metadata, VariantsListResponseResult(variants.toTypedArray())))
        }

    }
}