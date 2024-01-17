package net.maizegenetics.phgv2.brapi.api

import com.typesafe.config.ConfigFactory
import io.ktor.http.*
import io.ktor.server.application.*
import io.ktor.server.config.*
import io.ktor.server.response.*
import io.ktor.server.routing.*
import net.maizegenetics.phgv2.brapi.model.MetadataTokenPagination
import net.maizegenetics.phgv2.brapi.model.VariantsListResponse
import net.maizegenetics.phgv2.brapi.model.VariantsListResponseResult
import net.maizegenetics.phgv2.brapi.service.VariantsService

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

    }
}