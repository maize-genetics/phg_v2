package net.maizegenetics.phgv2.brapi.api

import io.ktor.server.http.content.*
import io.ktor.server.routing.*
import net.maizegenetics.phgv2.brapi.service.VariantSetsService
import java.io.File

/**
 * Method handles all REST messages coming to the BrAPI interface.  An entry
 * should be added to the "route" block for each endpoint which is supported.
 */
fun Routing.apiRoute() {

    route("/brapi/v2") {
        // this one is for testing
        baseRoute()

        // Add an api function for each endpoint, which should call a corresponding service
        // to process data for the endpoints.
        serverInfo()

        samples()

        // calls()
        // callSets()
        // referenceSets()
        // references()
        // samples()
        // studies()

        // alleleMatrix()

        variantSets()
        variants()

        // variantTables()

        staticFiles("/${VariantSetsService.allSamplesFileName}", File(VariantSetsService.allSamplesHvcf))

    }

}