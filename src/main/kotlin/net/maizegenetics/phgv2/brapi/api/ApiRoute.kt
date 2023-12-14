package net.maizegenetics.phgv2.brapi.api

import com.typesafe.config.ConfigFactory
import io.ktor.server.config.*
import io.ktor.server.routing.*

private val config = HoconApplicationConfig(ConfigFactory.load())

val defaultCallPageSize = config.property("callsPageSize").getString().toInt()
val defaultVariantsPageSize = config.property("variantsPageSize").getString().toInt()

/**
 * Method handles all REST messages coming to the BrAPI interface.  An entry
 * should be added to the "route" block for each endpoint which is supported.
 */
fun Routing.apiRoute() {

    route("/brapi/v2") {

        // Add an api function for each endpoint, which should call a corresponding service
        // to process data for the endpoints.
        serverInfo()

        // calls()
        // callSets()
        // referenceSets()
        // references()
        // samples()
        // studies()

        // alleleMatrix()

        // variantSets()
        // variants()

        // variantTables()

    }

}