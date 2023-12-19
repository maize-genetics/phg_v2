package net.maizegenetics.phgv2.brapi.api

import io.ktor.server.application.*
import io.ktor.server.response.*
import io.ktor.server.routing.*

fun Route.baseRoute() {

    get("/") {
        call.respondText("Hello Brapi User! Please send your requests to /brapi/v2/")
    }
}


