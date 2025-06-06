package net.maizegenetics.phgv2.brapi

import io.ktor.serialization.kotlinx.json.*
import io.ktor.server.application.*
import io.ktor.server.plugins.calllogging.CallLogging
import io.ktor.server.plugins.contentnegotiation.*
import io.ktor.server.plugins.defaultheaders.*
import io.ktor.server.routing.*
import kotlinx.serialization.json.Json
import net.maizegenetics.phgv2.brapi.api.apiRoute

/**
 * This is the main method that starts the PHG ktor web service.
 * It uses a Ktor HOCON application configuration file to get web and tiledb
 * location information.
 * That file is  src/main/kotlin/resources/application.conf
 *
 *   @author lcj34
 */

fun Application.module() {
    install(DefaultHeaders)
    install(CallLogging)
    // install(WebSockets)

    install(ContentNegotiation) {
        json(Json {
            prettyPrint = false
            isLenient = true
            encodeDefaults = true
        })
    }

    // Setup routing.  Individual endpoints create Kotlin Route extensions
    // to handle processing REST requests.

    routing {
        // this method routes brapi/v2/
        // Within apiRoute(), specific endpoint calls are handled
        apiRoute()
    }

}