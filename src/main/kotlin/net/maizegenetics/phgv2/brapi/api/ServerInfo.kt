package net.maizegenetics.phgv2.brapi.api

/**
 * This file returns information containing the list of brAPI endpoints supported by this server.
 */
import com.typesafe.config.ConfigFactory
import io.ktor.server.application.*
import io.ktor.server.config.*
import io.ktor.server.response.*
import io.ktor.server.routing.*
import net.maizegenetics.phgv2.brapi.model.*
import model.Metadata


private val config = HoconApplicationConfig(ConfigFactory.load())

val contactEmail = config.property("contactEmail").getString()

val documentationURL = config.property("documentationURL").getString()

val location = config.property("location").getString()

val organizationName = config.property("organizationName").getString()

val organizationURL = config.property("organizationURL").getString()

val serverDescription = config.property("serverDescription").getString()

val serverName = config.property("serverName").getString()


fun Route.serverInfo() {

    get("/serverinfo") {

        val endpointCalls = mutableListOf(

//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/allelematrix",
//                listOf(VersionsEnum._2)
//            ),

            Service(
                listOf(WSMIMEDataTypes.APPLICATION_JSON),
                listOf(MethodsEnum.GET),
                "/samples",
                listOf(VersionsEnum._0, VersionsEnum._1, VersionsEnum._2)
            ),
            Service(
                listOf(WSMIMEDataTypes.APPLICATION_JSON),
                listOf(MethodsEnum.GET),
                "/samples/{sampleDbId}",
                listOf(VersionsEnum._0, VersionsEnum._1, VersionsEnum._2)
            ),

            Service(
                listOf(WSMIMEDataTypes.APPLICATION_JSON),
                listOf(MethodsEnum.GET),
                "/serverinfo",
                listOf(VersionsEnum._0, VersionsEnum._1, VersionsEnum._2)
            ),
            Service(
                listOf(WSMIMEDataTypes.APPLICATION_JSON),
                listOf(MethodsEnum.GET),
                "/variants",
                listOf(VersionsEnum._0, VersionsEnum._1, VersionsEnum._2)
            ),

            
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/calls",
//                listOf(VersionsEnum._0, VersionsEnum._1, VersionsEnum._2)
//            ),
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/calls/{id}",
//                listOf(VersionsEnum._0, VersionsEnum._1, VersionsEnum._2)
//            ),
//
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/callsets",
//                listOf(VersionsEnum._0, VersionsEnum._1, VersionsEnum._2)
//            ),
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/callsets/{callSetDbId}",
//                listOf(VersionsEnum._0, VersionsEnum._1, VersionsEnum._2)
//            ),
//
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/references",
//                listOf(VersionsEnum._0, VersionsEnum._1, VersionsEnum._2)
//            ),
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/references/{referenceDbId}",
//                listOf(VersionsEnum._0, VersionsEnum._1, VersionsEnum._2)
//            ),
//
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/referencesets",
//                listOf(VersionsEnum._0, VersionsEnum._1, VersionsEnum._2)
//            ),
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/referencesets/{referenceSetDbId}",
//                listOf(VersionsEnum._0, VersionsEnum._1, VersionsEnum._2)
//            ),
//
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/studies",
//                listOf(VersionsEnum._0, VersionsEnum._1, VersionsEnum._2)
//            ),
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/studies/{studyDbId}",
//                listOf(VersionsEnum._0, VersionsEnum._1, VersionsEnum._2)
//            ),
//
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/variantTables",
//                listOf(VersionsEnum._1, VersionsEnum._2)
//            ),
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/variantTables/{variantTableDbId}",
//                listOf(VersionsEnum._1, VersionsEnum._2)
//            ),
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/variantTables/{variantTableDbId}/variants",
//                listOf(VersionsEnum._1, VersionsEnum._2)
//            ),
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/variantTables/{variantTableDbId}/samples",
//                listOf(VersionsEnum._1, VersionsEnum._2)
//            ),
//            Service(
//                listOf(WSMIMEDataTypes.APPLICATION_JSON),
//                listOf(MethodsEnum.GET),
//                "/variantTables/{variantTableDbId}/table",
//                listOf(VersionsEnum._1, VersionsEnum._2)
//            )

        )


        val serverInfo = ServerInfo(
            endpointCalls, contactEmail, documentationURL, location, organizationName,
            organizationURL, serverDescription, serverName
        )
        call.respond(ServerInfoResponse(Metadata(), serverInfo))

    }

}