package net.maizegenetics.phgv2.brapi.api

/**
 * This file returns information containing the list of brAPI endpoints supported by this server.
 */
import io.ktor.server.response.*
import io.ktor.server.routing.*
import net.maizegenetics.phgv2.brapi.model.*
import net.maizegenetics.phgv2.brapi.utilities.BrAPIConfig


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
            endpointCalls,
            BrAPIConfig.contactEmail,
            BrAPIConfig.documentationURL,
            BrAPIConfig.location,
            BrAPIConfig.organizationName,
            BrAPIConfig.organizationURL,
            BrAPIConfig.serverDescription,
            BrAPIConfig.serverName
        )
        call.respond(ServerInfoResponse(Metadata(), serverInfo))

    }

}