package net.maizegenetics.phgv2.brapi.utilities

import com.typesafe.config.ConfigFactory
import io.ktor.server.config.*

object BrAPIConfig {

    val tiledbURI: String
    val contactEmail: String
    val documentationURL: String
    val location: String
    val organizationName: String
    val organizationURL: String
    val serverDescription: String
    val serverName: String
    val defaultCallPageSize: Int
    val defaultVariantsPageSize: Int

    init {
        val config = HoconApplicationConfig(ConfigFactory.load())
        tiledbURI = config.property("TILEDB_URI").getString()
        contactEmail = config.property("contactEmail").getString()
        documentationURL = config.property("documentationURL").getString()
        location = config.property("location").getString()
        organizationName = config.property("organizationName").getString()
        organizationURL = config.property("organizationURL").getString()
        serverDescription = config.property("serverDescription").getString()
        serverName = config.property("serverName").getString()
        defaultCallPageSize = config.property("callsPageSize").getString().toInt()
        defaultVariantsPageSize = config.property("variantsPageSize").getString().toInt()
    }

}