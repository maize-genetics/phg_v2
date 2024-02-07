package net.maizegenetics.phgv2.brapi.utilities

import com.typesafe.config.ConfigFactory
import io.ktor.server.config.*
import org.apache.logging.log4j.LogManager

object BrAPIConfig {

    private val myLogger = LogManager.getLogger(BrAPIConfig::class.java)

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

        myLogger.info("BrAPIConfig:  TILEDB_URI: ${tiledbURI}")
        myLogger.info("BrAPIConfig:  contactEmail: ${contactEmail}")
        myLogger.info("BrAPIConfig:  documentationURL: ${documentationURL}")
        myLogger.info("BrAPIConfig:  location: ${location}")
        myLogger.info("BrAPIConfig:  organizationName: ${organizationName}")
        myLogger.info("BrAPIConfig:  organizationURL: ${organizationURL}")
        myLogger.info("BrAPIConfig:  serverDescription: ${serverDescription}")
        myLogger.info("BrAPIConfig:  serverName: ${serverName}")
        myLogger.info("BrAPIConfig:  callsPageSize: ${defaultCallPageSize}")
        myLogger.info("BrAPIConfig:  variantsPageSize: ${defaultVariantsPageSize}")

    }

}