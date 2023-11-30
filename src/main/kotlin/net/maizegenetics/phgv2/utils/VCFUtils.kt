@file:JvmName("VCFUtils")

package net.maizegenetics.phgv2.utils

import htsjdk.variant.vcf.VCFHeader
import org.apache.logging.log4j.LogManager

private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.utils.VCFUtils")

// Making Number a string as VCF allows for '.'
data class AltHeaderMetaData(
    val id: String, val description: String, val number: String, val source: String,
    val contig: String, val start: Int, val end: Int, val checksum: String, val refRange: String
)

/**
 * Helper function to parse out the ALT headers from the VCF file.
 *
 * We need to do a bit more involved parsing in this function as we cannot use the .getOtherHeaders() call from HTSJDK.
 * For some reason this only returns the first header when called, and we need all of them.
 * The workaround is that we can get all the metadata, filter out any that are not ALT then parse the ALT header using normal string parsing.
 * To make this easy, we just parse each piece of metadata into a key-value pair and then store in a map.
 */
fun parseALTHeader(header: VCFHeader): Map<String, AltHeaderMetaData> {

    // Need to turn the ALT File header into a Map<ID, AltHeaderMetaData>
    return header.metaDataInInputOrder
        .filter { it.key == "ALT" }
        .map {
            it.toString()
                .substringAfter("<")
                .substringBeforeLast(">")
        } // Keep the useful part of the ALT Tag
        .map { it.split(",") }
        .associate {
            val idsToValueMap = it.map { token -> token.split("=") }.associate { token -> token[0] to token[1] }
            // ID and Description are required fields by VCF spec, if these errors are thrown there is something wrong with the htsjdk library
            check(idsToValueMap.containsKey("ID")) { "ALT Header does not contain ID" }
            check(idsToValueMap.containsKey("Description")) { "ALT Header does not contain Description" }
            // These are optional header fields, so we check these in the unit test.
            check(idsToValueMap.containsKey("Number")) { "ALT Header does not contain Number" }
            check(idsToValueMap.containsKey("Source")) { "ALT Header does not contain Source" }
            check(idsToValueMap.containsKey("Contig")) { "ALT Header does not contain Contig" }
            check(idsToValueMap.containsKey("Start")) { "ALT Header does not contain Start" }
            check(idsToValueMap.containsKey("End")) { "ALT Header does not contain End" }
            check(idsToValueMap.containsKey("Checksum")) { "ALT Header does not contain Checksum" }
            check(idsToValueMap.containsKey("RefRange")) { "ALT Header does not contain RefRange" }

            idsToValueMap["ID"]!! to AltHeaderMetaData(
                idsToValueMap["ID"]!!,
                idsToValueMap["Description"]!!,
                idsToValueMap["Number"]!!,
                idsToValueMap["Source"]!!,
                idsToValueMap["Contig"]!!,
                idsToValueMap["Start"]!!.toInt(),
                idsToValueMap["End"]!!.toInt(),
                idsToValueMap["Checksum"]!!,
                idsToValueMap["RefRange"]!!
            )
        }

}