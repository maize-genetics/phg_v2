@file:JvmName("VCFUtils")

package net.maizegenetics.phgv2.utils

import htsjdk.variant.vcf.VCFAltHeaderLine
import htsjdk.variant.vcf.VCFHeader
import htsjdk.variant.vcf.VCFHeaderVersion
import net.maizegenetics.phgv2.api.SampleGamete
import org.apache.logging.log4j.LogManager

private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.utils.VCFUtils")

// Making Number a string as VCF allows for '.'
data class AltHeaderMetaData(
    val id: String, val description: String, val source: String, val sampleGamete: SampleGamete,
    val regions: List<Pair<Position, Position>>, val checksum: String, val refRange: String
) {
    fun sampleName() = sampleGamete.name
    fun gamete() = sampleGamete.gameteId
}

/**
 * Helper function to parse out the ALT headers from the VCF file.
 *
 * We need to do a bit more involved parsing in this function as we cannot use the .getOtherHeaders() call from HTSJDK.
 * For some reason this only returns the first header when called, and we need all of them.
 * The workaround is that we can get all the metadata, filter out any that are not ALT then parse the ALT header using normal string parsing.
 * To make this easy, we just parse each piece of metadata into a key-value pair and then store in a map.
 */
fun parseALTHeader(header: VCFHeader): Map<String, AltHeaderMetaData> {

    return header.metaDataInInputOrder.asSequence().filter { it.key == "ALT" }
        .map { it as VCFAltHeaderLine }
        .map { it.genericFields }
        .associateBy { it["ID"]!! }
        .map {
            check(it.value.containsKey("ID")) { "ALT Header does not contain ID" }
            check(it.value.containsKey("Description")) { "ALT Header does not contain Description" }
            // These are optional header fields, so we check these in the unit test.
            check(it.value.containsKey("Source")) { "ALT Header does not contain Source" }
            check(it.value.containsKey("SampleName")) { "ALT Header does not contain SampleName" }
            check(it.value.containsKey("Regions")) { "ALT Header does not contain Regions" }
            check(it.value.containsKey("Checksum")) { "ALT Header does not contain Checksum" }
            check(it.value.containsKey("RefRange")) { "ALT Header does not contain RefRange" }
            it.key to AltHeaderMetaData(
                it.value["ID"]!!,
                it.value["Description"]!!,
                it.value["Source"]!!,
                SampleGamete(it.value["SampleName"]!!, it.value["Gamete"]?.toInt() ?: 0),
                parseRegions(it.value["Regions"]!!),
                it.value["Checksum"]!!,
                it.value["RefRange"]!!
            )
        }.toList()
        .toMap()
}

/**
 * Function to parse the regions from the ALT header.
 */
fun parseRegions(regions: String): List<Pair<Position, Position>> {
    return regions.split(",").map { it.split(":") }.map {
        val positions = it[1].split("-").map { position -> position.toInt() }
        check(positions.size == 2) { "Region $it is not in the correct format.  It needs to be in the form: chr:stPos-endPos." }
        Pair(Position(it[0], positions[0]), Position(it[0], positions[1]))
    }
}

/**
 * Function to convert the ALT header metadata into a VCF header line.
 *
 * @param altHeaderData The ALT header metadata to convert.
 * @param altHeaderId The ID to use for the ALT header line.  If null, the ID from the ALT header metadata is used.
 */
fun altHeaderMetadataToVCFHeaderLine(altHeaderData: AltHeaderMetaData, altHeaderId: String? = null): VCFAltHeaderLine {

    return VCFAltHeaderLine(
        "<ID=${altHeaderId ?: altHeaderData.id}, " +
                "Description=\"haplotype data for line: ${altHeaderData.sampleName()}\">," +
                "Source=\"${altHeaderData.source}\",SampleName=\"${altHeaderData.sampleName()}\"," +
                "Regions=\"${altHeaderData.regions.joinToString(",") { "${it.first.contig}:${it.first.position}-${it.second.position}" }}\"," +
                "Checksum=\"Md5\",RefRange=\"${altHeaderData.refRange}\">",
        VCFHeaderVersion.VCF4_2
    )

}