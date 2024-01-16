package net.maizegenetics.phgv2.brapi.model

import kotlinx.serialization.Serializable

/**
 * A `VariantSet` represents the a a Variant Table.  It holds the dimensions of the table.
 * @param additionalInfo Additional arbitrary info
 * @param analysis The list of analysis which has been done to this VariantSet
 * @param availableFormats The list of VariantSetAvailable Formats which this VariantSet can be returned as
 * @param callSetCount The count of the number of callsets in this VariantSet
 * @param referenceSetDbId The reference Set DB id this Variant Set is using
 * @param studyDbId The study DB ID used to generate this VariantSet
 * @param variantCount The number of variants in this VariantSet
 * @param variantSetDbId the db id for this variant set.  In the PHG it is analogous to a method(haplotype/path)
 * @param variantSetName The name for this variant set name.
 */
@Serializable
data class VariantSet(
    val additionalInfo: Map<String, String> = mapOf(), val analysis: List<Analysis> = listOf(),
    val availableFormats: List<VariantSetAvailableFormats> = listOf(), val callSetCount: Int,
    val referenceSetDbId: String, val studyDbId: String, val variantCount: Int,
    val variantSetDbId: String, val variantSetName: String
)