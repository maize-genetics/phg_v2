package net.maizegenetics.phgv2.brapi.model

import kotlinx.serialization.Serializable
import utilities.OffsetDateTimeAsStringSerializer
import java.time.OffsetDateTime

/**
 * A `Variant` represents the a location where there is a variant across different Callsets in the variantset
 * @param additionalInfo Additional arbitrary info
 * @param alternateBases The list of alternative bases for a given variant.  These can be alleles or haplotypes
 * @param ciend the ending of a structural variant.  Currently is not used.
 * @param cipos the start of a structural variant.  Currently is not used.
 * @param created the date of when this record was created
 * @param end the end position of this variant record.
 * @param filtersApplied A boolean for  if filters have been applied or not to this variant
 * @param filtersFailed A list of all the filters which have been failed for this variant
 * @param filtersPassed A boolean for if any filters have been passed.
 * @param referenceBases A String of the reference bases for this Variant.
 * @param referenceName The name of the reference sequence.  This should match what is in references.
 * @param start The start position of this Variant.  In the case of a single position variant, it will just be where the variant is.  1-based inclusive like VCF.
 * @param svlen The structural variantion length of this variant.  Current is not used.
 * @param updated The Date of when this record was previously updated.
 * @param variantDbId The db id representing this variant
 * @param variantNames The list of names for this variant
 * @param variantSetDbId The list of variantset DB Ids which hold this variant.
 * @param variantType The type of variant this is.  Can be SNP or ReferenceRange.
 */
@Serializable
data class Variant(
    val additionalInfo: Map<String, String> = mapOf(),
    val alternateBases: List<String> = listOf(),
    val ciend: List<Int> = listOf(),
    val cipos: List<Int> = listOf(),
    @Serializable(with = OffsetDateTimeAsStringSerializer::class)
    val created: OffsetDateTime? = null,
    val end: Int, val filtersApplied: Boolean = false,
    val filtersFailed: List<String> = listOf(),
    val filtersPassed: Boolean = true,
    val referenceBases: String,
    val referenceName: String,
    val start: Int,
    val svlen: Int,
    @Serializable(with = OffsetDateTimeAsStringSerializer::class)
    val updated: OffsetDateTime? = null,
    val variantDbId: String,
    val variantNames: List<String> = listOf(),
    val variantSetDbId: List<String> = listOf(),
    val variantType: String
)
