package net.maizegenetics.phgv2.brapi.model

import kotlinx.serialization.Serializable

/**
 *
 * @param Atcontext
 * @param metadata
 * @param result
 */
@Serializable
data class VariantSingleResponse(
    val metadata: MetadataTokenPagination,
    val result: Variant?,
    var Atcontext: Context? = null
)