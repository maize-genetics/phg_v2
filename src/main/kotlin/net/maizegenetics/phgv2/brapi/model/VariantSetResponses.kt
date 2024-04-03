package net.maizegenetics.phgv2.brapi.model

import kotlinx.serialization.Serializable

@Serializable
data class VariantSetResponse(
    val metadata: Metadata,
    val result: VariantSet,
    var Atcontext: Context? = null
)

@Serializable
data class VariantSetsListResponse(
    val metadata: Metadata,
    val result: VariantSetsListResponseResult,
    var Atcontext: Context? = null
)

@Serializable
data class VariantSetsListResponseResult(val data: Array<VariantSet>)
