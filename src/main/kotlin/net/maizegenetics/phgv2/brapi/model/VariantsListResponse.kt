package net.maizegenetics.phgv2.brapi.model

import kotlinx.serialization.Serializable

@Serializable
data class VariantsListResponse(
    val metadata: MetadataTokenPagination,
    val result: VariantsListResponseResult,
    var Atcontext: Context? = null
)

@Serializable
data class VariantsListResponseResult(val data: Array<Variant>)