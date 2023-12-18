package net.maizegenetics.phgv2.brapi.model

import kotlinx.serialization.Serializable

/**
 * @param additionalInfo Additional arbitrary info
 * @param sampleName The name of the sample
 * @param sampleDescription Description of a sample MIAPPE V1.1 (DM-79) Sample description - Any information not captured by the other sample fields, including quantification, sample treatments and processing.
 * @param sampleDbId The ID which uniquely identifies a sample MIAPPE V1.1 (DM-76) Sample ID - Unique identifier for the sample.
 */
@Serializable
data class Sample(
    /* The name of the sample */
    var sampleName: String,
    var sampleDbId: String,
    /* Description of a sample  MIAPPE V1.1 (DM-79) Sample description - Any information not captured by the other sample fields, including quantification, sample treatments and processing. */
    var sampleDescription: String? = null,
    /* Additional arbitrary info */
    var additionalInfo: Map<String, String>? = null
)