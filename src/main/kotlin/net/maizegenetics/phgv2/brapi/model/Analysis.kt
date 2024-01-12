package model

import kotlinx.serialization.Serializable
import utilities.OffsetDateTimeAsStringSerializer
import java.time.OffsetDateTime

@Serializable
data class Analysis(
    val analysisDbId: String, val analysisName: String,
    @Serializable(with = OffsetDateTimeAsStringSerializer::class)
    val created: OffsetDateTime,
    val description: String, val software: List<String> = listOf(),
    val type: String,
    @Serializable(with = OffsetDateTimeAsStringSerializer::class)
    val updated: OffsetDateTime
)