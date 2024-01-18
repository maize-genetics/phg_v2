package net.maizegenetics.phgv2.brapi.model

import com.fasterxml.jackson.annotation.JsonValue
import kotlinx.serialization.Serializable

@Serializable
enum class DataFormatEnum(private val value: String) {
    DARTSEQ("DartSeq"), VCF("VCF"), HAPMAP("Hapmap"), TABULAR("tabular"), JSON("JSON");

    override fun toString(): String {
        return value
    }

    @JsonValue
    open fun getFormat(): String? {
        return value
    }

    companion object {
        fun fromValue(text: String): DataFormatEnum? {
            for (b in values()) {
                if (b.value == text) {
                    return b
                }
            }
            return null
        }
    }
}

@Serializable
enum class FileFormatEnum(private val value: String) {
    TEXT_CSV("text/csv"), TEXT_TSV("text/tsv"), APPLICATION_EXCEL("application/excel"),
    APPLICATION_ZIP("application/zip"), APPLICATION_JSON("application/json");

    override fun toString(): String {
        return value
    }

    @JsonValue
    open fun getFileFormat(): String? {
        return value
    }

    companion object {
        fun fromValue(text: String): FileFormatEnum? {
            for (b in values()) {
                if (b.value == text) {
                    return b
                }
            }
            return null
        }
    }
}

/**
 * A VariantSetAvailableFormats is holding all the information regarding the URL, file format and format of the data.
 *
 * @param dataFormat Format of the data
 * @param fileFormat Format of the file(csv, excel, txt)
 * @param fileURL URL of this file.  Can be the endpoint or file url.
 *
*/
@Serializable
data class VariantSetAvailableFormats (val dataFormat : DataFormatEnum, val fileFormat : FileFormatEnum, val fileURL : String)