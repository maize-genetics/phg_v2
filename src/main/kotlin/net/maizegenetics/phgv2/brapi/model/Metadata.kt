package net.maizegenetics.phgv2.brapi.model

import com.fasterxml.jackson.annotation.JsonProperty
import kotlinx.serialization.Serializable
import java.util.*

/**
 * Metadata
 */
@Serializable
class Metadata : MetadataBase() {
    @JsonProperty("pagination")
    var pagination: IndexPagination? = null

    override fun equals(o: Any?): Boolean {
        if (this === o) {
            return true
        }
        if (o == null || javaClass != o.javaClass) {
            return false
        }
        val metadata = o as Metadata
        return pagination == metadata.pagination &&
                super.equals(o)
    }

    override fun hashCode(): Int {
        return Objects.hash(pagination, super.hashCode())
    }

    override fun toString(): String {
        val sb = StringBuilder()
        sb.append("class Metadata {\n")
        sb.append("    ").append(toIndentedString(super.toString())).append("\n")
        sb.append("    pagination: ").append(toIndentedString(pagination)).append("\n")
        sb.append("}")
        return sb.toString()
    }

    /**
     * Convert the given object to string with each line indented by 4 spaces
     * (except the first line).
     */
     fun toIndentedString(o: Any?): String {
        return o?.toString()?.replace("\n", "\n    ") ?: "null"
    }
}
