package model

import com.fasterxml.jackson.annotation.JsonProperty
import kotlinx.serialization.Serializable
import java.util.*

/**
 *
 * @param currentPage The index number for the returned page of data. This should always match the requested page number or the default page (0).
 */
@Serializable
 class IndexPagination: BasePagination() {
        /* The index number for the returned page of data. This should always match the requested page number or the default page (0). */
        @JsonProperty("currentPage")
        var currentPage: kotlin.Int = 0


    override fun equals(o: Any?): Boolean {
        if (this === o) {
            return true
        }
        if (o == null || javaClass != o.javaClass) {
            return false
        }
        val indexPagination = o as IndexPagination
        return this.currentPage == indexPagination.currentPage &&
                super.equals(o)
    }

    override fun hashCode(): Int {
        return Objects.hash(currentPage, super.hashCode())
    }

    override fun toString(): String {
        val sb = StringBuilder()
        sb.append("class IndexPagination {\n")
        sb.append("    ").append(toIndentedString(super.toString())).append("\n")
        sb.append("    currentPage: ").append(toIndentedString(currentPage)).append("\n")
        sb.append("}")
        return sb.toString()
    }

//    /**
//     * Convert the given object to string with each line indented by 4 spaces
//     * (except the first line).
//     * LCJ - cannot get this one to work - says it cannot override BasePagination version as that is final
//     */
//     fun toIndentedString(o: Any?): String? {
//        return o?.toString()?.replace("\n", "\n    ") ?: "null"
//    }
}
