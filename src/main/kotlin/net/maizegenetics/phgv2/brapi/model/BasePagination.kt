package net.maizegenetics.phgv2.brapi.model

import kotlinx.serialization.Serializable
import java.util.*

/**
 *
 * @param pageSize The number of data elements returned, aka the size of the current page. If the requested page does not have enough elements to fill a page at the requested page size, this field should indicate the actual number of elements returned.
 * @param totalCount The total number of elements that are available on the server and match the requested query parameters.
 * @param totalPages The total number of pages of elements available on the server. This should be calculated with the following formula.   totalPages = CEILING( totalCount / requested_page_size)
 */
@Serializable
open class BasePagination(
        /* The number of data elements returned, aka the size of the current page. If the requested page does not have enough elements to fill a page at the requested page size, this field should indicate the actual number of elements returned. */
        var pageSize: Int = 20 // LCJ - I arbitrarily chose 20
        ,
        /* The total number of elements that are available on the server and match the requested query parameters. */
        var totalCount: Int? = null,
        /* The total number of pages of elements available on the server. This should be calculated with the following formula.   totalPages = CEILING( totalCount / requested_page_size) */
        var totalPages: Int? = null
) {
    override fun equals(o: Any?): Boolean {
        if (this === o) {
            return true
        }
        if (o == null || javaClass != o.javaClass) {
            return false
        }
        val basePagination = o as BasePagination
        return this.pageSize == basePagination.pageSize &&
                this.totalCount == basePagination.totalCount &&
                this.totalPages == basePagination.totalPages
    }

    override fun hashCode(): Int {
        return Objects.hash(pageSize, totalCount, totalPages)
    }

    override fun toString(): String {
        val sb = StringBuilder()
        sb.append("class BasePagination {\n")
        sb.append("    pageSize: ").append(toIndentedString(pageSize)).append("\n")
        sb.append("    totalCount: ").append(toIndentedString(totalCount)).append("\n")
        sb.append("    totalPages: ").append(toIndentedString(totalPages)).append("\n")
        sb.append("}")
        return sb.toString()
    }

    /**
     * Convert the given object to string with each line indented by 4 spaces
     * (except the first line).
     */
     fun toIndentedString(o: Any?): String? {
        return o?.toString()?.replace("\n", "\n    ") ?: "null"
    }
}
