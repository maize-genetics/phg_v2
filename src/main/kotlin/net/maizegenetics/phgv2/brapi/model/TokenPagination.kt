package net.maizegenetics.phgv2.brapi.model

import kotlinx.serialization.Serializable

/**
 * @param pageSize The number of data elements returned, aka the size of the currect page. If the requested page does not have enough elements to fill a page at the requested page size, this field should indicate the actual number of elements returned
 * @param currentPageToken The string token used to query the current page of data.
 * @param nextPageToken The string token used to query the next page of data.
 * @param prevPageToken The string token used to query the previous page of data.
 * @param totalCount The total number of elements that are available on the server and match the requested query parameters.
 * @param totalPages The total number of pages of elements available on the server. This should be calculated with the following formula.   totalPages = CEILING( totalCount / requested_page_size)
 */
@Serializable
data class TokenPagination (val pageSize: Int? = null, val currentPageToken: String? = null,
                            val nextPageToken: String? = null, val prevPageToken: String? = null,
                            val totalCount: Int? = null, val totalPages: Int? = null)
