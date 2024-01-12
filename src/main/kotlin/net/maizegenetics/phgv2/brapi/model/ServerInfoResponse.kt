package net.maizegenetics.phgv2.brapi.model

import kotlinx.serialization.Serializable
import model.Context
import model.Metadata
import java.util.*

@Serializable
data class ServerInfoResponse(
    val metadata: Metadata,
    val result: ServerInfo,
    val _atContext: Context? = null
) {

    override fun equals(serverInfoResponse: Any?): Boolean {
        if (this === serverInfoResponse) {
            return true
        }
        if (serverInfoResponse == null || serverInfoResponse !is ServerInfoResponse) {
            return false
        }

        return this._atContext == serverInfoResponse._atContext &&
                metadata == serverInfoResponse.metadata &&
                result == serverInfoResponse.result
    }

    override fun hashCode(): Int {
        return Objects.hash(_atContext, metadata, result)
    }

    override fun toString(): String {
        val sb = StringBuilder()
        sb.append("class ServerInfoResponse {\n")
        sb.append("    _atContext: ").append(toIndentedString(_atContext)).append("\n")
        sb.append("    metadata: ").append(toIndentedString(metadata)).append("\n")
        sb.append("    result: ").append(toIndentedString(result)).append("\n")
        sb.append("}")
        return sb.toString()
    }

    /**
     * Convert the given object to string with each line indented by 4 spaces
     * (except the first line).
     */
    private fun toIndentedString(o: Any?): String? {
        return o?.toString()?.replace("\n", "\n    ") ?: "null"
    }

}