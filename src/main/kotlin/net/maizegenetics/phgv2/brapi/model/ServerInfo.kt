package net.maizegenetics.phgv2.brapi.model

import kotlinx.serialization.Serializable
import java.util.*

@Serializable
class ServerInfo(
    var calls: MutableList<Service> = mutableListOf<Service>(), var contactEmail: String? = null,
    val documentationURL: String? = null, val location: String? = null, val organizationName: String? = null,
    val organizationURL: String? = null, val serverDescription: String? = null, val serverName: String? = null
) {

    fun addCallsItem(callsItem: Service) {
        calls.add(callsItem)
    }

    override fun equals(o: Any?): Boolean {
        if (this === o) {
            return true
        }
        if (o == null || o !is ServerInfo) {
            return false
        }
        val serverInfo: ServerInfo = o as ServerInfo
        return calls == serverInfo.calls &&
                contactEmail == serverInfo.contactEmail &&
                documentationURL == serverInfo.documentationURL &&
                location == serverInfo.location &&
                organizationName == serverInfo.organizationName &&
                organizationURL == serverInfo.organizationURL &&
                serverDescription == serverInfo.serverDescription &&
                serverName == serverInfo.serverName
    }

    override fun hashCode(): Int {
        return Objects.hash(
            calls,
            contactEmail,
            documentationURL,
            location,
            organizationName,
            organizationURL,
            serverDescription,
            serverName
        )
    }

    override fun toString(): String {
        val sb = StringBuilder()
        sb.append("class ServerInfo {\n")
        sb.append("    calls: ").append(toIndentedString(calls)).append("\n")
        sb.append("    contactEmail: ").append(toIndentedString(contactEmail)).append("\n")
        sb.append("    documentationURL: ").append(toIndentedString(documentationURL)).append("\n")
        sb.append("    location: ").append(toIndentedString(location)).append("\n")
        sb.append("    organizationName: ").append(toIndentedString(organizationName)).append("\n")
        sb.append("    organizationURL: ").append(toIndentedString(organizationURL)).append("\n")
        sb.append("    serverDescription: ").append(toIndentedString(serverDescription)).append("\n")
        sb.append("    serverName: ").append(toIndentedString(serverName)).append("\n")
        sb.append("}")
        return sb.toString()
    }

    private fun toIndentedString(o: Any?): String? {
        return o?.toString()?.replace("\n", "\n    ") ?: "null"
    }

}