package net.maizegenetics.phgv2.brapi.model

import kotlinx.serialization.Serializable
import java.util.*


/**
 * The JSON-LD Context is used to provide JSON-LD definitions to each field in a JSON object. By providing an array of context file urls, a BrAPI response object becomes JSON-LD compatible.    For more information, see https://w3c.github.io/json-ld-syntax/#the-context
 */
@Serializable
class Context : ArrayList<String?>() {
    override fun equals(o: Any?): Boolean {
        if (this === o) {
            return true
        }
        return !(o == null || javaClass != o.javaClass)
    }

    override fun hashCode(): Int {
        return Objects.hash(super.hashCode())
    }

    override fun toString(): String {
        val sb = StringBuilder()
        sb.append("class Context {\n")
        sb.append("    ").append(toIndentedString(super.toString())).append("\n")
        sb.append("}")
        return sb.toString()
    }

    /**
     * Convert the given object to string with each line indented by 4 spaces
     * (except the first line).
     */
    private fun toIndentedString(o: Any?): String {
        return o?.toString()?.replace("\n", "\n    ") ?: "null"
    }
}