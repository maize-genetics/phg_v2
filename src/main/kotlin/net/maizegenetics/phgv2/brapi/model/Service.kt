package net.maizegenetics.phgv2.brapi.model

import com.fasterxml.jackson.annotation.JsonValue
import kotlinx.serialization.Serializable
import java.util.*


@Serializable
enum class WSMIMEDataTypes(private val value: String) {
    APPLICATION_JSON("application/json"), TEXT_CSV("text/csv"), TEXT_TSV("text/tsv"), APPLICATION_FLAPJACK("application/flapjack");

    override fun toString(): String {
        return value
    }

    @JsonValue
    open fun getDataType(): String? {
        return value
    }

    companion object {
        fun fromValue(text: String): WSMIMEDataTypes? {
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
enum class MethodsEnum(private val value: String) {
    GET("GET"), POST("POST"), PUT("PUT"), DELETE("DELETE");

    override fun toString(): String {
        return value
    }

    @JsonValue
    open fun getMethods(): String? {
        return value
    }

    companion object {
        fun fromValue(text: String): MethodsEnum? {
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
enum class VersionsEnum(private val value: String) {
    _0("2.0"), _1("2.1"), _2("2.2");

    override fun toString(): String {
        return value
    }

    @JsonValue
    open fun getVersion(): String? {
        return value
    }

    companion object {
        fun fromValue(text: String): VersionsEnum? {
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
class Service(
    var dataTypes: List<WSMIMEDataTypes>? = null, var methods: List<MethodsEnum> = listOf(),
    var service: String? = null, var versions: List<VersionsEnum> = listOf()
) {

    override fun equals(o: Any?): Boolean {
        if (this === o) {
            return true
        }
        if (o == null || o !is Service) {
            return false
        }
        val service = o as Service
        return dataTypes == service.dataTypes &&
                methods == service.methods &&
                this.service == service.service &&
                versions == service.versions
    }

    override fun hashCode(): Int {
        return Objects.hash(dataTypes, methods, service, versions)
    }

    override fun toString(): String {
        val sb = StringBuilder()
        sb.append("class Service {\n")
        sb.append("    dataTypes: ").append(toIndentedString(dataTypes)).append("\n")
        sb.append("    methods: ").append(toIndentedString(methods)).append("\n")
        sb.append("    service: ").append(toIndentedString(service)).append("\n")
        sb.append("    versions: ").append(toIndentedString(versions)).append("\n")
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