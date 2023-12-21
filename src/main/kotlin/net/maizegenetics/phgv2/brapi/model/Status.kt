package model

import com.fasterxml.jackson.annotation.JsonCreator
import com.fasterxml.jackson.annotation.JsonProperty
import com.fasterxml.jackson.annotation.JsonValue
import kotlinx.serialization.Serializable
import java.util.*

/**
 * An array of status messages to convey technical logging information from the server to the client.
 * @param message A short message concerning the status of this request/response
 * @param messageType The logging level for the attached message
 */
@Serializable
data class Status (
        /* A short message concerning the status of this request/response */
        val message: String,
        /* The logging level for the attached message */
        val messageType: Status.MessageType

) {
    /**
     * The logging level for the attached message
     * Values: dEBUG,eRROR,wARNING,iNFO
     */
    enum class MessageType(val value: String){
        dEBUG("DEBUG"),
        eRROR("ERROR"),
        wARNING("WARNING"),
        iNFO("INFO");
    }
}
