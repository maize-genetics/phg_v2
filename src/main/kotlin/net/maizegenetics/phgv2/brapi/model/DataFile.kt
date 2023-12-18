package model

import kotlinx.serialization.Serializable

/**
 * A dataFile contains a URL and the relevant file metadata to represent a file
 * @param fileDescription A human readable description of the file contents
 * @param fileMD5Hash The MD5 Hash of the file contents to be used as a check sum
 * @param fileName The name of the file
 * @param fileSize The size of the file in bytes
 * @param fileType The type or format of the file. Preferably MIME Type.
 * @param fileURL The absolute URL where the file is located
 */
@Serializable
data class DataFile (
        /* The absolute URL where the file is located */
        val fileURL: kotlin.String
        ,
        /* A human readable description of the file contents */
        val fileDescription: kotlin.String? = null,
        /* The MD5 Hash of the file contents to be used as a check sum */
        val fileMD5Hash: kotlin.String? = null,
        /* The name of the file */
        val fileName: kotlin.String? = null,
        /* The size of the file in bytes */
        val fileSize: kotlin.Int? = null,
        /* The type or format of the file. Preferably MIME Type. */
        val fileType: kotlin.String? = null
) {
}
