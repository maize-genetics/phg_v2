package net.maizegenetics.phgv2.utils

import java.io.FileInputStream
import java.io.InputStream
import java.security.MessageDigest

/**
 * Allows user to specify the protocol, e.g. MD5, SHA-1, SHA-256
 *
 * @param filename filename to checksum
 * @param protocol protocol
 *
 * @return check sum
 */
fun getChecksum(filename: String, protocol: String = "MD5"): String {
    try {

        val inputStream: InputStream = FileInputStream(filename)
        val digester = MessageDigest.getInstance(protocol)
        val buffer = ByteArray(8192)
        var numOfBytesRead: Int
        while (inputStream.read(buffer).also { numOfBytesRead = it } > 0) {
            digester.update(buffer, 0, numOfBytesRead)
        }
        val hashValue = digester.digest()
        return convertBytesToHex(hashValue)
    } catch (ex: Exception) {
        println(ex.message)
        throw ex
    }

}

private fun convertBytesToHex(bytes: ByteArray): String {
    val builder = StringBuilder()
    for (i in bytes.indices) {
        builder.append(String.format("%02x", bytes[i].toInt() and 0xff))
    }
    return builder.toString()
}