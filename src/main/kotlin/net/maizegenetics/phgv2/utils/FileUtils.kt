package net.maizegenetics.phgv2.utils

import java.io.*
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

/**
 * Utilities for working with files and, potentially, other input and output streams.
 * Currently, the only methods implemented are for getting buffered readers and writers
 * for text files and gzipped text files.
 *
 * The code was ported from TASSEL. That code checked for http sources. That
 * feature was not ported because of a lack of a good source for testing and could be
 * added if needed.
 */

/**
 * Gets a BufferedReader for a text file. Checks for a .gz ending and
 * handles that as a gzip file. [bufSize] can be set but defaults to 8192.
 */
fun getBufferedReader(forFile: File, bufSize: Int = 8192): BufferedReader {
    return if (bufSize < 1) {
        getBufferedReader(forFile)
    } else if (forFile.name.endsWith(".gz")) {
        BufferedReader(InputStreamReader(GZIPInputStream(FileInputStream(forFile), bufSize)), bufSize)
    } else {
        BufferedReader(InputStreamReader(FileInputStream(forFile)), bufSize)
    }
}

/**
 * Gets a BufferedReader for a text file. Checks for a .gz ending and
 * handles that as a gzip file. [bufSize] can be set but defaults to 8192.
 */
fun getBufferedReader(filename: String, bufSize: Int = 8192): BufferedReader {
    return getBufferedReader(File(filename), bufSize)
}

/**
 * Gets a BufferedWriter for a text file. If the file name ends in .gz it
 * will create a gzipped file. If [append] is true and the file exists, the
 * content will be appended. Otherwise, the original file will be overwritten.
 */
fun getBufferedWriter(file: File, append: Boolean): BufferedWriter {
    return if (file.name.endsWith(".gz")) {
        BufferedWriter(OutputStreamWriter(GZIPOutputStream(FileOutputStream(file, append))))
    } else {
        BufferedWriter(OutputStreamWriter(FileOutputStream(file, append)))
    }
}

/**
 * Gets a BufferedWriter for a text file. If the file name ends in .gz it
 * will create a gzipped file. If [append] is true and the file exists, the
 * content will be appended. Otherwise, the original file will be overwritten.
 */
fun getBufferedWriter(filename: String, append: Boolean): BufferedWriter {
    return getBufferedWriter(File(filename), append)
}

