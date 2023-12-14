package net.maizegenetics.phgv2.utils

import java.io.*
import java.net.URL
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

fun getBufferedReader(forFile: File, bufSize: Int = 8192): BufferedReader {
    return if (bufSize < 1) {
        getBufferedReader(forFile)
    } else if (forFile.name.startsWith("http")) {
        if (forFile.name.endsWith(".gz")) {
            BufferedReader(InputStreamReader(GZIPInputStream(URL(forFile.name).openStream(), bufSize)), bufSize)
        } else {
            BufferedReader(InputStreamReader(URL(forFile.name).openStream()), bufSize)
        }
    } else if (forFile.name.endsWith(".gz")) {
        BufferedReader(InputStreamReader(GZIPInputStream(FileInputStream(forFile), bufSize)), bufSize)
    } else {
        BufferedReader(InputStreamReader(FileInputStream(forFile)), bufSize)
    }
}

fun getBufferedReader(filename: String, bufSize: Int = 8192): BufferedReader {
    return getBufferedReader(File(filename), bufSize)
}

fun getBufferedWriter(file: File, append: Boolean): BufferedWriter {
    return if (file.name.endsWith(".gz")) {
        BufferedWriter(OutputStreamWriter(GZIPOutputStream(FileOutputStream(file, append))))
    } else {
        BufferedWriter(OutputStreamWriter(FileOutputStream(file, append)))
    }
}

fun getBufferedWriter(filename: String, append: Boolean): BufferedWriter {
    return getBufferedWriter(File(filename), append)
}

