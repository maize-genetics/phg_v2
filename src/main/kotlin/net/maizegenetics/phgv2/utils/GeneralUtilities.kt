package net.maizegenetics.phgv2.utils

import java.io.BufferedInputStream
import java.lang.Exception
import java.util.logging.Logger

private val myLogger = Logger.getLogger("net.maizegenetics.phgv2.utils.GeneralUtilities")

// This function reads the output from a ProcessBuilder cammond.
// The input is a BufferedInputStream.  This is read one line at a time
// returning the results as a List<String>.   It is used to process output
// from tiledbvcf list --uri <uri> and potentially other commands.
fun inputStreamProcessing(tiledbList: BufferedInputStream): List<String> {
    val samples = mutableListOf<String>()  // this is the list of samples to be returned

    try {
        tiledbList.bufferedReader().use { br ->
            var line = br.readLine()
            while (line != null) {
                // skip blank lines.
                if (line != "") samples.add(line)
                line = br.readLine()
            }
        }
    } catch (exc: Exception) {
        myLogger.severe("Error reading tiledb list output: ${exc.message}")
    } finally {
        tiledbList.close()
    }
    return samples
}
