package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand

/**
 * generic command that pulls data from an AGC compressed file
 * User may specify a list of sample names, list of contigs and
 * a range of positions to pull from the file.
 *
 *
 */
class AgcGet : CliktCommand(help="Pull data from an AGC compressed file") {
    override fun run() {
        TODO("Not yet implemented")
    }
}