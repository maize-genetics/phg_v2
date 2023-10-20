package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.OptionGroup
import com.github.ajalt.clikt.parameters.groups.cooccurring
import com.github.ajalt.clikt.parameters.groups.provideDelegate
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.int
import org.apache.logging.log4j.LogManager

/**
 * generic command that pulls data from an AGC compressed file
 * User may specify a list of sample names, list of contigs and
 * a range of positions to pull from the file.
 *
 * If samplename is not specified, all will be pulled
 * If contig is not specified, all will be pulled
 * If range IS specified, contig and sample must be specified.
 *
 * The output for now will be streamed back to user.
 */

class RangeOptions: OptionGroup( help="Options controlling AGC get range start/end") {
    val start by option( )
        .int()
        .default(0).validate {
            require(it >= 0) {
                "--start must be greater than or equal to 0 "
            }
        }
    val end by option (help = "End position of range from which data will be pulled.  If specified, start must also be specified.")
        .int()
        .default(0).validate{
            require((it >= 0) && (it >= start)) {
                "--end must be greater than or equal to --start"
            }
        }
}
class AgcGet : CliktCommand(help="Pull data from an AGC compressed file") {

    private val myLogger = LogManager.getLogger(AgcGet::class.java)

    val dbPath by option(help = "Folder name where AGC compressed filelives.  Should be the same parent folder where tiledb datasets are stored.")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--db-path must not be blank"
            }
        }
    val sampleName by option (help = "Comma separated list of samples from which data will be pulled.  If not specified, all samples will be pulled.")
        .default("")

    val contig by option (help = "Comma separated list of contig names to pull.  If not specified, all contigs will be pulled.")
        .default("")
    val rangeOptions by RangeOptions().cooccurring()


    override fun run() {
        TODO("Not yet implemented")
    }

}

