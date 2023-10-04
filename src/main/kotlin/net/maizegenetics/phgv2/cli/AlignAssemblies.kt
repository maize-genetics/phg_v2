package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int

class AlignAssemblies : CliktCommand() {

    val gff by option(help = "Full path to the reference gff file")
        .required()

    val ref by option(help = "Full path to reference fasta file")
        .required()

    val assemblies by option(
        "-a",
        "--assemblies",
        help = "File containing list of assemblies to align, 1 per line, full path to file"
    )
        .required()

    val outputDir by option("-o", "--outputDir", help = "Directory where temporary and final files will be written")
        .required()

    val threads by option(help = "Number of threads to use for each assembly processed")
        .int()
        .default(1)

    val runs by option(
        help = "Number of assemblies to simultaneously process. " +
                "If you have 4 threads and 2 runs, then 2 assemblies will be processed at a time, each using 4 threads." +
                "The anchorwave application can take up to 50G per thread for each assembly processed, plus some overhead." +
                "Consider this memory factor when providing values for threadsPerRun and numRuns."
    )
        .int()
        .default(1)

    val refMaxAlignCov by option(help = "Anchorwave proali parameter R, indicating reference genome maximum alignment coverage.")
        .int()
        .default(1)

    val queryMaxAlignCov by option(help = "Anchorwave proali parameter Q, indicating query genome maximum alignment coverage.")
        .int()
        .default(1)

    override fun run() {
        // TBD
    }

}