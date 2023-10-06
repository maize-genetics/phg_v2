package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate

class FastaGenerator : CliktCommand() {

    val dbPath by option(help = "Tile DB URI")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--db-path must not be blank"
            }
        }

    val agcFile by option(help = "AGC index file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--agc-file must not be blank"
            }
        }

    val sampleName by option(help = "Sample name")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--sample-name must not be blank"
            }
        }

    val output by option("-o", "--output", help = "Name for output Fasta file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--output/-o must not be blank"
            }
        }

    override fun run() {
//        TODO("Not yet implemented")
    }
}