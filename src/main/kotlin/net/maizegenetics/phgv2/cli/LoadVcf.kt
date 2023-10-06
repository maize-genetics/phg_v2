package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate


class LoadVcf : CliktCommand() {

    val vcfDir by option(help = "VCF file directory")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--vcf-dir must not be blank"
            }
        }

    val dbPath by option(help = "Tile DB URI")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--db-path must not be blank"
            }
        }
    override fun run() {
//        TODO("Not yet implemented")
    }

}