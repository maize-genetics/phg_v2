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

    val dbPath by option(help = "Folder holding TileDB datasets")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--db-path must not be blank"
            }
        }

    val tempDir by option(help = "Folder where temporary files will be written")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--temp-dir must not be blank"
            }
        }

    override fun run() {
//        TODO("Not yet implemented")
    }

}