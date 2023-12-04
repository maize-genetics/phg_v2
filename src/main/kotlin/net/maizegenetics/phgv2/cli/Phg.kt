package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.subcommands
import com.github.ajalt.clikt.parameters.options.versionOption
import net.maizegenetics.phgv2.agc.AnnotateFasta
import net.maizegenetics.phgv2.utils.setupDebugLogging

class Phg : CliktCommand() {

    // Need an automated way to get version from build.gradle.kts
    private val version = "2.2.0"

    init {
        setupDebugLogging()
        versionOption(version)
    }

    override fun run() = Unit

}

fun main(args: Array<String>) = Phg()
    .subcommands(SetupEnvironment(), Initdb(),  CreateRanges(), AnnotateFasta(), AgcCompress(), AlignAssemblies(), CreateRefVcf(), CreateMafVcf(), LoadVcf(), ExportVcf(), CreateFastaFromHvcf())
    .main(args)
