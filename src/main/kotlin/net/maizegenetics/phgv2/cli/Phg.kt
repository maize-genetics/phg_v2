package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.subcommands
import com.github.ajalt.clikt.parameters.options.versionOption
import net.maizegenetics.phgv2.agc.AnnotateFastas
import net.maizegenetics.phgv2.utils.setupDebugLogging

class Phg : CliktCommand() {

    init {
        setupDebugLogging()

        // get version from version.properties file
        var majorVersion = 0
        var minorVersion = 0
        var patchVersion = 0
        var buildNumber = 0
        Phg::class.java.getResourceAsStream("/version.properties").bufferedReader().readLines().forEach {
            val (key, value) = it.split("=")
            when (key) {
                "majorVersion" -> majorVersion = value.toInt()
                "minorVersion" -> minorVersion = value.toInt()
                "patchVersion" -> patchVersion = value.toInt()
                "buildNumber" -> buildNumber = value.toInt()
            }
        }
        val version = "$majorVersion.$minorVersion.$patchVersion.$buildNumber"
        versionOption(version)
    }

    override fun run() = Unit

}

fun main(args: Array<String>) = Phg()
    .subcommands(SetupEnvironment(), Initdb(),  CreateRanges(), AnnotateFastas(), AgcCompress(), AlignAssemblies(), CreateRefVcf(), CreateMafVcf(), LoadVcf(), ExportVcf(), CreateFastaFromHvcf())
    .main(args)
