package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.subcommands
import com.github.ajalt.clikt.parameters.options.versionOption
import net.maizegenetics.phgv2.agc.PrepareAssemblies
import net.maizegenetics.phgv2.pathing.BuildKmerIndex
import net.maizegenetics.phgv2.pathing.FindPaths
import net.maizegenetics.phgv2.pathing.MapKmers
import net.maizegenetics.phgv2.simulation.GenerateReads
import net.maizegenetics.phgv2.utils.getBufferedReader
import net.maizegenetics.phgv2.utils.setupDebugLogging

class Phg : CliktCommand() {

    init {
        setupDebugLogging()

        // get version from version.properties file
        var majorVersion = 0
        var minorVersion = 0
        var patchVersion = 0
        var buildNumber = 0

        // Try to get the version.properties file from the jar file
        // If that fails, try to get it from the current working directory
        val reader = try {
            Phg::class.java.getResourceAsStream("/version.properties").bufferedReader()
        } catch (e: Exception) {
            val path = System.getProperty("user.dir")
            println("Getting version from: ${path}/version.properties")
            getBufferedReader("${path}/version.properties")
        }

        reader.readLines().forEach {
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
    .subcommands(
        SetupEnvironment(), Initdb(), CreateRanges(), PrepareAssemblies(), AgcCompress(), AlignAssemblies(),
        CreateAnchorwaveDotplot(), CreateRefVcf(), CreateMafVcf(), Gvcf2Hvcf(), LoadVcf(), ExportVcf(),
        BuildKmerIndex(), MapKmers(), FindPaths(), GenerateReads(), // Imputation
        CreateFastaFromHvcf(), MergeHvcfs(), CalcVcfMetrics(), StartServer // Utilities
    )
    .main(args)
