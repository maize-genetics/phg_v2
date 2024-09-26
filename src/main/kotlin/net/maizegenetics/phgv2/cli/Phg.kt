package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.context
import com.github.ajalt.clikt.core.subcommands
import com.github.ajalt.clikt.output.MordantHelpFormatter
import com.github.ajalt.clikt.parameters.options.versionOption
import net.maizegenetics.phgv2.pathing.BuildKmerIndex
import net.maizegenetics.phgv2.pathing.FindPaths
import net.maizegenetics.phgv2.pathing.MapKmers
import net.maizegenetics.phgv2.utils.phgVersion
import net.maizegenetics.phgv2.utils.setupDebugLogging

class Phg : CliktCommand() {

    init {
        setupDebugLogging()

        context {
            helpFormatter = { MordantHelpFormatter(it, showRequiredTag = true, showDefaultValues = true) }
        }

        val version = phgVersion()
        versionOption(version)
    }

    override fun run() = Unit

}

fun main(args: Array<String>) = Phg()
    .subcommands(
        SetupEnvironment(), Initdb(), CreateRanges(), PrepareAssemblies(), AgcCompress(), AlignAssemblies(), PrepareSlurmAlignFile(),
        CreateAnchorwaveDotplot(), CreateRefVcf(), CreateMafVcf(), Gvcf2Hvcf(), Hvcf2Gvcf(), LoadVcf(), ExportVcf(),
        BuildKmerIndex(), MapKmers(), FindPaths(), // Imputation
        CreateFastaFromHvcf(), ListSamples(), MergeHvcfs(), MergeGVCFs(), CalcVcfMetrics(), StartServer, // Utilities
        CompositeToHaplotypeCoords() // resequencing pipeline
    )
    .main(args)
