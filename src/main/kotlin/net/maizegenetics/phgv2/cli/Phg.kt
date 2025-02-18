package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.context
import com.github.ajalt.clikt.core.subcommands
import com.github.ajalt.clikt.output.MordantHelpFormatter
import com.github.ajalt.clikt.parameters.options.versionOption
import net.maizegenetics.phgv2.pathing.*
import net.maizegenetics.phgv2.pathing.ropebwt.ConvertRm2Ps4gFile
import net.maizegenetics.phgv2.pathing.ropebwt.ConvertRopebwt2Ps4gFile
import net.maizegenetics.phgv2.pathing.ropebwt.MapReads
import net.maizegenetics.phgv2.pathing.ropebwt.RopeBwtIndex
import net.maizegenetics.phgv2.utils.phgVersion
import net.maizegenetics.phgv2.utils.setupDebugLogging

/**
 * This class is the main class for the PHG command line interface.
 */
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
        BuildKmerIndex(), MapKmers(), FindPaths(), HapidSampleTable(), SampleHapidByRange(),
        RopeBwtIndex(), MapReads(), // Imputation
        ConvertRm2Ps4gFile(), ConvertRopebwt2Ps4gFile(), // PS4G File creations.
        CreateFastaFromHvcf(), ListSamples(), MergeHvcfs(), MergeGVCFs(), CalcVcfMetrics(), StartServer, ExtractEdgeReads(), //Utilities
        QcReadMapping(), ReadMappingCountQc(), PathsToGff(), // Utilities continued
        CompositeToHaplotypeCoords(), // resequencing pipeline
        InitHvcfArray(), LoadHvcf(), QueryHvcfArrays() // hvcf loading
        ,ConvertRm2Ps4gFile()
    )
    .main(args)
