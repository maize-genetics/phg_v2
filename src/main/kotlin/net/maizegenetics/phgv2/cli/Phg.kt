package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.context
import com.github.ajalt.clikt.core.subcommands
import com.github.ajalt.clikt.output.MordantHelpFormatter
import com.github.ajalt.clikt.parameters.options.versionOption
import net.maizegenetics.phgv2.pathing.*
import net.maizegenetics.phgv2.pathing.ropebwt.*
import net.maizegenetics.phgv2.scaffolding.ComputeScaffoldCorrelation
import net.maizegenetics.phgv2.scaffolding.ProcessMemFile
import net.maizegenetics.phgv2.scaffolding.ProcessPafFile
import net.maizegenetics.phgv2.scaffolding.ProcessSamFile
import net.maizegenetics.phgv2.scaffolding.ReportMemAlignments
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
        BuildKmerIndex(), MapKmers(), FindPaths(), HapidSampleTable(), SampleHapidByRange(),
        RopeBwtIndex(), RopeBwtChrIndex(),MapReads(), ImputationMetrics(), // Imputation
        BuildSplineKnots(),ConvertRm2Ps4gFile(), ConvertRopebwt2Ps4gFile(), ConvertVcf2Ps4gFile(), // PS4G File creations.
        CreateFastaFromHvcf(), ListSamples(), MergeHvcfs(), MergeGVCFs(), CalcVcfMetrics(), StartServer, ExtractEdgeReads(), //Utilities
        QcReadMapping(), ReadMappingCountQc(), PathsToGff(), // Utilities continued
        CompositeToHaplotypeCoords(), // resequencing pipeline
        InitHvcfArray(), LoadHvcf(), QueryHvcfArrays(), // hvcf loading
        UpdateHvcfSpec(), // hvcf updating
        ProcessPafFile(), ProcessSamFile(), ProcessMemFile(), ReportMemAlignments(),ComputeScaffoldCorrelation() //Scaffolding  TODO delete and move to better place

    )
    .main(args)
