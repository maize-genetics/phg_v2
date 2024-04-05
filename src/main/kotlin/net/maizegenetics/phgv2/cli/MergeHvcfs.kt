package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.enum
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.SymbolicAllele
import net.maizegenetics.phgv2.api.exportMultiSampleHVCF
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * Merge multiple HVCF files into a single HVCF file.
 * Uses HaplotypeGraph to merge the HVCF files.
 * And exports the merged HVCF file using ExportHaplotypeGraph.exportMultiSampleHVCF().
 */
class MergeHvcfs : CliktCommand(help = "Merge multiple HVCF files into a single HVCF file.") {

    private val myLogger = LogManager.getLogger(MergeHvcfs::class.java)

    val inputDir by option(help = "Full path to input HVCF file directory")
        .required()

    val idFormat by option(help = "ID format for the HVCF files")
        .enum<SymbolicAllele>()
        .default(SymbolicAllele.RANGE_SAMPLE_GAMETE)

    val referenceFile by option(help = "Full path to reference fasta file")

    val outputFile by option(help = "Full path to output HVCF file")
        .required()

    val rangeBedfile by option(help = "Full path to range bedfile")

    override fun run() {

        val inputFiles = File(inputDir)
            .walk()
            .filter {
                it.isFile && (it.name.endsWith(".h.vcf") || it.name.endsWith(".h.vcf.gz") ||
                        it.name.endsWith(".hvcf") || it.name.endsWith(".hvcf.gz"))
            }
            .map { it.absolutePath }
            .toList()

        val graph = HaplotypeGraph(inputFiles)

        exportMultiSampleHVCF(graph, outputFile, referenceFile, idFormat, rangeBedfile)

    }

}