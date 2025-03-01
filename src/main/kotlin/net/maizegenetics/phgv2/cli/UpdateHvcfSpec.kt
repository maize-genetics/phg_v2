package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.SymbolicAllele
import net.maizegenetics.phgv2.api.exportMultiSampleHVCF
import org.apache.logging.log4j.LogManager

class UpdateHvcfSpec : CliktCommand(help = "Update HVCF file to latest specification.") {

    private val myLogger = LogManager.getLogger(UpdateHvcfSpec::class.java)

    val inputFile by option(help = "Full path to input HVCF file")
        .required()

    val outputFile by option(help = "Full path to updated output HVCF file")
        .required()

    val referenceFile by option(help = "Full path to reference fasta file")

    override fun run() {

        logCommand(this)

        val graph = HaplotypeGraph(listOf(inputFile))

        val idFormat = if (isMd5Checksum(graph.hapIds(graph.ranges().first()).first())) {
            SymbolicAllele.CHECKSUM
        } else {
            SymbolicAllele.RANGE_SAMPLE_GAMETE
        }

        exportMultiSampleHVCF(graph, outputFile, referenceFile, idFormat)

    }

    /**
     * Checks if the input is a valid MD5 checksum.
     */
    private fun isMd5Checksum(input: String): Boolean {
        // Regex to match a valid MD5 checksum (32 hexadecimal characters)
        val md5Regex = Regex("^[a-fA-F0-9]{32}$")
        return md5Regex.matches(input)
    }

}

fun main() {
    val inputFile = "/Users/terry/git/phgv2/data/test/smallseq/LineA.h.vcf"
    val outputFile = "/Users/terry/git/phgv2/data/test/smallseq/LineA.h.updated.vcf"
    val refFile = "/Users/terry/git/phgv2/data/test/smallseq/Ref.fa"
    val result = UpdateHvcfSpec().test(
        "--input-file $inputFile --output-file $outputFile --reference-file $refFile"
    )
}