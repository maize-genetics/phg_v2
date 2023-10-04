package net.maizegenetics.phgv2.cli

import biokotlin.featureTree.Gene
import biokotlin.featureTree.Genome
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.choice
import com.github.ajalt.clikt.parameters.types.int
import java.io.File

/**
 * A [CliktCommand] class for generating reference ranges
 *
 * Reference ranges are parsed from an input GFF file and converted to a BED file
 */
class CreateRanges: CliktCommand(help="Create BED file of reference ranges from GFF file") {
    val gff by option(help = "GFF file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--gff must not be blank"
            }
        }
    val boundary by option(help = "Reference range boundaries")
        .choice("gene", "cds")
        .default("gene")
    val pad by option(help = "Number of base pairs to flank regions")
        .int()
        .default(0)
    val output by option("-o", "--output", help = "Name for output BED file")

    /**
     * Identifies minimum and maximum positions for a list of genes
     *
     * @param genes A list of type `Gene`
     * @param boundary A type of boundary. Either `gene` or `cds`
     * @param pad Number of base-pairs to flank boundary regions
     *
     * @return A list of type `Pair<Int, Int>`
     */
    fun idMinMaxBounds(genes: List<Gene>, boundary: String, pad: Int): List<Pair<Int, Int>> {
        val boundMinMax = when (boundary) {
            "gene" -> genes.map {gene ->
                Pair(gene.start, gene.end)
            }
            "cds" -> {
                genes.map { gene ->
                    gene.children
                        .filter { transcript ->
                            transcript.attributes.containsKey("canonical_transcript")
                        }
                        .flatMap { it.children }
                        .filter { it.type().name == "CODING_SEQUENCE" }
                        .map { Pair(it.start, it.end) }
                        .let { ranges ->
                            Pair(ranges.minOf { it.first }, ranges.maxOf { it.second })
                        }
                }
            }
            else -> throw Exception("Undefined boundary")
        }.map { (first, second) ->
            val modifiedFirst = if (first - pad < 0) 1 else first - pad
            Pair((modifiedFirst - 1), (second + pad) - 1)
        }

        return(boundMinMax)
    }

    /**
     * Generates a list of BED-formatted string rows
     *
     * @param bounds A list of type `Pair<Int, Int>`. Generated from [idMinMaxBounds]`.
     * @param genes A list of type `Gene`
     * @param delimiter A column delimiter for a BED file. Defaults to `\t`
     * @param featureId Identifier key for value to pull from ID field of GFF. Defaults to `"ID"`
     *
     * @return A list of type `String`
     */
    fun generateBedRows(
        bounds: List<Pair<Int, Int>>,
        genes: List<Gene>,
        delimiter: String = "\t",
        featureId: String = "ID"
    ): List<String> {
        val bedLinesToPrint = genes.mapIndexed { index, gene ->
            listOf(
                gene.seqid,
                bounds[index].first,
                bounds[index].second,
                gene.attributes[featureId],
                if (gene.score.isNaN()) 0 else gene.score,
                gene.strand
            ).joinToString(delimiter)
        }

        return(bedLinesToPrint)
    }

    override fun run() {
        val genome = Genome.fromGFF(gff)
        val genes = genome.genes()

        val boundMinMax = idMinMaxBounds(genes, boundary, pad)

        val bedLinesToPrint = generateBedRows(boundMinMax, genes)

        if(output!= null) {
            File(output).bufferedWriter().use { output ->
                bedLinesToPrint.forEach {
                    output.write("$it\n")
                }
            }
        }
    }
}
