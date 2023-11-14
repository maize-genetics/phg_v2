package net.maizegenetics.phgv2.cli

import biokotlin.featureTree.Feature
import biokotlin.featureTree.Genome
import biokotlin.featureTree.Strand
import biokotlin.seq.NucSeq
import biokotlin.seqIO.NucSeqIO
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.choice
import com.github.ajalt.clikt.parameters.types.int
import java.io.File

/**
 * Data class to hold the information needed to output a BedRecord.
 */
data class BedRecord(val contig: String, val start : Int, val end: Int, val name: String, val score : Int, val strand: String )

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

    val makeOnlyGenic by option(help = "Set this option to have create-ranges only make genic/cds regions without the in-between regions")
        .flag()

    val referenceFile by option(help = "Full path to the reference fasta file for filling in the intergenic/interCDS regions.  If supplied a Bed region will be created between the last GFF region and the end of the chromosome. If not, this region will be left off..")
        .default("")

    /**
     * Identifies minimum and maximum positions for a list of genes
     *
     * @param genes A list of type `Gene`
     * @param boundary A type of boundary. Either `gene` or `cds`
     * @param pad Number of base-pairs to flank boundary regions
     *
     * @return A list of type `Pair<Int, Int>`
     */
    fun idMinMaxBounds(genes: List<Feature>, boundary: String, pad: Int): List<Pair<Int, Int>> {
        val boundMinMax = when (boundary) {
            "gene" -> genes.map {gene ->
                Pair(gene.start, gene.end)
            }
            "cds" -> {
                genes.map { gene ->
                    gene.children
                        .filter { transcript ->
                            transcript.allAttributes().containsKey("canonical_transcript")
                        }
                        .flatMap { it.children }
                        .filter { it.type == "CDS" }
                        .map { Pair(it.start, it.end) }
                        .let { ranges ->
                            if(ranges.isNotEmpty())
                                Pair(ranges.minOf { it.first }, ranges.maxOf { it.second })
                            else
                                Pair(0,0)
                        }
                }.filter { it.first != 0 && it.second != 0 }
            }
            else -> throw Exception("Undefined boundary")
        }.map { (first, second) ->
           // val modifiedFirst = if (first - pad < 0) 1 else first - pad
           // Pair((modifiedFirst - 1), (second + pad) )//Subtracting one from the start to do the conversion from GFF 1-based inclusive inclusive to BED's 0-based inclusive exclusive

            // LCJ - may not want to do it this way.  May want to fix the overlaps/embedded
            // as we go, which means creating a RangeMap and calling addRange()
            // return the original start and end positions without padding, but subtract 1 from the start position to convert from 1-based to 0-based
            Pair((first - 1), (second) )//Subtracting one from the start to do the conversion from GFF 1-based inclusive inclusive to BED's 0-based inclusive exclusive
        }

        return(boundMinMax)
    }


    /**
     * Function to generate a list of BedRecord objects based on the input GFF file features.
     *
     * @param bounds A list of type `Pair<Int, Int>`. Generated from [idMinMaxBounds]`.
     * @param genes A list of type `Gene`
     * @param featureId Identifier key for value to pull from ID field of GFF. Defaults to `"ID"`
     *
     * @return A list of type `BedRecord`
     *
     */
    fun generateBedRecords(bounds: List<Pair<Int, Int>>,
                        genes: List<Feature>,
                        featureId: String = "ID"
    ) : List<BedRecord> {
        return genes.mapIndexed { index, gene ->
            BedRecord(
                gene.seqid,
                bounds[index].first,
                bounds[index].second,
                gene.attribute(featureId).first(),
                if (gene.score?.isNaN() != false) 0 else gene?.score?.toInt()?:0,
                when(gene.strand) {
                    Strand.PLUS -> "+"
                    Strand.MINUS -> "-"
                    else -> "."
                }
            )
        }
    }

    /**
     * Convert the BedRecord objects into a list of strings for output in BED format
     */
    fun convertBedRecordsIntoOutputStrings(records:List<BedRecord>) : List<String> {
        return records.map { record ->
            "${record.contig}\t${record.start}\t${record.end}\t${record.name}\t${record.score}\t${record.strand}"
        }
    }

    /**
     * Function to fill in regions between the genic bed records with intergenic regions.
     * This will also make a region for the start and end of each chromosome.
     */
    fun fillIntergenicRegions(bedRecords: List<BedRecord>, refSeq: Map<String, NucSeq>) : List<BedRecord> {
        val bedRecordsFilled = mutableListOf<BedRecord>()

        //make an initial BED record for the first region
        val firstRecord = bedRecords[0]
        if(firstRecord.start >= 1) {
            val firstIntergenicRecord = BedRecord(firstRecord.contig, 0, firstRecord.start, "intergenic_${firstRecord.contig}:0-${firstRecord.start}", 0, "+")
            bedRecordsFilled.add(firstIntergenicRecord)
        }
        bedRecordsFilled.add(firstRecord)

        //loop through the rest of the records and fill in the intergenic regions
        for(i in 1 until bedRecords.size) {
            val currentRecord = bedRecords[i]
            val previousRecord = bedRecords[i-1]
            if(currentRecord.contig != previousRecord.contig) {
                //fill in the end of chrom region, and make a new region for the beginning of the next chrom
                val lengthOfLastChrom = refSeq[previousRecord.contig]?.seq()?.length?:0
                if(previousRecord.end < lengthOfLastChrom) {
                    val lastIntergenicRecord = BedRecord(previousRecord.contig, previousRecord.end, lengthOfLastChrom, "intergenic_${previousRecord.contig}:${previousRecord.end}-${lengthOfLastChrom}", 0, "+")
                    bedRecordsFilled.add(lastIntergenicRecord)
                }

                //add in new beginning of next chrom
                if(currentRecord.start >= 1) {
                    val firstIntergenicRecord = BedRecord(currentRecord.contig, 0, currentRecord.start, "intergenic_${currentRecord.contig}:0-${currentRecord.start}", 0, "+")
                    bedRecordsFilled.add(firstIntergenicRecord)
                }
            }
            else {
                val intergenicRecord = BedRecord(
                    currentRecord.contig,
                    previousRecord.end,
                    currentRecord.start,
                    "intergenic_${currentRecord.contig}:${previousRecord.end}-${currentRecord.start}",
                    0,
                    "+"
                )
                bedRecordsFilled.add(intergenicRecord)
            }
            bedRecordsFilled.add(currentRecord)
        }
        //Add a final intergenic if need be
        val lastRecord = bedRecords.last()
        val lengthOfLastChrom = refSeq[lastRecord.contig]?.seq()?.length?:0
        if(lastRecord.end < lengthOfLastChrom) {
            val lastIntergenicRecord = BedRecord(lastRecord.contig, lastRecord.end, lengthOfLastChrom, "intergenic_${lastRecord.contig}:${lastRecord.end}-${lengthOfLastChrom}", 0, "+")
            bedRecordsFilled.add(lastIntergenicRecord)
        }

        return bedRecordsFilled
    }

    override fun run() {
        val genome = Genome.fromFile(gff)

        val genes = genome.iterator().asSequence().filter { it.type == "gene" }.toList()

        // boundMinMax should now be bounds without the padding
        val boundMinMax = idMinMaxBounds(genes, boundary, pad)

        // New method gets called that takes the padding and creates new bounds.
        // It will :
        //  a. toss embedded regions.
        //  b. merge regions that overlap
        //  c. add padding to remaining regions.
        //    1.  if there are not enough bps between regions to allow for total padding,
        //        the bps that do exist will be split between the 2 regions.

        val bedRecords = generateBedRecords(boundMinMax, genes)

        val filledInBedRecords = if(makeOnlyGenic) bedRecords else fillIntergenicRegions(bedRecords, NucSeqIO(referenceFile).readAll())

        val bedLinesToPrint = convertBedRecordsIntoOutputStrings(filledInBedRecords)

        if(output!= null) {
            File(output).bufferedWriter().use { output ->
                bedLinesToPrint.forEach {
                    output.write("$it\n")
                }
            }
        }
    }
}
