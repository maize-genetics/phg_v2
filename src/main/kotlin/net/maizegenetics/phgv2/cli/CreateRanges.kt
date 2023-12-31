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
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import com.google.common.collect.TreeRangeMap
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.addRange
import net.maizegenetics.phgv2.utils.createFlankingList
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
        .validate {
            require(it.isNotBlank()) {
                "--reference-file must not be blank"
            }
        }


    /**
     * Identifies minimum and maximum positions for a list of genes
     *
     * @param genes A list of type `Gene`
     * @param boundary A type of boundary. Either `gene` or `cds`
     * @param pad Number of base-pairs to flank boundary regions
     *
     * Initially creates a `Triple<Position, Position, String>`
     * The "String" in the Triple is a concatenated list of Gene ID.  Will be 1 or more,
     * based on whethere overlapping genes were merged.
     *
     * This Triple is further processed to create a RangeMap<Position, String>
     * of merged genes:  embedded genes are tossed, overlapping genes are merged.
     */
    fun idMinMaxBounds(refSeq: Map<String, NucSeq>, genes: List<Feature>, boundary: String): RangeMap<Position, String> {
        val boundMinMax = when (boundary) {
            "gene" -> genes.map {gene ->
                Triple(Position(gene.seqid,gene.start), Position(gene.seqid,gene.end), gene.attribute("ID").first())
            }
            "cds" -> {
                genes.map { gene ->
                    val canonical = gene.children.filter { transcript -> transcript.allAttributes().containsKey("canonical_transcript") }
                    // if there is no canonical transcript, use the first transcript
                    val transcripts = canonical.ifEmpty { listOf<Feature>(gene.children.first()) }
                    transcripts
                        .flatMap { it.children }
                        .filter { it.type == "CDS" }
                        .map {
                            Triple(Position(gene.seqid,it.start), Position(gene.seqid,it.end), gene.attribute("ID").first()) }// ok to here on changes
                        .let { ranges ->
                            // Ranges is List<Triple<Position, Position, String>>
                            // We want the min and max from all the cds's that are marked canonical
                            // for this gene.
                            if(ranges.isNotEmpty()){
                                val id = gene.attribute("ID").first()
                                Triple(Position(gene.seqid,ranges.minOf { it.first.position }), Position(gene.seqid,ranges.maxOf { it.second.position }),id)
                            }
                            else
                                Triple(Position("0",0),Position("0",0),"")
                        }
                }.filter { it.first.position != 0 && it.second.position != 0 }
            }
            else -> throw Exception("Undefined boundary")
        }.map { (first, second,third) ->
            // GFF is 1-based.  The BED file is 0-based, so we need to subtract 1 from the start position
            Triple(Position(first.contig,first.position-1), Position(second.contig,second.position),third)
        }

        // Resolve overlaps/embedded regions
        val geneRange: RangeMap<Position, String> = TreeRangeMap.create()
        for (bound in boundMinMax) {
            val boundRange = Range.closed(bound.first, bound.second)
            // addRange will merge overlaps, toss embedded regions,
            // and concat gene names when ranges are merged.
            addRange(geneRange,boundRange,bound.third)
        }
        return geneRange
    }


    /**
     * Function to generate a list of BedRecords based on a RangeMap of Position.
     * Because overlapping genes were merged, the score and strand values of the
     * BedRecord can not be determined and are set to 0 and "." respectively.
     * Also, because of this, the number of bed records does not necessarily match
     *  the number of genes in the GFF file.  THere may be fewer if overlaps were merged.
     */
    fun generateBedRecords(ranges: RangeMap< Position,String>): List<BedRecord> {

        val bedRecordList = mutableListOf<BedRecord>()
        // Are these sorted?
        ranges.asMapOfRanges().forEach{ range ->
            bedRecordList.add(BedRecord(range.key.lowerEndpoint().contig, range.key.lowerEndpoint().position, range.key.upperEndpoint().position, range.value, 0, "."))
        }
        return bedRecordList
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

        val refSeq = NucSeqIO(referenceFile).readAll()
        // Determine the boundaries with flanking for the gene or CDS regions
        val geneRange = idMinMaxBounds(refSeq, genes, boundary)

        // Add flanking regions if specified by user
        // var flankingRange:RangeMap<Position,String> = TreeRangeMap.create()
        val flankingRange = if (pad > 0) createFlankingList(geneRange, pad, refSeq) else geneRange

        // Create bed records from the ranges
        val bedRecords = generateBedRecords(flankingRange)

        val filledInBedRecords = if(makeOnlyGenic) bedRecords else fillIntergenicRegions(bedRecords, refSeq)

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
