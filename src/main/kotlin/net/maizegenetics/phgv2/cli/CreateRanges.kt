package net.maizegenetics.phgv2.cli

import biokotlin.featureTree.Feature
import biokotlin.featureTree.Genome
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
 * This will be used as records for a BED file, so the start and end are 0-based, with start-inclusive, end-exclusive.
 */
data class BedRecord(val contig: String, val start : Int, val end: Int, val name: String, val score : Int, val strand: String )

/**
 * A [CliktCommand] class for generating reference ranges
 *
 * Reference ranges are parsed from an input GFF file and converted to a BED file.
 * Required parameters are the gff file and the reference file.  The reference file is used to fill in intergenic regions.
 * Overlapping genes/CDS regions are merged, embedded regions dropped.
 *
 * Optional parameters include the boundary type (gene or cds), the number of base pairs to flank the regions, and the output file name.
 * Additionally, users may specify a minimun range size.  Any regions falling shorter than the minimum size will be joined
 * with the shortest adjacent range.   If the user specifies the `makeOnlyGenic` flag, the output will only include genic and CDS regions
 * and minimum range size will be ignored.
 *
 */
class CreateRanges: CliktCommand(help="Create a BED file of reference ranges from a GFF file") {
    val gff by option(help = "GFF file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--gff must not be blank"
            }
        }
    val boundary by option(help = "Reference range boundaries: either gene or cds")
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

    val rangeMinSize by option(help = "Minimum size of a range to be included in the output BED file,  If not specified, a default of 500 bps will be used.")
        .int()
        .default(500)

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
            // GFF is 1-based. We keep 1-based, end-inclusive notation until we make the bed record
            Triple(Position(first.contig,first.position), Position(second.contig,second.position),third)
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
     * Function to merge ranges that are less than rangeMinSize in length
     * This function will merge the current record with either the previous or the next record,
     * whichever is shorter.
     */
    fun mergeShortRanges(bedRecords: List<BedRecord>, minSize:Int) : List<BedRecord> {
        val bedRecordsMerged = mutableListOf<BedRecord>()

        val firstRecord = bedRecords[0]
        bedRecordsMerged.add(firstRecord)

        var previousRecord = bedRecords[0]
        //loop through the rest of the records, checking for short ranges to merge
        var idx = 1
        while (idx  < bedRecords.size) {
            val currentRecord = bedRecords[idx]

            if(currentRecord.contig != previousRecord.contig) {
                //add in new beginning of next chrom
                bedRecordsMerged.add(currentRecord)
                previousRecord = currentRecord
                idx++
            }
            else {
                // check size of current record.  If it is less than rangeMinSize, merge it with the previous record
                // by changing the end of the previous record to the end of the current record.
                // Drop the previous record from the list, and add the merged record to the list
                // We also check the previous record size.  This takes care of an edge case where the
                // first region in a contig is shorter than minSize.  This will also merge the first and second records.
                val currentRecordSize = currentRecord.end - currentRecord.start
                val prevRecordSize = previousRecord.end - previousRecord.start
                val nextRecordSize = if ((idx < bedRecords.size-1) && (currentRecord.contig == bedRecords[idx+1].contig)) bedRecords[idx+1].end - bedRecords[idx+1].start else -1

                // twoBackSize is used when the previous record is too short.  We want to merge to the shorter of
                // the 2 records flanking it.  If it is the first record in a contig, we merge with the current record.
                // By setting the default to one greater than the currentRecordSize we ensure that we merge with the
                // current record if the previous record is too short and it has no prior on this contig.
                val mergedListSize = bedRecordsMerged.size
                val twoBackSize = if (mergedListSize > 1 && bedRecordsMerged[mergedListSize-2].contig == previousRecord.contig) bedRecordsMerged[mergedListSize-2].end - bedRecordsMerged[mergedListSize-2].start else currentRecordSize+1

                if (prevRecordSize < minSize && twoBackSize < currentRecordSize ){
                    // Merge the previous record with the one before it.
                    // We need to drop both the previous record and the one before it, then
                    // add the newly merged record to the merged list.
                    // We hit this condition if we have several short regions in a row. These regions are
                    // merged to each other, but the combined length is still shorter than minSize.
                    // When we finally get to a longer region (currentRecord), we determine the previous region's
                    // predecessor is shorter than the current region, so we merge left.
                    val nameToAdd = "${bedRecordsMerged[mergedListSize-2].name},${previousRecord.name}"
                    val mergedRecord = BedRecord(
                        previousRecord.contig,
                        bedRecordsMerged[mergedListSize-2].start,
                        previousRecord.end,
                        nameToAdd,
                        0,
                        "+")

                    // Need to drop the last 2 records, and add the merged record
                    bedRecordsMerged.removeAt(bedRecordsMerged.size-1)
                    bedRecordsMerged.removeAt(bedRecordsMerged.size-1)
                    bedRecordsMerged.add(mergedRecord)
                    previousRecord = mergedRecord
                    // reset index as we merged the previous record with its previous record
                    // and have not yet processed the current record
                    idx--
                }
                // we merge left if:
                // 1.  previous record is too short, merge that with the current record
                // 2.  current record is too short AND there is no next record in this contig, or the next record is longer than the previous record
                else if ((currentRecordSize < minSize && (nextRecordSize < 0 || nextRecordSize >= prevRecordSize)) || prevRecordSize < minSize) {
                    val nameToAdd = "${previousRecord.name},${currentRecord.name}"
                    val mergedRecord = BedRecord(
                        currentRecord.contig,
                        previousRecord.start,
                        currentRecord.end,
                        nameToAdd,
                        0,
                        "+"
                    )
                    // Drop the previous record from the list
                    bedRecordsMerged.removeAt(bedRecordsMerged.size-1)
                    bedRecordsMerged.add(mergedRecord)
                    previousRecord = mergedRecord
                } else if (currentRecordSize < minSize) {
                    // currentRecordSize is shorter than minSize, there is a next record and its size is shorter than the
                    // previous record size.  In this case, we merge the current record with the next record.
                    val nameToAdd = "${currentRecord.name},${bedRecords[idx+1].name}"
                    val mergedRecord = BedRecord(
                        currentRecord.contig,
                        currentRecord.start,
                        bedRecords[idx+1].end,
                        nameToAdd,
                        0,
                        "+"
                    )
                    bedRecordsMerged.add(mergedRecord)
                    // move past the record we merged into so when we increment idx at the bottom we are
                    // the correct record. (ie, idx needs to be incremented twice since we already processed the next record)
                    idx++
                    previousRecord = mergedRecord
                } else {
                    // add the current record unchanged to the list
                    bedRecordsMerged.add(currentRecord)
                    previousRecord = currentRecord
                }
                idx++
            }
        }

        return bedRecordsMerged
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
            // gffs were 1-based, end inclusive. Bed format is 0-based, end exclusive. So we have to subtract 1 from the lower endpoint now
            bedRecordList.add(BedRecord(range.key.lowerEndpoint().contig, range.key.lowerEndpoint().position-1, range.key.upperEndpoint().position, range.value, 0, "."))
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
     * The trick is when the first record of a contig is too short.  It gets put on
     * the list automatically and we don't start checking size of regions until we process
     * the second region in the chromosome.
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
                // if previous record and current record are not directly adjacent, add an intergenic record
                if (previousRecord.end != currentRecord.start) {
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
            }
            //Either way, add the current record
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

        // Merge ranges that are not rangeMinSize in length
        val bedRecordsMerged = if (!makeOnlyGenic && rangeMinSize > 0) mergeShortRanges(filledInBedRecords,rangeMinSize) else filledInBedRecords

        val bedLinesToPrint = convertBedRecordsIntoOutputStrings(bedRecordsMerged)

        if(output!= null) {
            File(output).bufferedWriter().use { output ->
                bedLinesToPrint.forEach {
                    output.write("$it\n")
                }
            }
        }
    }

}
