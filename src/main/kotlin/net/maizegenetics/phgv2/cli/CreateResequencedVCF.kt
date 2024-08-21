package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import com.google.common.collect.TreeRangeMap
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.SAMSequenceRecord
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeader
import htsjdk.variant.vcf.VCFHeaderLine
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.parseALTHeader
import org.apache.logging.log4j.LogManager
import java.io.File
import kotlin.math.abs

/**
 * this command takes a pathing hvcf file and a vcf file created by aligning reads to a composite genome and creates a
 * vcf file in hapotype coordinates.
 * The incoming variantVCF file is in composite genome coordinates.
 * The incoming path hvcf file is in reference coordinates.
 * The output file will be in haplotype coordinates, calculated from the Regions information
 * in the path hvcf file's ALT header information, combined with the reference coordinates from the path hvcf entries.
 *
 *
 */

class CreateResequencedVCF: CliktCommand(help = "Create g.vcf file for a PHG pathing h.vcf using data from existing PHG created g.vcf files")  {
    private val myLogger = LogManager.getLogger(CreateResequencedVCF::class.java)
    val pathHvcf by option(help = "Full path to the hvcf file created by the find-paths command ")
        .required()

    val variantVcf by option(help = "Full path to the vcf file created by aligning the reads to the composite genome.  " +
            "\nThe composite genome should have been created by running create-fasta-from-hvcf om the path hvcf file.")
        .required()

    val sampleName by option(help="Sample name to use in the output VCF file")
        .required()

    val outputFile by option(help = "Full path to output VCF file")
        .required()

    override fun run() {
        val (hapidToLength,pathingVCFRangeMap) = processPathVCF(pathHvcf)
        val hapVCs = createHaplotypeVariantContexts(variantVcf, pathingVCFRangeMap)
        writeHaplotypeVCF(variantVcf,hapidToLength, hapVCs, sampleName,outputFile)

    }

    fun reseqVCFHeader(reader:VCFFileReader, hapidToLength:Map<String,Int>,sampleNames:List<String>): VCFHeader {
        // create the headers
        // From the variantVCF, which is a composite genome VCF file, get all the headers
        // EXCEPT the contig headers.  This includes "other" e.g. fileformat=VCFv4.2,
        // ##DeepVariant_version=1.6.0 and lines the show commands run to filter the vcf, etc, which
        // indicate how the vcf was made.  We want to preserve all of these headers, just not the contigs
        // as those now become the hapids.
        val header = reader.fileHeader
        val newHeaderLines = mutableSetOf<VCFHeaderLine>()
        newHeaderLines.addAll(header.filterLines)
        newHeaderLines.addAll(header.formatHeaderLines)
        newHeaderLines.addAll(header.infoHeaderLines)
        newHeaderLines.addAll(header.otherHeaderLines)
        val vcfHeader = VCFHeader(newHeaderLines, sampleNames)

        // Add the hapids from the pathVCF
        val sequenceRecordList = hapidToLength.keys.map { SAMSequenceRecord(it,
            hapidToLength[it]!!) }
        vcfHeader.setSequenceDictionary(SAMSequenceDictionary(sequenceRecordList))
        return vcfHeader
    }

    // Function to write the haplotype VCF file
    fun writeHaplotypeVCF(variantVCF:String, hapidToLength:Map<String,Int>,hapVCs:List<VariantContext>,sampleName:String, outputFileName:String) {
        // setup the writer
        val writer = VariantContextWriterBuilder()
            .unsetOption(Options.INDEX_ON_THE_FLY)
            .setOutputFile(File(outputFileName))
            .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
            .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
            .build()

        //read the original variant vcf files, create the new headers from the original
        val reader = VCFFileReader(File(variantVCF), false)
        val vcfHeader = reseqVCFHeader(reader,hapidToLength,listOf(sampleName))
        writer.writeHeader(vcfHeader)

        // Add the variant lines, created and adjusted for the haplotype coordinates
        for (vc in hapVCs) {
            writer.add(vc)
        }
        writer.close()
    }

    // Take a variant VCF file and a RangeMap of the regions that are in the pathing hvcf file.
    // The variant VCF file is in composite genome coordinates.
    // THe pathingVCFRangeMap has been created in composite genome coordinates.
    // Write a new vcf file with the haplotype id in the CHROM field and the haplotype-relative
    // positions in the POS field.
    fun createHaplotypeVariantContexts(variantVcf: String, pathingVCFRangeMap: RangeMap<Position, String>):List<VariantContext> {
        val reader = VCFFileReader(File(variantVcf), false)
        val records = reader.iterator().asSequence().toList()
        val haplotypeVariantContexts = mutableListOf<VariantContext>()
        for (record in records) {
            val chrom = record.contig
            val pos = record.start
            val range = Position(chrom,pos)
            val entry = pathingVCFRangeMap.getEntry(range)

            if (entry != null) {
                // create a new record with the haplotype id in the CHROM field
                // and the position in the POS field
                // and the rest of the fields as they are in the record
                val hapStart = entry.key.lowerEndpoint().position
                val hapId = entry.value
                val newPos = pos - hapStart
                // Create a VariantContext record that has the haplotype id in the CHROM field
                // and the newPos in the POS field
                // and the rest of the fields as they are in the record

                // Create a new VariantContext with the updated CHROM and POS
                val newRecord = VariantContextBuilder(record)
                    .chr(hapId)   // Set the new CHROM value
                    .start(newPos.toLong()) // Set the new POS value
                    .stop(newPos.toLong()) // Adjust stop position accordingly - SNP, so same as start
                    .make()

                haplotypeVariantContexts.add(newRecord)

            }
        }
        return haplotypeVariantContexts
    }

    // Take a pathing hvcf file and create a RangeMap of the regions that are in the pathing hvcf file.
    // This range map needs to be in composite genome coordinates.
    // We want the RangeMap of Positions->haplotypeID, but we also need a map of
    // haplotypeID to the length of the haplotype.  THis is used for the SequenceDictionary
    // when writing the new VCF file.  IE, these are the ##contig headers
    //
    fun processPathVCF(pathHvcf: String): Pair<Map<String,Int>,RangeMap<Position, String>> {
        // get the header information from the path hvcf file
        val reader = VCFFileReader(File(pathHvcf), false)
        val altHeaders = parseALTHeader(reader.header)
        // You want altHeaders: map of alt header id to len(AltHeaderMetaData.regions)
        // create a map of altHeader.id to the sum of the length of the regions in the altHeaderMetaData

        val headerIDtoRegionLengthMap = mutableMapOf<String,Int>()
        for (altHeader in altHeaders) {
            // Why are there some haplotypes with size=1?  In some cases, the start > end - for reverse strand?
            val regionLength = altHeader.value.regions.map { abs(it.second.position - it.first.position) + 1 }.sum()
            headerIDtoRegionLengthMap[altHeader.key] = regionLength
        }


        val hvcfRecords = reader.iterator().asSequence().toList()
        // Walk through the path hvcf records, create a map of chrom to List<ALT> which is a map
        // of the value in the CHROM field to a List of all the values in the ALT field that
        // show up for that chrom.  If the value in the ALT field is ".", skip it.  The LIst
        // must be in the order that the ALT values show up in the file.  These are the haplotypes
        // whose sequence make up the composite genome at each chromosome.
        val chromToAltMap = mutableMapOf<String, List<String>>()
        for (record in hvcfRecords) {
            val chrom = record.contig
            // set val alt to the first value in the ALT field, and to "." if there is no value
            var alt = if (record.alternateAlleles != null && record.alternateAlleles.size > 0) record.alternateAlleles.get(0).displayString else "."
            if (alt == ".") {
                continue
            }

            // remove a leading "<" and trailing ">" from the alt value
            if (alt.startsWith("<") && alt.endsWith(">")) {
                alt = alt.substring(1,alt.length-1)
            }
            if (chromToAltMap.containsKey(chrom)) {
                chromToAltMap[chrom] = chromToAltMap[chrom]!!.plus(alt)
            } else {
                chromToAltMap[chrom] = listOf(alt)
            }
        }

        // Craete a map of haplotype IDs to their positions in the composite genome.
        val hapidToGenomePositions = createHapPositionMap(chromToAltMap,headerIDtoRegionLengthMap)

        return Pair(headerIDtoRegionLengthMap,hapidToGenomePositions)
    }

    fun createHapPositionMap(chromToAltMap:Map<String,List<String>>,headerIdToRegionLengthMap:Map<String,Int>):RangeMap<Position,String> {

        // Create a RangeMap of POS to haplotypeID, where POS is a Position object with 1-based start and end
        // coordinates.  The start coordinate is the sum of the lengths of the haplotypes that come before it.
        // The end coordinate is the start coordinate plus the length of the haplotype.
        // This gives us the positions in the composite genome that each haplotype covers.
        val hapidToGenomePositions:RangeMap<Position,String> = TreeRangeMap.create()
        for (chrom in chromToAltMap.keys) {
            var accumulated = 0
            val haplotypes = chromToAltMap[chrom]!!
            for (hap in haplotypes) { // Does this process them in order?
                // create a range from accumulated to accumulated + headerIDtoRegionLengthMap[hap]
                // put this range in the hapidToGenomePositions map with the hap as the value
                val hapLen = headerIdToRegionLengthMap[hap]!!
                if (hapLen <= 0) {
                    myLogger.error("Haplotype $hap has length 0")
                    continue
                }
                val range = Range.closed(Position(chrom,accumulated+1), Position(chrom,accumulated+1 + hapLen))
                hapidToGenomePositions.put(range,hap)
                accumulated += hapLen
            }
        }
        return hapidToGenomePositions
    }
}