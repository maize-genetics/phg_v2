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
 * this command takes a pathing hvcf file and a vcf file created by aligning reads to a composite genome
 * (the latter from a program e.g. DeepVariant or Octopus) and creates a
 * vcf file in haplotype coordinates.
 * The incoming variantVCF file is in composite genome coordinates.
 * The incoming path hvcf file is in reference coordinates.
 * The output file will be in haplotype coordinates, calculated from the Regions information
 * in the path hvcf file's ALT header information, combined with the reference coordinates from the path hvcf entries.
 *
 *
 */
data class HaplotypeData(val id: String, val refContig: String, val refStart: Int, val hapLen: Int)

class CompositeToHaplotypeCoords: CliktCommand(help = "Create g.vcf file for a PHG pathing h.vcf using data from existing PHG created g.vcf files")  {
    private val myLogger = LogManager.getLogger(CompositeToHaplotypeCoords::class.java)
    val pathHvcf by option(help = "Full path to the hvcf file created by the find-paths command ")
        .required()

    val variantVcf by option(help = "Full path to the vcf file created by aligning the reads to the composite genome, then running a variant caller e.g. DeepVariant.  " +
            "\nThe composite genome should have been created by running create-fasta-from-hvcf om the path hvcf file.")
        .required()

    val sampleName by option(help="Sample name to use in the output VCF file")
        .required()

    val outputFile by option(help = "Full path to output VCF file")
        .required()

    override fun run() {

        logCommand(this)

        val (hapidToLength,pathingVCFRangeMap) = processPathVCF(pathHvcf)
        // The pathingVCFRangeMap is to identify the haplotype in which the
        // variantVCF's POS is located.  The haplotype ID is used in the CHROM field
        // of the new VCF file.
        myLogger.info("Calling createHaplotypeVariantContexts")
        val hapVCs = createHaplotypeVariantContexts(variantVcf, pathingVCFRangeMap)
        // hapid length is needed when creating the "contig" headers in the VCF file
        // we need the length of each.
        myLogger.info("Calling writeHaplotypeVCF")
        writeHaplotypeVCF(variantVcf,hapidToLength, hapVCs, sampleName,outputFile)
        // create new outputfile which is the parent of the outputFile parameter with
        // a file name of chromLengths.txt
        val chromLengthsFile = File(outputFile).parent + "/chromLengths.txt"

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
        // These are the ##contig headers
        val sequenceRecordList = hapidToLength.keys.map { SAMSequenceRecord(it,
            hapidToLength[it]!!) }
        vcfHeader.setSequenceDictionary(SAMSequenceDictionary(sequenceRecordList))
        return vcfHeader
    }

    // Function to write the haplotype VCF file
    // The sample name here is the sample name we want used for the NEW haplotype based
    // vcf file.  This is the name that will show up in the sample column of the VCF file.
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
        myLogger.info("Haplotype VCF file written to: $outputFileName")
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
        myLogger.info("\ncreateHaplotypeVariantContexts: Number of records in the variant VCF file: ${records.size}")
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
                // If this haplotype is the first one for the chrom, and it has positions 1-25
                // and the variant is at position 20, then the new position is 20-1+1 = 20
                // We are off-by-1 if we don't add the 1 at the end.
                // or, if the haplotype starts at position 57 in the file, and the SNP is at
                // position 60, then the new position is 60-57+1 = 4 in the haplotype.
                val newPos = pos - hapStart +1
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
        val chromToHapData = mutableMapOf<String, MutableList<HaplotypeData>>()
        val haplotypeVariants = reader.iterator().asSequence().toList()

        val headerIDtoRegionLengthMap = mutableMapOf<String,Int>() // This is the length of the haplotype
        // should only be 1 sample in this file.
        val sampleName = reader.header.sampleNamesInOrder[0]

        haplotypeVariants.filter{it.getGenotype(sampleName).getAllele(0).displayString.replace("<","").replace(">","") != ""}
            .map{
                val hapId = it.getGenotype(sampleName).getAllele(0).displayString.replace("<","").replace(">","")
                // If the hapId from the variants is a non-blank value that is not in the ALT header, throw an exception
                check(altHeaders.containsKey(hapId)) { "Haplotype ID $hapId not found in ALT Header" }
                val altMetaData = altHeaders[hapId]

                val regions =  altMetaData!!.regions
                //Add the haplotype data to the list of haplotype sequences for this chromosome
                val chromList = chromToHapData.getOrDefault(it.contig, mutableListOf())
                val hapLen = regions.map{abs(it.second.position-it.first.position)+1}.sum()
                chromList.add(HaplotypeData(hapId,it.contig,it.start,hapLen))
                chromToHapData[it.contig] = chromList
                headerIDtoRegionLengthMap[hapId] = hapLen
            }

        // Craete a map of haplotype IDs to their positions in the composite genome.
        val hapidToGenomePositions = createHapPositionMap(chromToHapData)
        return Pair(headerIDtoRegionLengthMap,hapidToGenomePositions)
    }

    /**
     * This method splits the chromosome into the haplotypes that make it up.
     * We want to know where each haplotype starts and ends in the composite genome for that chromosome.
     */
    fun createHapPositionMap(chromToAltMap:Map<String,List<HaplotypeData>>):RangeMap<Position,String> {

        // Create a RangeMap of POS to haplotypeID, where POS is a Position object with 1-based start and end
        // coordinates.  The start coordinate is the sum of the lengths of the haplotypes that come before it.
        // The end coordinate is the start coordinate plus the length of the haplotype.
        // This gives us the positions in the composite genome that each haplotype covers.
        val hapidToGenomePositions:RangeMap<Position,String> = TreeRangeMap.create()

        for (chrom in chromToAltMap.keys) {
            var accumulated = 0
            // Need the haplotypes sorted by refStart to ensure the order they appear
            // in the composite fasta for this chromosome is correct
            val haplotypes = chromToAltMap[chrom]!!.sortedBy {it.refStart} // get all the haplotypes for this chrom

            for (hap in haplotypes) {
                // create a range from accumulated to accumulated + headerIDtoRegionLengthMap[hap]
                // put this range in the hapidToGenomePositions map with the hap as the value
                // Positions will be 1-based as the vcf file is 1-based
                accumulated // start at 1 past the last position
                val range = Range.closed(Position(chrom,accumulated+1), Position(chrom,accumulated + hap.hapLen))
                hapidToGenomePositions.put(range,hap.id)
                accumulated += hap.hapLen
            }
        }
        return hapidToGenomePositions
    }

}