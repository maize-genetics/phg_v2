package net.maizegenetics.phgv2.cli

import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.utils.AltHeaderMetaData
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.parseALTHeader
import net.maizegenetics.phgv2.utils.retrieveAgcContigs
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter


data class HaplotypeSequence(val id: String, val sequence: String, val refRangeId: String, val refContig: String, val refStart: Int, val refEnd: Int,
                             val asmRegions:List<Pair<Position, Position>>)

/**
 * Class to create either a composite or a haplotype fasta file from an input hvcf file or from TileDB directly.
 *
 * As mentioned in the terminology section of the Documentation, a composite fasta file is a fasta file that contains
 * all the haplotypes for a given contig concatenated together by contig.  This pseudo genome can be used for rare
 * allele finding.
 *
 * A Haplotype fasta file is where we output each haplotype as a separate fasta entry.  This can be used for read
 * mapping purposes and imputation purposes or simple haplotype sequence retrieval.
 *
 * The class can be used in two ways.  The first is to create a fasta file from a hvcf file.  The second is to create a
 * fasta file from TileDB directly.  The first method is useful for creating a fasta file from a hvcf file that is
 * already created.  The second method is useful for creating a fasta file directly from a TileDB database to avoid
 * having to create an intermediate hvcf file first.
 *
 *
 * Note as of now, the class only supports creating a fasta file from a hvcf file.  The TileDB functionality will be
 * added in the future.
 */
class CreateFastaFromHvcf : CliktCommand( help = "Create a fasta file from a hvcf file/TileDB") {

    val dbPath by option(help = "Tile DB URI")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--db-path must not be blank"
            }
        }

    val output by option("-o", "--output", help = "Name for output Fasta file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--output/-o must not be blank"
            }
        }

    val fastaType by option("--fasta-type", help = "Type of fasta exported.  Can be either composite or haplotype")
        .default("")
        .validate {
            require(it == "composite" || it == "haplotype") {
                "--fasta-type must be either composite or haplotype"
            }
        }

    val hvcfFile by option("--hvcf-file", help = "Path to hVCF file. Data will be pulled directly from this file instead of querying TileDB")
        .default("")

    val hvcfDir by option("--hvcf-dir", help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        .default("")

    /**
     * Function to build the Fasta file from the HVCF and the agc record.
     * Right now it does not support pulling from TileDB, but will in the future.
     */
    fun buildFastaFromHVCF(dbPath: String, outputFile: String, fastaType:String, hvcfDir: String ,hvcfFile : String) {
        if(hvcfDir != "") {
            //Loop through the directory and figure out which files are hvcf files
            val hvcfFiles = File(hvcfDir).listFiles { file -> file.extension == "vcf" || file.name.endsWith("vcf.gz") }
            //Loop through each file and run the buildFastaFromHVCF function
            BufferedWriter(FileWriter(outputFile)).use { output ->
                writeSequences(output,hvcfFiles.flatMap { processSingleHVCF(VCFFileReader(it,false), dbPath) }
                    .associate { it.id to it }
                    .values.toList(), fastaType)
            }
        }
        else  if(hvcfFile == "") {
            //Load in the TileDB
            TODO("TileDB VCF Reader Not implemented yet.  Please run with --hvcf-file")
        } else {
            //Load in the HVCF
            val hvcfFile = VCFFileReader(File(hvcfFile),false)
            BufferedWriter(FileWriter(outputFile)).use { output ->
                val records = processSingleHVCF(hvcfFile, dbPath)
                writeSequences(output,records, fastaType)
            }
        }
    }

    fun processSingleHVCF(vcfFileReader: VCFFileReader, dbPath: String) : List<HaplotypeSequence> {
        //extract out the haplotype sequence boundaries for each haplotype from the hvcf
        val altHeaderMap = parseALTHeader(vcfFileReader.header)

        val samples = vcfFileReader.header.sampleNamesInOrder
        val hvcfRecords = vcfFileReader.iterator().asSequence().toList()

        return samples.flatMap { sample -> createHaplotypeSequences(dbPath, sample, hvcfRecords, altHeaderMap) }
    }

    fun writeSequences(outputWriter: BufferedWriter, haplotypeSequences: List<HaplotypeSequence>, fastaType: String) {
        if(fastaType == "composite")
            writeCompositeSequence(outputWriter, haplotypeSequences)
        else if(fastaType == "haplotype") {
            writeHaplotypeSequence(outputWriter, haplotypeSequences)
        }
    }

    /**
     * Function to create haplotype Sequences for each of the haplotype variants in the hvcf
     * Currently, sampleName is from a single HVCF file, and is the sample from which the haplotype sequences will be
     * extracted.
     *
     * In this function,
     *    "sampleName" parameter is the sampleName from the headerline of the hvcf file
     *    "hapSampleName" is the samplename associated with the haplotype sequence, and this information
     *           is pulled from the ALT header line using the hapid as an index.
     */
    fun createHaplotypeSequences(dbPath:String, sampleName: String, haplotypeVariants: List<VariantContext>, altHeaders: Map<String, AltHeaderMetaData>): List<HaplotypeSequence> {
        val rangesAndOtherInfo = haplotypeVariants.filter { it.hasGenotype(sampleName) }.map {
            val hapId = it.getGenotype(sampleName).getAllele(0).displayString.replace("<","").replace(">","")
            check(altHeaders.containsKey(hapId)) { "Haplotype ID $hapId not found in ALT Header" }
            val altMetaData = altHeaders[hapId]
            val hapSampleName = altMetaData!!.sampleName
            //Need to subtract 1 from start as it uses 0 based format
            val regions =  altMetaData!!.regions
            val queryRanges = mutableListOf<String>()
            val displayRanges = mutableListOf<String>()
            for(region in regions) {
                if(region.first.position-1 > region.second.position-1) {
                    queryRanges.add("${region.first.contig}@${hapSampleName}:${region.second.position-1}-${region.first.position-1}")
                    displayRanges.add("${hapSampleName}@${region.first.contig}:${region.second.position-1}-${region.first.position-1}")
                }
                else {
                    queryRanges.add("${region.first.contig}@${hapSampleName}:${region.first.position - 1}-${region.second.position - 1}")
                    displayRanges.add("${hapSampleName}@${region.first.contig}:${region.first.position-1}-${region.second.position-1}")
                }
            }
            Triple(queryRanges, displayRanges, HaplotypeSequence(hapId, "", altMetaData.refRange, it.contig, it.start, it.end, regions))
        }

        //Create a list of ranges we need to extract.
        //Go through and split them into sets of non-identical ranges
        //Pull all the non-identical ranges in a multithreaded fashion
        //Store in Map<Range,Sequence>
        val ranges = rangesAndOtherInfo.flatMap { it.first }
        val seqs = retrieveAgcContigs(dbPath,ranges)

        return rangesAndOtherInfo.map { it.third.copy(sequence = buildHapSeq(seqs, it.second,it.third)) }
    }

    /**
     * Function to build the haplotype sequence based on the list of display regions and the given haplotype sequence object.
     * The sequence has already extracted out of AGC and stored in the seqs map.
     * The incoming "seqs" parameter has the key as a Pair(sampleName,displayRegion), where display region could
     * be just a contig, or a contig with ranges:  e.g. "chr1" or "chr1:100-200".  SampleName is the sample from which
     * the sequence was extracted.
     */
    fun buildHapSeq(seqs: Map<Pair<String,String>,NucSeq> , displayRegions : List<String>, hapSeqObjects: HaplotypeSequence) : String {
        // hapSeqRegions is the HapltoypeSEquence object's list of regions.
        // This was passed in as part of the Triple created in the calling method.
        // Because of the way these were created, the displayRegions and the hapSeqRegions should be in the same order.
        val hapSeqRegions = hapSeqObjects.asmRegions

        // This gets all the sequences from all the regions in the list,
        // and joins them to a string with no separator.  The string that is
        // returned is the sequence for the haplotype.
        return displayRegions.mapIndexed{ idx, currentDisplayRegion ->
            val currentHapSeqRegion = hapSeqRegions[idx]

            // The displayRegions are of the form: sampleName@contig:stPos-endPos
            val sampleName = currentDisplayRegion.split("@")[0]
            val region = currentDisplayRegion.split("@")[1]

            val seq = seqs[Pair(sampleName,region)]!!

            //Check to see if we have an inverted sub region based on the currentHapSeqRegion
            if(currentHapSeqRegion.first.position > currentHapSeqRegion.second.position) {
                //If so we need to reverse compliment the sequence
                seq.reverse_complement().seq()
            }
            else {
                seq.seq()
            }
        }.joinToString("")

    }

    /**
     * Function to output composite contig sequences to a fasta file.
     * Here a composite sequence is a pseudo genome where each haplotype sequence for a given chromosome is concatenated together.
     * Note no Ns are added between the haplotype sequences.
     */
    fun writeCompositeSequence(outputFileWriter: BufferedWriter, haplotypeSequences: List<HaplotypeSequence>) {
        //group the sequences by chromosome
        val sequencesByChr = haplotypeSequences.groupBy { it.refContig }
        for(chr in sequencesByChr.keys.sorted()) {
            outputFileWriter.write(">$chr\n")
            //sort the sequences by startPos
            val sequencesByStartPos = sequencesByChr[chr]!!.sortedBy { it.refStart }
            //merge and output the sequences
            sequencesByStartPos.map { it.sequence }
                .joinToString("")
                .chunked(80)//Chunking into 80 character lines
                .forEach { outputFileWriter.write(it + "\n") }
        }
    }

    /**
     * Function to output haplotype sequences to a fasta file.  Here each haplotype is exported as its own fasta record without concatenating things together.
     * This is almost identical to how fastas were exported in the original version of the pipeline.
     */
    fun writeHaplotypeSequence(outputFileWriter: BufferedWriter, haplotypeSequences: List<HaplotypeSequence>, exportFullIdLine : Boolean = true) {
        for(hapSeq in haplotypeSequences) {
            outputFileWriter.write(">${hapSeq.id}")
            if(exportFullIdLine) {
                outputFileWriter.write(" Ref_Range_Id=${hapSeq.refRangeId} " +
                "Ref_Contig=${hapSeq.refContig} Ref_Start=${hapSeq.refStart} Ref_End=${hapSeq.refEnd} " +
                        "Asm_Regions=${hapSeq.asmRegions.joinToString(",") { "${it.first.contig}:${it.first.position}-${it.second.position}" }}")
            }
            outputFileWriter.write("\n")
            hapSeq.sequence
                .chunked(80)//Chunking into 80 character lines
                .forEach { outputFileWriter.write(it + "\n") }
        }
    }

    override fun run() {
        buildFastaFromHVCF(dbPath, output, fastaType, hvcfDir, hvcfFile)
    }
}