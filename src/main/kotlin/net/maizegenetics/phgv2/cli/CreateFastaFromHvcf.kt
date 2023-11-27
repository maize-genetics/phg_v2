package net.maizegenetics.phgv2.cli

import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeader
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.retrieveAgcContigs
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter


//Making Number a string as VCF allows for '.'
data class AltHeaderMetaData(val id: String, val description:String, val number: String, val source:String,
                             val regions:List<Pair<Position, Position>>, val checksum:String, val refRange:String)
data class HaplotypeSequence(val id: String, val sequence: String, val refRangeId: String, val refContig: String, val refStart: Int, val refEnd: Int,
                             val asmRegions:List<Pair<Position, Position>>)

/**
 * Class to create either a composite or a haplotype fasta file from an input hvcf file or from TileDB directly.
 *
 * As mentioned in the terminology section of the Documentation, a composite fasta file is a fasta file that contains
 * all of the haplotypes for a given contig concatenated together by contig.  This pseudo genome can be used for rare
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
     * Helper function to parse out the ALT headers from the VCF file.
     *
     * We need to do a bit more involved parsing in this function as we cannot use the .getOtherHeaders() call from HTSJDK.
     * For some reason this only returns the first header when called and we need all of them.
     * The work around is that we can get all the metadata, filter out any that are not ALT then parse the ALT header using normal string parsing.
     * To make this easy, we just parse each piece of metadata into a key-value pair and then store in a map.
     */
    fun parseALTHeader(header: VCFHeader) : Map<String, AltHeaderMetaData> {
        //Need to turn the ALT File header into a Map<ID,AltHeaderMetaData>
        return header.metaDataInInputOrder
            .filter { it.key == "ALT" }
            .map { it.toString()
                .substringAfter("<")
                .substringBeforeLast(">")
            } //Keep the useful part of the ALT Tag
            .map { it.split(",") }
            .associate {
                val idsToValueMap = it.map { token -> token.split("=") }.associate { token -> token[0] to token[1] }
                //ID and Description are required fields by VCF spec, if these errors are thrown there is something wrong with the htsjdk library
                check(idsToValueMap.containsKey("ID")) { "ALT Header does not contain ID" }
                check(idsToValueMap.containsKey("Description")) { "ALT Header does not contain Description" }
                //These are optional header fields so we check these in the unit test.
                check(idsToValueMap.containsKey("Number")) { "ALT Header does not contain Number" }
                check(idsToValueMap.containsKey("Source")) { "ALT Header does not contain Source" }
                check(idsToValueMap.containsKey("Regions")) { "ALT Header does not contain Regions" }
                check(idsToValueMap.containsKey("Checksum")) { "ALT Header does not contain Checksum" }
                check(idsToValueMap.containsKey("RefRange")) { "ALT Header does not contain RefRange" }

                idsToValueMap["ID"]!! to AltHeaderMetaData(idsToValueMap["ID"]!!,idsToValueMap["Description"]!!,idsToValueMap["Number"]!!,idsToValueMap["Source"]!!,
                    parseRegions(idsToValueMap["Regions"]!!),idsToValueMap["Checksum"]!!,idsToValueMap["RefRange"]!!)
            }
    }

    /**
     * Function to parse the regions from the ALT header.
     */
    fun parseRegions(regions: String) : List<Pair<Position, Position>> {
        return regions.split(",").map { it.split(":") }.map {
            val positions = it[1].split("-").map { position -> position.toInt() }
            check(positions.size == 2) { "Region ${it} is not in the correct format.  It needs to be in the form: chr:stPos-endPos." }
            Pair(Position(it[0],positions[0]), Position(it[0],positions[1]))
        }
    }

    /**
     * Function to create haplotype Sequences for each of the haplotype variants in the hvcf
     */
    fun createHaplotypeSequences(dbPath:String, sampleName: String, haplotypeVariants: List<VariantContext>, altHeaders: Map<String, AltHeaderMetaData>): List<HaplotypeSequence> {
        val rangesAndOtherInfo = haplotypeVariants.filter { it.hasGenotype(sampleName) }.map {
            val hapId = it.getGenotype(sampleName).getAllele(0).displayString.replace("<","").replace(">","")
            check(altHeaders.containsKey(hapId)) { "Haplotype ID $hapId not found in ALT Header" }
            val altMetaData = altHeaders[hapId]
            //Need to subtract 1 from start as it uses 0 based format
            val regions =  altMetaData!!.regions
            val queryRanges = mutableListOf<String>()
            val displayRanges = mutableListOf<String>()
            for(region in regions) {
                if(region.first.position-1 > region.second.position-1) {
                    queryRanges.add("${region.first.contig}@${sampleName}:${region.second.position-1}-${region.first.position-1}")
                    displayRanges.add("${region.first.contig}:${region.second.position-1}-${region.first.position-1}")
                }
                else {
                    queryRanges.add("${region.first.contig}@${sampleName}:${region.first.position - 1}-${region.second.position - 1}")
                    displayRanges.add("${region.first.contig}:${region.first.position-1}-${region.second.position-1}")
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
     * The sequence is already extracted out of AGC and stored in the seqs map.
     */
    fun buildHapSeq(seqs: Map<String,NucSeq> ,displayRegions : List<String>, hapSeqObjects: HaplotypeSequence) : String {
        val hapSeqRegions = hapSeqObjects.asmRegions

        return displayRegions.mapIndexed{ idx, currentDisplayRegion ->
            val currentHapSeqRegion = hapSeqRegions[idx]

            val seq = seqs[currentDisplayRegion]!!

            //Means it is the first region
            if(currentHapSeqRegion.first.position > currentHapSeqRegion.second.position) {
                //Means it needs to be reverse complemented
                seq.reverse_complement().seq()
            }
            else {
                seq.seq()
            }
        }.joinToString()

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