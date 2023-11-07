package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.utils.AltHeaderMetaData
import net.maizegenetics.phgv2.utils.parseALTHeader
import net.maizegenetics.phgv2.utils.retrieveAgcContigs
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter


data class HaplotypeSequence(val id: String, val sequence: String, val refRangeId: String, val refContig: String, val refStart: Int, val refEnd: Int,
                             val asmContig : String, val asmStart: Int, val asmEnd: Int)

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

    val hvcfFile by option("--hvcf-file", help = "hVCF file to import instead of hitting TileDB")
        .default("")

    /**
     * Function to build the Fasta file from the HVCF and the agc record.
     * Right now it does not support pulling from TileDB, but will in the future.
     */
    fun buildFastaFromHVCF(dbPath: String, outputFile: String, fastaType:String, hvcfFile : String) {
        val vcfFileReader = if(hvcfFile == "") {
            //Load in the TileDB
            TODO("TileDB VCF Reader Not implemented yet.  Please run with --hvcf-file")
        } else {
            //Load in the HVCF
            VCFFileReader(File(hvcfFile),false)
        }

        //extract out the haplotype sequence boundaries for each haplotype from the hvcf
        val altHeaderMap = parseALTHeader(vcfFileReader.header)

        val samples = vcfFileReader.header.sampleNamesInOrder
        val hvcfRecords = vcfFileReader.iterator().asSequence().toList()
        //loop through each sample and walk through the hvcf
        for(sample in samples) {
            val haplotypeSequences = createHaplotypeSequences(dbPath, sample,hvcfRecords,altHeaderMap)
            //write out the haplotype sequence to a fasta file
            if(fastaType == "composite")
                writeCompositeSequence(outputFile, haplotypeSequences)
            else if(fastaType == "haplotype") {
                writeHaplotypeSequence(outputFile, haplotypeSequences)
            }
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
            val range = "${altMetaData!!.contig}@${sampleName}:${altMetaData!!.start-1}-${altMetaData!!.end-1}"
            val outputDisplayName = "${altMetaData!!.contig}:${altMetaData!!.start-1}-${altMetaData!!.end-1}"
            Triple(range, outputDisplayName, HaplotypeSequence(hapId, "", altMetaData.refRange, it.contig, it.start, it.end, altMetaData!!.contig, altMetaData!!.start, altMetaData!!.end))
        }

        val ranges = rangesAndOtherInfo.map { it.first }
        val seqs = retrieveAgcContigs(dbPath,ranges)

        return rangesAndOtherInfo.map { it.third.copy(sequence = seqs[it.second]!!.seq()) }
    }

    /**
     * Function to output composite contig sequences to a fasta file.
     * Here a composite sequence is a pseudo genome where each haplotype sequence for a given chromosome is concatenated together.
     * Note no Ns are added between the haplotype sequences.
     */
    fun writeCompositeSequence(outputFile: String, haplotypeSequences: List<HaplotypeSequence>) {
        BufferedWriter(FileWriter(outputFile)).use { output ->
            //group the sequences by chromosome
            val sequencesByChr = haplotypeSequences.groupBy { it.refContig }
            for(chr in sequencesByChr.keys.sorted()) {
                output.write(">$chr\n")
                //sort the sequences by startPos
                val sequencesByStartPos = sequencesByChr[chr]!!.sortedBy { it.refStart }
                //merge and output the sequences
                sequencesByStartPos.map { it.sequence }
                    .joinToString("")
                    .chunked(80)//Chunking into 80 character lines
                    .forEach { output.write(it + "\n") }
            }
        }
    }

    /**
     * Function to output haplotype sequences to a fasta file.  Here each haplotype is exported as its own fasta record without concatenating things together.
     * This is almost identical to how fastas were exported in the original version of the pipeline.
     */
    fun writeHaplotypeSequence(outputFile: String, haplotypeSequences: List<HaplotypeSequence>, exportFullIdLine : Boolean = true) {
        BufferedWriter(FileWriter(outputFile)).use { output ->
            for(hapSeq in haplotypeSequences) {
                output.write(">${hapSeq.id}")
                if(exportFullIdLine) {
                    output.write(" Ref_Range_Id=${hapSeq.refRangeId} " +
                    "Ref_Contig=${hapSeq.refContig} Ref_Start=${hapSeq.refStart} Ref_End=${hapSeq.refEnd} " +
                            "Asm_Contig=${hapSeq.asmContig} Asm_Start=${hapSeq.asmStart} Asm_End=${hapSeq.asmEnd}")
                }
                output.write("\n")
                hapSeq.sequence
                    .chunked(80)//Chunking into 80 character lines
                    .forEach { output.write(it + "\n") }
            }

        }
    }

    override fun run() {
        buildFastaFromHVCF(dbPath, output, fastaType,hvcfFile)
    }
}