package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeader
import net.maizegenetics.phgv2.utils.retrieveAgcContigs
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter


//Making Number a string as VCF allows for '.'
data class AltHeaderMetaData(val id: String, val description:String, val number: String, val source:String,
                             val contig:String, val start:Int, val end:Int, val checksum:String, val refRange:String)
data class HaplotypeSequence(val id: String, val sequence: String, val refRangeId: String, val refContig: String, val refStart: Int, val refEnd: Int,
                             val asmContig : String, val asmStart: Int, val asmEnd: Int)
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
                check(idsToValueMap.containsKey("ID")) { "ALT Header does not contain ID" }
                check(idsToValueMap.containsKey("Description")) { "ALT Header does not contain Description" }
                check(idsToValueMap.containsKey("Number")) { "ALT Header does not contain Number" }
                check(idsToValueMap.containsKey("Source")) { "ALT Header does not contain Source" }
                check(idsToValueMap.containsKey("Contig")) { "ALT Header does not contain Contig" }
                check(idsToValueMap.containsKey("Start")) { "ALT Header does not contain Start" }
                check(idsToValueMap.containsKey("End")) { "ALT Header does not contain End" }
                check(idsToValueMap.containsKey("Checksum")) { "ALT Header does not contain Checksum" }
                check(idsToValueMap.containsKey("RefRange")) { "ALT Header does not contain RefRange" }

                idsToValueMap["ID"]!! to AltHeaderMetaData(idsToValueMap["ID"]!!,idsToValueMap["Description"]!!,idsToValueMap["Number"]!!,idsToValueMap["Source"]!!,
                    idsToValueMap["Contig"]!!,idsToValueMap["Start"]!!.toInt(),idsToValueMap["End"]!!.toInt(),idsToValueMap["Checksum"]!!,idsToValueMap["RefRange"]!!)
            }
    }

    fun createHaplotypeSequences(dbPath:String, sampleName: String, haplotypeVariants: List<VariantContext>, altHeaders: Map<String, AltHeaderMetaData>): List<HaplotypeSequence> {
        return haplotypeVariants.filter { it.hasGenotype(sampleName) }.map {
            val hapId = it.getGenotype(sampleName).getAllele(0).displayString.replace("<","").replace(">","")
            check(altHeaders.containsKey(hapId)) { "Haplotype ID $hapId not found in ALT Header" }
            val altMetaData = altHeaders[hapId]
            //Need to subtract 1 from start as it uses 0 based format
            val ranges = listOf("${altMetaData!!.contig}@${sampleName}:${altMetaData!!.start-1}-${altMetaData!!.end-1}")
            val outputDisplayName = "${altMetaData!!.contig}:${altMetaData!!.start-1}-${altMetaData!!.end-1}"
            val seqs = retrieveAgcContigs(dbPath,ranges)
            HaplotypeSequence(hapId, seqs[outputDisplayName]!!.seq(), altMetaData.refRange, it.contig, it.start, it.end, altMetaData!!.contig, altMetaData!!.start, altMetaData!!.end)
        }
    }

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