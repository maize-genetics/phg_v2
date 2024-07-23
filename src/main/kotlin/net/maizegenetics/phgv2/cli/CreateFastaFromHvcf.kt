package net.maizegenetics.phgv2.cli

import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.convert
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.utils.*
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter


data class HaplotypeSequence(val id: String, val sequence: String, val refRangeId: String, val refContig: String, val refStart: Int, val refEnd: Int,
                             val asmRegions:List<Pair<Position, Position>>)

/**
 * Class to create either a composite or a haplotype fasta file from an input hvcf file or from TileDB directly.
 *
 * Users MAY specify either a folder with multiple hvcf files, or a single hvcf file.  Output files will be
 * written to a specified directory, with names based on the hvcf file name.
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

sealed class HvcfInput {
    data class HvcfFile(val hvcfFile: String): HvcfInput()
    data class HvcfDir(val hvcfDir: String): HvcfInput()
}
class CreateFastaFromHvcf : CliktCommand( help = "Create a FASTA file from a h.vcf file or TileDB directly") {

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        .default("")

    val outputDir by option("-o", "--output-dir", help = "Directory where output fasta Files will be written")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--output-dir/-o must not be blank"
            }
        }

    val fastaType by option("--fasta-type", help = "Type of fasta exported.  Can be either composite or haplotype")
        .default("")
        .validate {
            require(it == "composite" || it == "haplotype") {
                "--fasta-type must be either composite or haplotype"
            }
        }

    // Users may input either a single hvcf file or a directory with multiple hvcf files
    // Currently one of the two options must be provided but this may change when we support
    // processing directly from TileDB
    val hvcfInput: HvcfInput? by mutuallyExclusiveOptions<HvcfInput>(
        option("--hvcf-file", help = "Path to hVCF file. Data will be pulled directly from this file instead of querying TileDB")
            .convert { HvcfInput.HvcfFile(it) },
        option("--hvcf-dir", help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
            .convert { HvcfInput.HvcfDir(it) }
    ).single() // cannot provide both hvcf-file and hvcf-dir at the same time
        .required() // one of the two options must be provided

    val condaEnvPrefix by option (help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
        .default("")

    // Pre-compile the Regex pattern - used when creating the output fasta file names
    val HVCF_PATTERN = Regex("""(\.hvcf|\.h\.vcf|\.hvcf\.gz|\.h\.vcf\.gz)$""")

    /**
     * Function to build the Fasta file from the HVCF and the agc record.
     * Right now it does not support pulling from TileDB, but will in the future.
     */
    fun buildFastaFromHVCF(dbPath: String, outputDir: String, fastaType:String, hvcfInput:HvcfInput?, condaEnvPrefix:String) {

        if (hvcfInput is HvcfInput.HvcfDir) {
            // Loop through the directory and figure out which files are hvcf files
            // The gvcf and hvcf files may be in the same folder, so verify specific extension
            val hvcfFiles = File(hvcfInput.hvcfDir).listFiles { file ->
                HVCF_PATTERN.containsMatchIn(file.name)
            }

            // Loop through each file and run the processSingleHVCF function
            hvcfFiles?.forEach { hvcfFile ->
                val outputFileName = File(hvcfFile.name.replace(HVCF_PATTERN, ".fa"))
                    .name.replace(".fa", "_${fastaType}.fa")
                val outputFile = "$outputDir/$outputFileName"

                BufferedWriter(FileWriter(outputFile)).use { output ->
                    writeSequences(output,
                        processSingleHVCF(VCFFileReader(hvcfFile, false), dbPath, condaEnvPrefix),
                        fastaType)
                }
            }

        } else if (hvcfInput is HvcfInput.HvcfFile){
            // Load in the HVCF
            val hvcfFileReader = VCFFileReader(File(hvcfInput.hvcfFile), false)
            val outputFileName = File(File(hvcfInput.hvcfFile).name.replace(HVCF_PATTERN, ".fa"))
                .name.replace(".fa", "_${fastaType}.fa")
            val outputFile = "$outputDir/$outputFileName"

            BufferedWriter(FileWriter(outputFile)).use { output ->
                val records = processSingleHVCF(hvcfFileReader, dbPath, condaEnvPrefix)
                writeSequences(output, records, fastaType)
            }
        }  else {
            // Load in the TileDB
            TODO("TileDB VCF Reader Not implemented yet.  Please run with --hvcf-file or --hvcf-dir")
        }

    }

    fun processSingleHVCF(vcfFileReader: VCFFileReader, dbPath: String,condaEnvPrefix:String) : List<HaplotypeSequence> {
        //extract out the haplotype sequence boundaries for each haplotype from the hvcf
        val altHeaderMap = parseALTHeader(vcfFileReader.header)

        val samples = vcfFileReader.header.sampleNamesInOrder
        val hvcfRecords = vcfFileReader.iterator().asSequence().toList()

        return samples.flatMap { sample -> createHaplotypeSequences(dbPath, sample, hvcfRecords, altHeaderMap,condaEnvPrefix) }
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
    fun createHaplotypeSequences(dbPath:String, sampleName: String, haplotypeVariants: List<VariantContext>, altHeaders: Map<String, AltHeaderMetaData>,condaEnvPrefix:String): List<HaplotypeSequence> {
        val chromToAltEntryData = mutableMapOf<String, MutableList<Triple<MutableList<String>,MutableList<String>,HaplotypeSequence>>>()
        val hapSeqList = mutableListOf<HaplotypeSequence>()
        // Filter out the haplotype variants that do not have the sampleName
        // Also filter out the haplotype variants that have a blank haplotype ID
        // FindPaths includes missing haplotypes in the h.vcf file.  These are represented via "." in the h.vcf file and appear
        // as empty string in the variable "hapId" below.  This is not an error
        haplotypeVariants.filter { it.hasGenotype(sampleName) }
            .filter{it.getGenotype(sampleName).getAllele(0).displayString.replace("<","").replace(">","") != ""}
            .map {
            val hapId = it.getGenotype(sampleName).getAllele(0).displayString.replace("<","").replace(">","")
            // If the hapId from the variants is a non-blank value that is not in the ALT header, throw an exception
            check(altHeaders.containsKey(hapId)) { "Haplotype ID $hapId not found in ALT Header" }
            val altMetaData = altHeaders[hapId]
            val hapSampleName = altMetaData!!.sampleName()
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
            //Add the haplotype sequence to the list of haplotype sequences for this chromosome
            val chromList = chromToAltEntryData.getOrDefault(it.contig, mutableListOf())
            chromList.add(Triple(queryRanges, displayRanges, HaplotypeSequence(hapId, "", altMetaData.refRange, it.contig, it.start, it.end, regions)))
            chromToAltEntryData[it.contig] = chromList
        }

        // The full set of ranges in a single query is too much for both ProcessBuilder
        // and for AGC.  Loop over the chromosomes, grabbing sequence from AGC a chrom at a time.
        // This returns a list of HaplotypeSequence objects.
        for(chrom in chromToAltEntryData.keys) {
            val ranges = chromToAltEntryData[chrom]!!.flatMap { it.first }
            val seqs = retrieveAgcContigs(dbPath,ranges,condaEnvPrefix)
            // add to the hapSeqList an updated HaplotypeSequence object containing the sequence built from AGC contigs
             hapSeqList.addAll(chromToAltEntryData[chrom]!!.map { it.third.copy(sequence = buildHapSeq(seqs, it.second,it.third)) })
        }
        return hapSeqList
    }

    /**
     * Function to build the haplotype sequence based on the list of display regions and the given haplotype sequence object.
     * The sequence has already extracted out of AGC and stored in the seqs map.
     * The incoming "seqs" parameter has the key as a Pair(sampleName,displayRegion), where display region could
     * be just a contig, or a contig with ranges:  e.g. "chr1" or "chr1:100-200".  SampleName is the sample from which
     * the sequence was extracted.
     */
    fun buildHapSeq(seqs: Map<Pair<String,String>,NucSeq> , displayRegions : List<String>, hapSeqObjects: HaplotypeSequence) : String {
        // hapSeqRegions is the HaplotypeSequence object's list of regions.
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
        // Check the dbPath and set it to the current working directory if it is not provided
        val dbPath = if (dbPath.isBlank()) {
            System.getProperty("user.dir")
        } else {
            dbPath
        }

        // Verify the tiledbURI
        // If it doesn't an exception will be thrown
        val validDB = verifyURI(dbPath,"hvcf_dataset",condaEnvPrefix)

        buildFastaFromHVCF(dbPath, outputDir, fastaType, hvcfInput, condaEnvPrefix)
    }
}