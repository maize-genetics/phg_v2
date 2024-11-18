package net.maizegenetics.phgv2.utils

import biokotlin.seq.NucSeq
import biokotlin.util.bufferedReader
import com.google.common.collect.Range
import com.google.common.collect.RangeSet
import com.google.common.collect.TreeRangeSet
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.SAMSequenceRecord
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.*
import org.apache.logging.log4j.LogManager
import java.io.File
import java.nio.file.Files
import java.nio.file.Paths
import java.security.MessageDigest
import java.security.NoSuchAlgorithmException
import java.util.*


private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.utils.VariantLoadingUtils")

// We could use the BioKotlin SeqPosition, but it is heavier than we need
// as it contains a NucSeqRecord to get the chromosome name, and that includes sequence
// It needs to be Comparable to be used in a RangeSet
data class Position (val contig: String, val position: Int) : Comparable<Position> {
    override fun compareTo(other: Position): Int {
        if (this.contig == other.contig) {
            return this.position.compareTo(other.position)
        }
        return this.contig.compareTo(other.contig)
    }

    override fun toString(): String {
        return "$contig:$position"
    }
}

/**
 * Function to write out the Variant Contexts to a file.
 */
fun exportVariantContext(sampleName: String, variantContexts: List<VariantContext>, outputFileName: String,
                             refGenomeSequence: Map<String,NucSeq>, altHeaderLines:Set<VCFHeaderLine>) {
    val writer = VariantContextWriterBuilder()
        .unsetOption(Options.INDEX_ON_THE_FLY)
        .setOutputFile(File(outputFileName))
        .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
        .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
        .build()

    val header = createGenericHeader(listOf(sampleName),altHeaderLines)
    addSequenceDictionary(header, refGenomeSequence)
    writer.writeHeader(header)
    for(variant in variantContexts) {
        writer.add(variant)
    }

    writer.close()
}

fun createGenericHeader(taxaNames: List<String>, altLines:Set<VCFHeaderLine>): VCFHeader {
    val headerLines = createGenericHeaderLineSet() as MutableSet<VCFHeaderLine>
    headerLines.addAll(altLines)
    return VCFHeader(headerLines, taxaNames)
}


fun createGenericHeaderLineSet(): Set<VCFHeaderLine> {
    val headerLines: MutableSet<VCFHeaderLine> = HashSet()
    headerLines.add(VCFFormatHeaderLine("AD", 3, VCFHeaderLineType.Integer, "Allelic depths for the ref and alt alleles in the order listed"))
    headerLines.add(
        VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Read Depth (only filtered reads used for calling)")
    )
    headerLines.add(VCFFormatHeaderLine("GQ", 1, VCFHeaderLineType.Integer, "Genotype Quality"))
    headerLines.add(VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"))
    headerLines.add(
        VCFFormatHeaderLine("PL", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification")
    )
    headerLines.add(VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Total Depth"))
    headerLines.add(VCFInfoHeaderLine("NS", 1, VCFHeaderLineType.Integer, "Number of Samples With Data"))
    headerLines.add(VCFInfoHeaderLine("AF", 3, VCFHeaderLineType.Integer, "Allele Frequency"))
    headerLines.add(VCFInfoHeaderLine("END", 1, VCFHeaderLineType.Integer, "Stop position of the interval"))
    // These last 4 are needed for assembly g/hvcfs, but not for reference.  I am keeping them in as header
    // lines but they will only be added to the data lines if they are present in the VariantContext.
    headerLines.add(VCFInfoHeaderLine("ASM_Chr", 1, VCFHeaderLineType.String, "Assembly chromosome"))
    headerLines.add(VCFInfoHeaderLine("ASM_Start", 1, VCFHeaderLineType.Integer, "Assembly start position"))
    headerLines.add(VCFInfoHeaderLine("ASM_End", 1, VCFHeaderLineType.Integer, "Assembly end position"))
    headerLines.add(VCFInfoHeaderLine("ASM_Strand", 1, VCFHeaderLineType.String, "Assembly strand"))
    return headerLines
}

/**
 * function to bgzip and create tabix index on a gvcf file
 *
 * This method uses the -f option on both bgzip and tabix to overwrite any existing files
 * This handles the case where there already exists a bgzipped file, however
 * it is still required to have the non-bzipped file present.
 * @return the filename of the bgzipped gvcf
 */
fun bgzipAndIndexGVCFfile(gvcfFileName: String): String {

    try {
        // First bgzip the file

        // bgzip the file - needed to create index
        // use the -f option to overwrite any existing file
        myLogger.info("bgzipping  file ${gvcfFileName}")
        val gvcfGzippedFile = gvcfFileName + ".gz"
        var builder = ProcessBuilder("conda","run","-n","phgv2-conda",
            "bgzip", "-f", gvcfFileName)

        var process = builder.start()
        var error: Int = process.waitFor()
        if (error != 0) {
            myLogger.warn("\nERROR $error creating bgzipped  version of file: $gvcfFileName")
            throw IllegalStateException("bgzipAndIndexGVCFfile: error trying to bgzip file ${gvcfFileName}: ${error}")
        }

        // File has been gzipped, now index it.
        // Use the -f option to overwrite any existing index
        // We will use bcftools to create the csi index
        // ORiginal PHG used tabix, we wnat csi indexes to allow for large genomes e.g wheat.
        // TileDB supports .csi indexed files.
        builder = ProcessBuilder("conda","run","-n","phgv2-conda","bcftools", "index", "-c",gvcfGzippedFile)
        process = builder.start()
        error = process.waitFor()
        if (error != 0) {
            myLogger.warn("\nERROR $error creating tabix indexed  version of file: $gvcfGzippedFile")
            throw IllegalStateException("bgzipAndIndexGVCFfile: error trying to run bcftools index -c on file ${gvcfGzippedFile}: ${error}")
        }
        return gvcfGzippedFile
    } catch (exc:Exception) {
        throw IllegalStateException("bgzipAndIndexGVCFfile: error bgzipping and/or indexing file ${gvcfFileName}: Message:${exc.message}")
    }

}

/**
 * Function to add a sequence Dictionary based on the reference genome.  This uses the loaded genome to get the lengths.
 */
fun addSequenceDictionary(vcfheader : VCFHeader, refGenomeSequence: Map<String,NucSeq>) {

    val sequenceRecordList = refGenomeSequence.keys.map { SAMSequenceRecord(it,
        refGenomeSequence[it]!!.size()) }

    vcfheader.setSequenceDictionary(SAMSequenceDictionary(sequenceRecordList))
}

// Method to verify anchors for genome intervals have no overlapping positions
// Overlapping intervals are not supported in the PHG
fun verifyIntervalRanges(intervalFile: String): Set<String> {
    val overlappingPositions: MutableSet<String> = HashSet() // overlaps to be returned
    val intervalRanges: RangeSet<Position> = TreeRangeSet.create()
    // Read the anchor file, store to RangeSet, check for overlaps as you add
    // Store overlapping anchors to a Set to be returned to calling method
    myLogger.info("verifyIntervalRanges: checking file ${intervalFile} for overlaps")
    try {
        bufferedReader(intervalFile).use { br ->
            var curLine = br.readLine()
            while (curLine != null) {
                if (curLine.uppercase(Locale.getDefault()).contains("CHROMSTART")) continue
                val tokens =
                    //curLine.split("\\t".toRegex()).dropLastWhile { it.isEmpty() }.toTypedArray()
                    curLine.split("\t").dropLastWhile{ it.isEmpty() }.toTypedArray()
                val chrom = tokens[0]
                val interval =
                    Range.closedOpen(
                        Position(
                            chrom,
                            tokens[1].toInt()
                        ), Position(chrom, tokens[2].toInt())
                    )
                if (intervalRanges.intersects(interval)) {
                    overlappingPositions.add(curLine)
                }
                intervalRanges.add(interval)
                curLine = br.readLine()
            }
        }
    } catch (exc: Exception) {
        throw IllegalArgumentException("verifyIntervalRanges :  error reading intervals file: " + exc.message)
    }
    return overlappingPositions
}

/**
 * Helper method to create a Reference Range VariantContext for assemblies.  The DP value
 * is defaulted to 0 for assemblies.  If this is not set, -1 is used as default in
 * GenotypeBuilder.   That causes assembly problems down the line when storing the
 * value as a byte in a long.
 *
 * The AD (ref depth) is consistent with what BioKotlin uses for refDepth in MAFToGVCF
 * The PL (Phred-scaled likelihoods) are consistent with what BioKotlin uses for refPL in MAFToGVCF
 * The Depth is set to 30 to match BioKotlin's MAFToGVCF totalDP variable
 *
 * @param refSequence
 * @param assemblyTaxon
 * @param refRangeStart
 * @param refRangeEnd
 * @param asmStart
 * @param asmEnd
 * @return
 */
fun createRefRangeVC(refSequence: Map<String,NucSeq>, assemblyTaxon: String, refRangeStart: Position, refRangeEnd: Position,
    asmStart: Position?, asmEnd: Position?, asmStrand:String): VariantContext {
    // -1 added to refRangeStart.position to get the correct base in call below
    val firstRefAllele = Allele.create(refSequence[refRangeStart.contig]!!.get(refRangeStart.position-1).toString(),true)
    val gt = GenotypeBuilder().name(assemblyTaxon).alleles(Arrays.asList(firstRefAllele)).DP(30).AD(intArrayOf(30, 0)).PL(intArrayOf(0,90,90)).make()
    check(refRangeStart.position <= refRangeEnd.position) { "createRefRangeVC - start position greater than end: start=" +
                refRangeStart.position + " end=" + refRangeEnd.position
    }
    val vcb = VariantContextBuilder()
        .chr(refRangeStart.contig)
        .start(refRangeStart.position.toLong())
        .stop(refRangeEnd.position.toLong())
        .attribute("END", refRangeEnd.position)
        .alleles(Arrays.asList(firstRefAllele, Allele.NON_REF_ALLELE))
        .genotypes(gt)

    // Add assembly coordinates as attributes
    // If the asmStart and asmEnd position are null, nothing will be added.
    val vcbWithASMAnnos = addASMCoordsToVariantContextBuilder(vcb, asmStart, asmEnd,asmStrand)
    return vcbWithASMAnnos.make()
}

/**
 * Helper method to create a SNP Variant context for assemblies.  The DP value
 * is defaulted to 0 for assemblies.  If this is not set, -1 is used as default in
 * GenotypeBuilder.   That causes assembly problems down the line when storing the
 * value as a byte in a long.
 * @param assemblyTaxon
 * @param startPosition
 * @param endPosition
 * @param calls
 * @param asmStart
 * @param asmEnd
 * @return
 */
fun createSNPVC(assemblyTaxon: String, startPosition: Position, endPosition: Position, calls: Pair<String, String>,
                asmStart: Position?, asmEnd: Position?,asmStrand: String?): VariantContext {
    val refCall = Allele.create(calls.first, true)
    val altCall = Allele.create(calls.second, false)
    check(startPosition.position <= endPosition.position) {
        "createSNPVC - start postion greater than end: start=" +
                startPosition.position + " end=" + endPosition.position
    }

    //Need to add AD for Alt >0 here so that the API will work correctly.  Otherwise it is treated as missing as it thinks AD = 0,0.
    // When coming from an assembly it should always use the ALT in a SNP pos
    val gt = GenotypeBuilder().name(assemblyTaxon).alleles(Arrays.asList(altCall)).DP(30).AD(intArrayOf(0, 30, 0)).make()
    val vcb = VariantContextBuilder()
        .chr(startPosition.contig)
        .start(startPosition.position.toLong())
        .stop(endPosition.position.toLong())
        .alleles(Arrays.asList(refCall, altCall, Allele.NON_REF_ALLELE))
        .genotypes(gt)

    // Add assembly coordinates as attributes
    val vcbWithASMAnnos = addASMCoordsToVariantContextBuilder(vcb, asmStart, asmEnd,asmStrand)
    return vcbWithASMAnnos.make()
}

// Creates a VariantContext record for a Haplotype vcf file
fun createHVCFRecord(assemblyTaxon: String, startPosition: Position, endPosition: Position, calls: Pair<String, String>): VariantContext {
    val refCall = Allele.create(calls.first, true)
    val altCall = Allele.create(symbolicAllele(calls.second), false)
    check(startPosition.position <= endPosition.position) {
        "createHVCFRecord - start position greater than end: start=" +
                startPosition.position + " end=" + endPosition.position
    }

    //Need to add AD for Alt >0 here so that the API will work correctly.  Otherwise it is treated as missing as it thinks AD = 0,0.
    // When coming from an assembly it should always use the ALT in a SNP pos
    val gt = GenotypeBuilder().name(assemblyTaxon).alleles(Arrays.asList(altCall)).DP(2).AD(intArrayOf(0, 2, 0)).make()
    val vcb = VariantContextBuilder()
        .chr(startPosition.contig)
        .start(startPosition.position.toLong())
        .stop(endPosition.position.toLong())
        .attribute("END", endPosition.position)
        .alleles(Arrays.asList(refCall, altCall)) // no NON_REF allele for h.vcf files
        .genotypes(gt)


    return vcb.make()
}

fun createDiploidHVCFRecord(sampleName: String, startPosition: Position, endPosition: Position, calls: List<String?>, refAlleleStr: String): VariantContext {
    val refCall = Allele.create(refAlleleStr, true)
    val sampleAlleles = calls.map { if (it == null) Allele.NO_CALL else Allele.create(symbolicAllele(it!!), false) }

    val distinctAltCalls = sampleAlleles.filter{ !it.isNoCall }.distinct()


    check(startPosition.position <= endPosition.position) {
        "createDiploidHVCFRecord: start position greater than end for ${startPosition.contig}: ${startPosition.position} - ${endPosition.position}"
    }

    val alleleList = mutableListOf(refCall)
    alleleList.addAll(distinctAltCalls)
    val gt = GenotypeBuilder().name(sampleName).alleles(sampleAlleles).phased(true).make()
    val vcb = VariantContextBuilder()
        .chr(startPosition.contig)
        .start(startPosition.position.toLong())
        .stop(endPosition.position.toLong())
        .attribute("END", endPosition.position)
        .alleles(alleleList)
        .genotypes(gt)


    return vcb.make()
}

// Symbolic alleles for VariantContext records must be surrounded
// by <> characters.  This method adds them.
fun symbolicAllele(allele: String): String {
    return "<${allele}>"
}
fun addASMCoordsToVariantContextBuilder(vcb: VariantContextBuilder, asmStart: Position?, asmEnd: Position?, asmStrand:String?) :VariantContextBuilder {
    var newVCB = vcb
    //Only add in the Assembly start and end if they exist
    if (asmStart != null && asmEnd != null) {
        // Set the asm coordinates as VC record attributes
        check(asmStart.contig == asmEnd.contig) { "createRefRangeVC - assembly start and end contigs do not match: start=" +
                asmStart.contig + " end=" + asmEnd.contig
        }
        newVCB = vcb.attribute("ASM_Chr", asmStart.contig)
            .attribute("ASM_Start", asmStart.position)
            .attribute("ASM_End", asmEnd.position)
            .attribute("ASM_Strand", asmStrand)
    }
    return newVCB
}

// Create an Md5 checksum.  While this method will create a check sum from
// any number of methods, for PHG we expect it to always be Md5
fun getChecksumForString(seq: String, protocol: String="Md5"): String {
    // from https://www.mkyong.com/java/java-md5-hashing-example/
    try {
        val md = MessageDigest.getInstance(protocol)
        md.update(seq.toByteArray())
        val byteData = md.digest()
        // convert the byte to hex format
        val sb = StringBuffer()
        for (idx in byteData.indices) {
            sb.append(Integer.toString((byteData[idx].toInt() and 0xff) + 0x100, 16).substring(1))
        }
        return sb.toString()
    } catch (exc: NoSuchAlgorithmException) {
        myLogger.warn("getChecksumForString: problem getting checksum: " + exc.message)
        throw IllegalStateException("CheckSum: getChecksumForString: error: " + exc.message)
    }
}

/**
 * This function verifies the path given for the tiledb datasets is valid
 *  If it is not, an exception is thrown and the user is directed to create
 *  the dataset using the Initdb command.
 *  The uri should either be gvcf_dataset or hvcf_dataset
 *  The user determines the parent folder name where these datasets live
 *  The actual tiledb dataset names are constant and are either gvcf_dataset or hvcf_dataset
 */

fun verifyURI(dbPath:String,uri:String,condaEnvPrefix:String): Boolean {
    // Check that the user supplied db folder exists
    check(File(dbPath).exists()) { "Folder $dbPath does not exist - please send a valid path that indicates the parent folder for your tiledb datasets." }

    // Check if the dataset exists
    val dataset = "${dbPath}/${uri}"
    val datasetPath = Paths.get(dataset)

    if (File(dataset).exists() && Files.isRegularFile(datasetPath)) {
        throw IllegalArgumentException("URI ${dataset}is a file, not a tiledb dataset folder.  The parent folder must not contain any files/folders named gvcf_dataset or hvcf_dataset that is not a tiledb created URI")
    }

    // Create tne temp folder if it doesn't exist
    // This will be used to write output files from ProcessBuilder commands
    // called elsewhere in this class
    val tempDir = "${dbPath}/log"
    Files.createDirectories(Paths.get(tempDir))

    if (File(dataset).exists()  && Files.isDirectory(Paths.get(dataset))){
        // check if is a tiledb dataset
        var command = if (condaEnvPrefix.isNotBlank()) mutableListOf("conda","run","-p",condaEnvPrefix,"tiledbvcf","stat","--uri",dataset)
        else mutableListOf("conda","run","-n","phgv2-conda","tiledbvcf","stat","--uri",dataset)

        var builder = ProcessBuilder(command)
        var redirectOutput = tempDir + "/tiledb_statURI_output.log"
        var redirectError = tempDir + "/tiledb_statURI_error.log"
        builder.redirectOutput( File(redirectOutput))
        builder.redirectError( File(redirectError))

        // verify if the output.log contains "Version"
        // if not, then the URI is not a tiledbvcf URI
        myLogger.info("begin Command:" + builder.command().joinToString(" "))

        try {
            var process = builder.start()
            var error = process.waitFor()
            if (error != 0) {
                myLogger.error("LoadTiledbH tiledb stat returned error code $error")
                throw IllegalArgumentException("Error: URI is not a tiledb URI folder created via the tiledb create command: ${error}")
            }

        } catch (exc: java.lang.Exception) {
            myLogger.error("Error: could not run tiledb stat on ${uri}.")
            throw IllegalArgumentException("Error running ProcessBuilder to stat tiledb URI: ${exc}")
        }

        myLogger.info("Using  TileDB datasets created in folder $dbPath.")
        return true
    } else {
        myLogger.info("TileDB datasets not found in folder $dbPath. Either send a valid dbPath folder variable, or run InitDB to create the datasets in the specified folder.")
        throw IllegalArgumentException("TileDB datasets not found in folder $dbPath. Please run InitDB to create the datasets.")
    }
}