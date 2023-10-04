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
import java.io.File
import java.security.MessageDigest
import java.security.NoSuchAlgorithmException
import java.util.*
import java.util.logging.Logger

private val myLogger = Logger.getLogger("net.maizegenetics.pangenome.db_loading.VariantLoadingUtils")

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
        var builder = ProcessBuilder(
            "bgzip", "-f", gvcfFileName)

        var process = builder.start()
        var error: Int = process.waitFor()
        if (error != 0) {
            myLogger.warning("\nERROR $error creating bgzipped  version of file: $gvcfFileName")
            throw IllegalStateException("bgzipAndIndexGVCFfile: error trying to bgzip file ${gvcfFileName}: ${error}")
        }

        // File has been gzipped, now index it.
        // Use the -f option to overwrite any existing index
        // We will use bcftools to create the csi index
        // ORiginal PHG used tabix, we wnat csi indexes to allow for large genomes e.g wheat.
        // TileDB supports .csi indexed files.
        val tabixFile = gvcfGzippedFile + ".csi"
        builder = ProcessBuilder("bcftools", "index", "-c",gvcfGzippedFile)
        process = builder.start()
        error = process.waitFor()
        if (error != 0) {
            myLogger.warning("\nERROR $error creating tabix indexed  version of file: $gvcfGzippedFile")
            throw IllegalStateException("bgzipAndIndexGVCFfile: error trying to run bcftools index -c on file ${gvcfGzippedFile}: ${error}")
        }
        return gvcfGzippedFile
    } catch (exc:Exception) {
        throw IllegalStateException("bgzipAndIndexGVCFfile: error bgzipping and/or indexing file ${gvcfFileName}")
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
    println("verifyIntervalRanges: checking file ${intervalFile} for overlaps")
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
 * @param refSequence
 * @param assemblyTaxon
 * @param refRangeStart
 * @param refRangeEnd
 * @param asmStart
 * @param asmEnd
 * @return
 */
fun createRefRangeVC(refSequence: Map<String,NucSeq>, assemblyTaxon: String, refRangeStart: Position, refRangeEnd: Position,
    asmStart: Position, asmEnd: Position): VariantContext {
    val firstRefAllele = Allele.create(refSequence[refRangeStart.contig]!!.get(0).toString())
    val gt = GenotypeBuilder().name(assemblyTaxon).alleles(Arrays.asList(firstRefAllele)).DP(2).AD(intArrayOf(2, 0)).make()
    check(refRangeStart.position <= refRangeEnd.position) { "createRefRangeVC - start position greater than end: start=" +
                refRangeStart.position + " end=" + refRangeEnd.position
    }
    var vcb = VariantContextBuilder()
        .chr(refRangeStart.contig)
        .start(refRangeStart.position.toLong())
        .stop(refRangeEnd.position.toLong())
        .attribute("END", refRangeEnd.position)
        .alleles(Arrays.asList(firstRefAllele, Allele.NON_REF_ALLELE))
        .genotypes(gt)

    // Add assembly coordinates as attributes
    if (asmStart != null && asmEnd != null) {
        // Set the asm coordinates as VC record attributes
        vcb = vcb.attribute("ASM_Start", asmStart.position)
        vcb = vcb.attribute("ASM_End", asmEnd.position)
    }
    return vcb.make()
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
                asmStart: Position, asmEnd: Position): VariantContext {
    val refCall = Allele.create(calls.first, true)
    val altCall = Allele.create(calls.second, false)
    check(startPosition.position <= endPosition.position) {
        "createSNPVC - start postion greater than end: start=" +
                startPosition.position + " end=" + endPosition.position
    }

    //Need to add AD for Alt >0 here so that the API will work correctly.  Otherwise it is treated as missing as it thinks AD = 0,0.
    // When coming from an assembly it should always use the ALT in a SNP pos
    val gt = GenotypeBuilder().name(assemblyTaxon).alleles(Arrays.asList(altCall)).DP(2).AD(intArrayOf(0, 2, 0)).make()
    var vcb = VariantContextBuilder()
        .chr(startPosition.contig)
        .start(startPosition.position.toLong())
        .stop(endPosition.position.toLong())
        .alleles(Arrays.asList(refCall, altCall, Allele.NON_REF_ALLELE))
        .genotypes(gt)


    //Only add in the Assembly start and end if they exist
    if (asmStart != null && asmEnd != null) {
        // Set the asm coordinates as VC record attributes
        vcb = vcb.attribute("ASM_Start", asmStart.position)
        vcb = vcb.attribute("ASM_End", asmEnd.position)
    }
    return vcb.make()
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
        println("getChecksumForString: problem getting checksum: " + exc.message)
        throw IllegalStateException("CheckSum: getChecksumForString: error: " + exc.message)
    }
}