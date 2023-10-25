package net.maizegenetics.phgv2.cli

import biokotlin.genome.*
import biokotlin.seq.NucSeq
import biokotlin.seqIO.NucSeqIO
import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFAltHeaderLine
import htsjdk.variant.vcf.VCFHeaderLine
import htsjdk.variant.vcf.VCFHeaderVersion
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.createHVCFRecord
import net.maizegenetics.phgv2.utils.exportVariantContext
import net.maizegenetics.phgv2.utils.getChecksumForString
import java.io.File

class CreateMafVcf : CliktCommand(help = "Create gVCF and hVCF from Anchorwave MAF files") {

    val bed by option(help = "BED file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--bed must not be blank"
            }
        }

    val reference by option(help = "Reference FASTA file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--reference must not be blank"
            }
        }

    val mafDir by option(help = "MAF file directory")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--maf-dir must not be blank"
            }
        }

    val outputDir by option("-o", "--output-dir", help = "Name for output VCF file Directory")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--output-dir/-o must not be blank"
            }
        }

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--db-path must not be blank"
            }
        }

    /**
     * Function to create the ASM hVCF and gVCF.
     * It will first use Biokotlin to build the gVCF and then will use the BED file to extract out the hVCF information.
     */
    fun createASMHvcfs(dbPath: String, bedFileName: String, referenceFileName: String, mafDirName: String, outputDirName: String) {
        //load the bed file into some data structure
//        val ranges = bedfileToSRangeSet(bedFileName,referenceFileName)
        val ranges = loadRanges(bedFileName)
        val refGenomeSequence = buildRefGenomeSeq(referenceFileName)

        //loop through the maf files in mafDirName and getGVCFVariantsFromMafFile
        File(mafDirName).walk().filter { !it.isHidden && !it.isDirectory }
            .filter { it.extension == "maf" }
            .forEach {
                val sampleName = it.nameWithoutExtension //This will likely need to change in the future
                val gvcfVariants = getGVCFVariantsFromMafFile(referenceFileName, it.absolutePath, it.nameWithoutExtension)

                exportVariantContext(sampleName, gvcfVariants, "${outputDirName}/${it.nameWithoutExtension}.g.vcf",refGenomeSequence, setOf())


                val asmHeaderLines = mutableSetOf<VCFHeaderLine>()
                //convert the GVCF records into hvcf records
                val hvcfVariants = convertGVCFToHVCF(sampleName, ranges, gvcfVariants, refGenomeSequence, dbPath, asmHeaderLines)
                //export the hvcfRecords
                exportVariantContext(sampleName, hvcfVariants, "${outputDirName}/${it.nameWithoutExtension}.h.vcf",refGenomeSequence, asmHeaderLines)
            }

    }

    fun loadRanges(bedFileName: String) : List<Pair<Position, Position>> {
        return bufferedReader(bedFileName).readLines().map { line ->
            val lineSplit = line.split("\t")
            val chrom = lineSplit[0]
            val start = lineSplit[1].toInt()+1
            val end = lineSplit[2].toInt()
            Pair(Position(chrom,start),Position(chrom,end))
        }.sortedBy { it.first }
    }

    //Function to load in the reference using Biokotlin
    fun buildRefGenomeSeq(referenceFileName: String) : Map<String, NucSeq> {
        return NucSeqIO(referenceFileName).readAll()
    }

    /**
     * Simple Wrapper Function to call biokotlin's MAFToGVCF.getVariantContextsfromMAF
     */
    fun getGVCFVariantsFromMafFile(referenceFileName: String, mafFileName : String, sampleName: String, fillGaps: Boolean = false, sortMaf: Boolean = true) : List<VariantContext> {
        return MAFToGVCF().getVariantContextsfromMAF(
            mafFileName,
            referenceFileName,
            sampleName,
            fillGaps,
            sortMaf
        )
    }

//    fun convertGVCFToHVCF(bedRanges : SRangeSet, gvcfVariants: List<VariantContext>) : List<VariantContext> {
    fun convertGVCFToHVCF(sampleName: String, bedRanges : List<Pair<Position,Position>>, gvcfVariants: List<VariantContext>,
                          refGenomeSequence : Map<String, NucSeq>, agcArchiveName: String, asmHeaders: MutableSet<VCFHeaderLine>) : List<VariantContext> {
        // group the gvcfVariants by contig
        val gvcfVariantsByContig = gvcfVariants.groupBy { it.contig }

        val bedRegionsByContig = bedRanges.groupBy { it.first.contig }

        return gvcfVariantsByContig.keys.filter { bedRegionsByContig.containsKey(it) }.flatMap { convertGVCFToHVCFForChrom(sampleName, bedRegionsByContig[it]!!, refGenomeSequence, agcArchiveName, gvcfVariantsByContig[it]!!, asmHeaders) }
    }

    fun convertGVCFToHVCFForChrom(sampleName: String, bedRanges: List<Pair<Position,Position>>, refGenomeSequence: Map<String, NucSeq>, agcArchiveName: String, variantContexts: List<VariantContext>, asmHeaders: MutableSet<VCFHeaderLine> ) : List<VariantContext> {

        /**
         * Loop through the bed file
         * Loop through the gvcf records as well
         *
         * We need to determine if our BED region overlaps with the gvcf record
         * To do this we need to collect the gvcf records into their corresponding bed Regions
         * Then from those collected regions, we take the first and last ones and resize based on the bed regions to get the asm_Start and asm_End
         * Then extract the sequence out of the AGC archive and md5 hash it
         * Then call the createHVCFRecord with this information
         */
        val outputVariants = mutableListOf<VariantContext>()
        var currentVariantIdx = 0
        for(region in bedRanges) {
            val regionStart = region.first.position
            val regionEnd = region.second.position
            val regionChrom = region.first.contig
            val tempVariants = mutableListOf<VariantContext>()

            check(regionChrom in refGenomeSequence.keys) { "Chromosome $regionChrom not found in reference" }

            //Need to subtract here as the Biokotlin NucSeq is 0 based
            val refRangeSeq = refGenomeSequence[regionChrom]!![regionStart-1..regionEnd-1]

            while (currentVariantIdx < variantContexts.size) {
                val currentVariant = variantContexts[currentVariantIdx]

                //check different cases for the variant
                //If variant is fully contained in Bed region add to temp list and increment currentVariantIdx
                //If variant is partially contained in Bed region add to temp list do not increment as we need to see if the next bed also overlaps
                //If variant is not contained in Bed region, skip and do not increment as we need to see if the next bed overlaps
                if(bedRegionContainedInVariant(region, currentVariant)) {
                    outputVariants.add(convertGVCFRecordsToHVCF(sampleName, region, refRangeSeq, agcArchiveName, listOf(currentVariant), asmHeaders))
                    tempVariants.clear()
                    break
                }
                if(variantFullyContained(region, currentVariant)) {
                    //This is the case where the variant is completely contained within the region
                    tempVariants.add(currentVariant)
                    currentVariantIdx++
                }
                else if(variantPartiallyContainedStart(region,currentVariant)) {
                    tempVariants.add(currentVariant)
                    break
                }
                else if(variantPartiallyContainedEnd(region, currentVariant)) {
                    tempVariants.add(currentVariant)
                    currentVariantIdx++
                }
                else if(variantAfterRegion(region, currentVariant)) {
                    //write out what is in tempVariants
                    if(tempVariants.isNotEmpty()) {
                        outputVariants.add(
                            convertGVCFRecordsToHVCF(
                                sampleName,
                                region,
                                refRangeSeq,
                                agcArchiveName,
                                tempVariants,
                                asmHeaders
                            )
                        )
                        tempVariants.clear()

                    }
                    //move up Bed region
                    break
                }
                else { //this is the case if the Variant is behind the BED region
                    //move up Variant
                    currentVariantIdx++
                }
            }
            if(tempVariants.isNotEmpty()) {
                outputVariants.add(
                    convertGVCFRecordsToHVCF(
                        sampleName,
                        region,
                        refRangeSeq,
                        agcArchiveName,
                        tempVariants,
                        asmHeaders
                    )
                )
                tempVariants.clear()

            }
        }

        return outputVariants
    }

    fun bedRegionContainedInVariant(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
        return variant.contig == region.first.contig && variant.start <= region.first.position && variant.end >= region.second.position
    }
    fun variantFullyContained(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
        return variant.contig == region.first.contig && variant.start >= region.first.position && variant.end <= region.second.position
    }

    fun variantPartiallyContainedStart(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
        return variant.contig == region.first.contig &&
                variant.start >= region.first.position &&
                variant.start <= region.second.position &&
                variant.end > region.second.position
    }

    fun variantPartiallyContainedEnd(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
        return variant.contig == region.first.contig &&
                variant.end <= region.second.position &&
                variant.end >= region.first.position &&
                variant.start < region.first.position
    }
    fun variantAfterRegion(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
        return variant.contig == region.first.contig && variant.start > region.second.position
    }

    fun convertGVCFRecordsToHVCF(sampleName: String, region: Pair<Position,Position>, refRangeSeq: NucSeq, agcArchiveName: String ,variants: List<VariantContext>, asmHeaders: MutableSet<VCFHeaderLine> ) : VariantContext {
        //Take the first and the last variantContext
        val firstVariant = variants.first()
        val lastVariant = variants.last()

        //val check strandedness of the variants
        val strand = firstVariant.getAttributeAsString("ASM_Strand","+")
        check(strand == lastVariant.getAttributeAsString("ASM_Strand","+")) { "Strand of first and last variantContexts do not match" }
        //Resize the first and last variantContext ASM start and end based on the regions
        var newASMStart = resizeVariantContext(firstVariant, region.first.position, strand)
        if(newASMStart == -1) {
            newASMStart = if(strand == "+") firstVariant.getAttributeAsInt("ASM_Start",region.first.position)
                            else firstVariant.getAttributeAsInt("ASM_End",region.first.position)
        }

        var newASMEnd = resizeVariantContext(lastVariant, region.second.position, strand)
        if(newASMEnd == -1) {
            newASMEnd = if(strand == "+") lastVariant.getAttributeAsInt("ASM_End",region.second.position)
                            else lastVariant.getAttributeAsInt("ASM_Start",region.second.position)
        }

        //Extract out the sequence for the assembly haplotype
        val assemblyHaplotypeSeq:String = extractSequenceFromHaplotypes(firstVariant.sampleNames.first(), firstVariant.getAttributeAsString("ASM_Chr",""), newASMStart, newASMEnd, agcArchiveName)
        //md5 hash the assembly sequence
        val assemblyHaplotypeHash = getChecksumForString(assemblyHaplotypeSeq)
        //md5 has the refSequence
        val refSeqHash = getChecksumForString(refRangeSeq.toString())

        //create the asmHeader lines
        asmHeaders.add(
            VCFAltHeaderLine(
                "<ID=${assemblyHaplotypeHash}, Description=\"haplotype data for line: ${sampleName}\">,Number=9,Source=\"${agcArchiveName}\",Contig=\"${region.first.contig}\",Start=\"${region.first.position}\",End=\"${region.second.position}\",Asm_Contig=\"${firstVariant.getAttributeAsString("ASM_Chr","")}\",Asm_Start=\"${newASMStart}\",Asm_End=\"${newASMEnd}\",Checksum=\"Md5\",RefRange=\"${refSeqHash}\">",
                VCFHeaderVersion.VCF4_2
            )
        )

        //build a variant context of the HVCF with the hashes
        return createHVCFRecord(sampleName, region.first, region.second, Pair(refRangeSeq[0].toString(), assemblyHaplotypeHash))
    }

    fun extractSequenceFromHaplotypes(sampleName:String, chrom: String, start: Int, end: Int, agcArchiveName: String) : String {
        //will be replaced by Lynn's call using sampleName, chrom, start, end
        val assemblySeqs = buildRefGenomeSeq(agcArchiveName)
        val assemblySeq = assemblySeqs[chrom]!!
        return assemblySeq[start-1..end-1].seq() //Need to subtract 1 from both boundaries as NucSeq is 0 based
    }



    /**
     * Function will return -1 if unable to resize the variantcontext due to its type(mostly an INDEL)
     * If the requested position is outside of the current variants coordinates it will return the ASM_Start for + strand and ASM_End for - strand
     */
    fun resizeVariantContext(variant: VariantContext, position: Int, strand : String) : Int {
        //check to see if the variant is either a RefBlock or is a SNP with equal lengths
        return if(isVariantResizable(variant)) {
            when {
                position < variant.start && strand == "+" -> variant.getAttributeAsInt("ASM_Start",variant.start)
                position < variant.start && strand == "-" -> variant.getAttributeAsInt("ASM_End",variant.end)
                position > variant.end && strand == "+" -> variant.getAttributeAsInt("ASM_End",variant.end)
                position > variant.end && strand == "-" -> variant.getAttributeAsInt("ASM_Start",variant.start)
                strand == "+" -> {
                    val offset = position - variant.start
                    variant.getAttributeAsInt("ASM_Start",variant.start) + offset
                }
                strand == "-" -> {
                    val offset = position - variant.start
                    variant.getAttributeAsInt("ASM_End",variant.end) - offset
                }
                else -> -1
            }
        } else {
            -1
        }
    }

    fun isVariantResizable(variant: VariantContext) : Boolean {
        return when {
            variant.getReference().baseString.length == 1 && variant.end - variant.start > 0 -> true //refBlock
            variant.reference.baseString.length == variant.getAlternateAllele(0).baseString.length -> true //This covers both SNPs and multiallelic polymorphisms
            else -> false
        }
    }

    override fun run() {
        createASMHvcfs(dbPath, bed, reference, mafDir, outputDir)
    }

}