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
import net.maizegenetics.phgv2.utils.*
import java.io.File

data class HVCFRecordMetadata(val sampleName: String, val refSeq : String = "", val asmSeq : String = "",
                              val refContig : String, val refStart: Int, val refEnd: Int,
                              val asmRegions: List<Pair<Position,Position>>)
data class DisplayRegion(val contig: String, val start: Int, val end: Int)


class CreateMafVcf : CliktCommand(help = "Create g.vcf and h.vcf files from AnchorWave MAF files") {

    val bed by option(help = "BED file with entries that define the haplotype boundaries")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--bed must not be blank"
            }
        }

    val referenceFile by option(help = "Path to local Reference FASTA file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--reference-file must not be blank"
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
     * if [twoGvcfs] is true, then the output will be split into two gvcf files, one for each gamete.
     */
    fun createASMHvcfs(dbPath: String, bedFileName: String, referenceFileName: String, mafDirName: String, outputDirName: String, twoGvcfs:Boolean=false) {
        //load the bed file into some data structure
//        val ranges = bedfileToSRangeSet(bedFileName,referenceFileName)
        val ranges = loadRanges(bedFileName)
        println("CreateASMHvcfs: calling buildRefGenomeSeq")
        val refGenomeSequence = buildRefGenomeSeq(referenceFileName)

        //loop through the maf files in mafDirName and getGVCFVariantsFromMafFile
        File(mafDirName).walk().filter { !it.isHidden && !it.isDirectory }
            .filter { it.extension == "maf" }
            .forEach {
                println("CreateASMHvcfs: processing ${it.absolutePath}")
                val sampleName = it.nameWithoutExtension //This will likely need to change in the future
                val gvcfVariants = getGVCFVariantsFromMafFile(refGenomeSequence, it.absolutePath, it.nameWithoutExtension, twoGvcfs=twoGvcfs)
                //export the gvcfRecords
                if (gvcfVariants.size == 1){
                    println("createASMHvcfs: gvcfVariants.size == 1")
                    val sampleName = gvcfVariants.keys.first()
                    val variants = gvcfVariants.values.first()
                    println("createASMHvcfs: processing sampleName = $sampleName")
                    exportVariantContext(sampleName, variants, "${outputDirName}/${it.nameWithoutExtension}.g.vcf",refGenomeSequence, setOf())
                    bgzipAndIndexGVCFfile("${outputDirName}/${it.nameWithoutExtension}.g.vcf")

                    val asmHeaderLines = mutableMapOf<String,VCFHeaderLine>()
                    //convert the GVCF records into hvcf records
                    println("createASMHvcfs: calling convertGVCFToHVCF for $sampleName")
                    val hvcfVariants = convertGVCFToHVCF(dbPath,sampleName, ranges, variants, refGenomeSequence, dbPath, asmHeaderLines)
                    val asmHeaderSet = asmHeaderLines.values.toSet()
                    //export the hvcfRecords
                    println("createASMHvcfs: calling exportVariantContext for $sampleName")
                    exportVariantContext(sampleName, hvcfVariants, "${outputDirName}/${it.nameWithoutExtension}.h.vcf",refGenomeSequence, asmHeaderSet)
                    //bgzip the files
                    bgzipAndIndexGVCFfile("${outputDirName}/${it.nameWithoutExtension}.h.vcf")
                } else if (gvcfVariants.size == 2) {
                    println("createASMHvcfs: gvcfVariants.size == 2")
                    val gvcfOutput = "${outputDirName}/${it.nameWithoutExtension}.g.vcf"
                    val outputNames = MAFToGVCF().twoOutputFiles(gvcfOutput)
                    gvcfVariants.entries.forEachIndexed { index, (name, variants) ->
                        val outputFile =
                            exportVariantContext(name, variants, outputNames[index], refGenomeSequence, setOf())
                        bgzipAndIndexGVCFfile(outputNames[index])
                        val asmHeaderLines = mutableMapOf<String,VCFHeaderLine>()
                        //convert the GVCF records into hvcf records
                        val hvcfVariants = convertGVCFToHVCF(dbPath,sampleName, ranges, variants, refGenomeSequence, dbPath, asmHeaderLines)
                        val asmHeaderSet = asmHeaderLines.values.toSet()
                        //export the hvcfRecords
                        exportVariantContext(sampleName, hvcfVariants, "${outputDirName}/${it.nameWithoutExtension}.h.vcf",refGenomeSequence, asmHeaderSet)
                        //bgzip the files
                        bgzipAndIndexGVCFfile("${outputDirName}/${it.nameWithoutExtension}.h.vcf")
                    }

                }

            }

    }

    /**
     * Simple function to load a BED file in.  This will be replaced by a lightweight Biokotlin ranges class eventually.
     *
     * This will sort in alphabetical order first then will check if there are numbers in the chromosome name and will
     * sort those numerically. This means that chr10 will come after chr2.
     */
    fun loadRanges(bedFileName: String) : List<Pair<Position, Position>> {
        return bufferedReader(bedFileName).readLines().map { line ->
            val lineSplit = line.split("\t")
            val chrom = lineSplit[0]
            val start = lineSplit[1].toInt()+1
            val end = lineSplit[2].toInt()
            Pair(Position(chrom,start),Position(chrom,end))
        }.sortedWith(compareBy(SeqRangeSort.alphaThenNumberSort) { positionRange:Pair<Position,Position> -> positionRange.first.contig})
    }

    //Function to load in the reference using Biokotlin
    fun buildRefGenomeSeq(referenceFileName: String) : Map<String, NucSeq> {
        return NucSeqIO(referenceFileName).readAll()
    }

    /**
     * Simple Wrapper Function to call biokotlin's MAFToGVCF.getVariantContextsfromMAF
     * The parameters to BioKotlin's (version 0.10) getVariantContextsfromMAF are:
     * 1. [mafFileName] - the name of the MAF file
     * 2. [refGenomeSequence] - the reference genome sequence
     * 3. [sampleName] - the name of the sample
     * 4. [fillGaps] - whether to fill in the gaps in the MAF file, default is false
     * 5. [twoGvcfs] - whether to split the output into two gvcf files, default false
     * 6. [outJustGT] - whether to only output the GT field, default false
     * 7. [outputType] - either gvcf or vcf (BioKotlin doesn't currently have h.vcf, default is gvcf)
     *
     */
    fun getGVCFVariantsFromMafFile(refSeq: Map<String,NucSeq>, mafFileName : String, sampleName: String, fillGaps: Boolean = false, twoGvcfs: Boolean = false) : Map<String,List<VariantContext>> {
        return MAFToGVCF().getVariantContextsfromMAF(
            mafFileName,
            refSeq,
            sampleName,
            fillGaps,
            twoGvcfs
        )
    }

    /**
     * Function to convert a GVCF file into an HCVF file
     */
    fun convertGVCFToHVCF(dbPath: String,sampleName: String, bedRanges : List<Pair<Position,Position>>, gvcfVariants: List<VariantContext>,
                          refGenomeSequence : Map<String, NucSeq>, agcArchiveName: String, asmHeaders: MutableMap<String,VCFHeaderLine>) : List<VariantContext> {
        // group the gvcfVariants by contig
        val gvcfVariantsByContig = gvcfVariants.groupBy { it.contig }

        val bedRegionsByContig = bedRanges.groupBy { it.first.contig }


        println("in convertGVCFToHVCF: sort and call converGVCFToHVCFForChrom")
        return gvcfVariantsByContig.keys
            .sortedWith(compareBy(SeqRangeSort.alphaThenNumberSort){ name:String -> name}) //Need to do a sort here as we need to make sure we process the chromosomes in
            .filter { bedRegionsByContig.containsKey(it) }
            .flatMap { convertGVCFToHVCFForChrom(dbPath, sampleName, bedRegionsByContig[it]!!, refGenomeSequence, agcArchiveName, gvcfVariantsByContig[it]!!, asmHeaders) }
    }

    fun convertGVCFToHVCFForChrom(dbPath: String, sampleName: String, bedRanges: List<Pair<Position,Position>>, refGenomeSequence: Map<String, NucSeq>, agcArchiveName: String, variantContexts: List<VariantContext>, asmHeaders: MutableMap<String,VCFHeaderLine> ) : List<VariantContext> {
        
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
        println("in convertGVCFToHVCFForChrom: bedRanges.size = ${bedRanges.size}")
        val outputVariantMetadata = mutableListOf<HVCFRecordMetadata>()
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
                    outputVariantMetadata.add(
                        convertGVCFRecordsToHVCFMetaData(
                            sampleName,
                            region,
                            refRangeSeq,
                            listOf(currentVariant)
                        )
                    )
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
                        outputVariantMetadata.add(
                            convertGVCFRecordsToHVCFMetaData(sampleName,
                                region,
                                refRangeSeq,
                                tempVariants
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
                outputVariantMetadata.add(convertGVCFRecordsToHVCFMetaData(
                    sampleName,
                    region,
                    refRangeSeq,
                    tempVariants
                ))
                tempVariants.clear()
            }
        }

        val metaDataWithSequence = addSequencesToMetaData(dbPath, outputVariantMetadata)
        val outputVariants = convertMetaDataToHVCFContexts(metaDataWithSequence, asmHeaders, dbPath)

        return outputVariants
    }

    /**
     * Function to see if the BED region is fully contained within a VariantContext
     * Indels are left-justified
     * Bed:       |---|
     * Var: |--------------|
     */
    fun bedRegionContainedInVariant(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
        val end = if (variant.type == VariantContext.Type.SYMBOLIC) variant.end else variant.start
        return variant.contig == region.first.contig && variant.start <= region.first.position && end >= region.second.position
    }

    /**
     * Function to see if the VariantContext is fully contained within a BED region
     * Indels are left-justified
     * Bed: |--------------|
     * Var:       |---|
     */
    fun variantFullyContained(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
        val end = if (variant.type == VariantContext.Type.SYMBOLIC) variant.end else variant.start
        return variant.contig == region.first.contig && variant.start >= region.first.position && end <= region.second.position
    }

    /**
     * Function to see if the start of the variant is partially contained in the BED region
     * Indels are left-justified
     * Bed: |--------------|
     * Var:              |---|
     */
    fun variantPartiallyContainedStart(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
        val end = if (variant.type == VariantContext.Type.SYMBOLIC) variant.end else variant.start
        return variant.contig == region.first.contig &&
                variant.start >= region.first.position &&
                variant.start <= region.second.position &&
                end > region.second.position
    }

    /**
     * Function to see if the end of the variant is partially contained in the BED region
     * Indels are left-justified
     * Bed:      |--------------|
     * Var:    |---|
     */
    fun variantPartiallyContainedEnd(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
        val end = if (variant.type == VariantContext.Type.SYMBOLIC) variant.end else variant.start
        return variant.contig == region.first.contig &&
                end <= region.second.position &&
                end >= region.first.position &&
                variant.start < region.first.position
    }

    /**
     * Function to see if the variant is after the bed region
     * Bed: |--------------|
     * Var:                   |---|
     */
    fun variantAfterRegion(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
        return variant.contig == region.first.contig && variant.start > region.second.position
    }


    /**
     * Function to extract all the needed information out of the ASM gVCF record and put them into HVCFRecordMetadata objects
     * This will first try to resize the positions based on the ref start position and then will extract out all the other information.
     */
    fun convertGVCFRecordsToHVCFMetaData(sampleName: String, region: Pair<Position,Position>, refRangeSeq: NucSeq, variants: List<VariantContext> ) : HVCFRecordMetadata {
        //Take the first and the last variantContext
        val firstVariant = variants.first()
        val lastVariant = variants.last()

        //val check strandedness of the variants
        val firstStrand = firstVariant.getAttributeAsString("ASM_Strand","+")
        //LCJ - the line below should be removed?  This was the first exception from scale-testing Zack was fixing,
        // ie it is ok if these strands are not the same -handling inversion
        // check(strand == lastVariant.getAttributeAsString("ASM_Strand","+")) { "Strand of first and last variantContexts do not match" }

        val lastStrand = lastVariant.getAttributeAsString("ASM_Strand","+")
        //Resize the first and last variantContext ASM start and end based on the regions
        var newASMStart = resizeVariantContext(firstVariant, region.first.position, firstStrand)
        if(newASMStart == -1) {
            newASMStart = if(firstStrand == "+") firstVariant.getAttributeAsInt("ASM_Start",region.first.position)
                else firstVariant.getAttributeAsInt("ASM_End",region.first.position)
        }

        var newASMEnd = resizeVariantContext(lastVariant, region.second.position, lastStrand)
        if(newASMEnd == -1) {
            newASMEnd = if(lastStrand == "+") lastVariant.getAttributeAsInt("ASM_End",region.second.position)
                else lastVariant.getAttributeAsInt("ASM_Start",region.second.position)
        }

        val regions = buildNewAssemblyRegions(newASMStart,newASMEnd,variants)


        return HVCFRecordMetadata(sampleName=sampleName, refSeq = refRangeSeq.toString(), asmSeq = "",
            refContig = region.first.contig, refStart = region.first.position, refEnd = region.second.position,
            regions)

    }

    /**
     * Function to build the new assembly region coordinates based on the new Start and end and the list of VariantContexts
     * Any consecutive regions should be merged together so we do not make the eventual string too long
     * The output will be a List<Pair<Position,Position>> which will be the new coordinates for the all the  assembly regions
     */
    fun buildNewAssemblyRegions(newStart: Int, newEnd: Int, variants: List<VariantContext>) : List<Pair<Position,Position>> {
        val variantsConverted = variants.map { convertVariantContextToPositionRange(it) }

        //resize the first and last position based on the strand
        val resizedFirst = resizePositionRange(variantsConverted.first(),newStart,true)
        val resizedLast = resizePositionRange(variantsConverted.last(),newEnd,false)

        //merge the first and last with the rest of the variants
        val mergedVariants = mutableListOf<Pair<Position,Position>>()
        if(variantsConverted.size == 1) {
            val resizedFirstAndLast = resizePositionRange(resizePositionRange(variantsConverted.first(), newStart, true), newEnd, false)
            mergedVariants.add(resizedFirstAndLast)
        } else {
            mergedVariants.add(resizedFirst) // add the first variant
            if (variantsConverted.size > 2) { // add any variants in the middle
                mergedVariants.addAll(variantsConverted.subList(1, variantsConverted.size - 1))
            }
            mergedVariants.add(resizedLast) // add the last variant
        }

        //merge the consecutive regions
        val mergedConsecutiveVariants = mergeConsecutiveRegions(mergedVariants)

        return mergedConsecutiveVariants
    }

    /**
     * Strand aware function to merge together consecutive assembly regions.  This is done to reduce the number of entries in the hvcf alt header.
     */
    fun mergeConsecutiveRegions(variants: List<Pair<Position,Position>>) : List<Pair<Position,Position>> {
        val mergedConsecutiveVariants = mutableListOf<Pair<Position,Position>>()
        var currentStart = variants.first().first
        var currentEnd = variants.first().second
        for(i in 1 until variants.size) {
            val nextStart = variants[i].first
            val nextEnd = variants[i].second

            //Check to see if the next region is only 1 bp long.  If so we need to check both normal and inverted boundaries and extend if it makes sense, if not reset the currentStart and currentEnd
            if(nextStart.position == nextEnd.position) {
                //Check to see if the next region is on the + strand
                //Using nextStart here as it equals nextEnd
                if(nextStart.position == currentEnd.position + 1 || nextStart.position == currentEnd.position -1) {
                    currentEnd = nextStart
                }
                else {
                    mergedConsecutiveVariants.add(Pair(currentStart,currentEnd))
                    currentStart = nextStart
                    currentEnd = nextEnd
                }
            }
            else if(currentEnd < currentStart && nextEnd < nextStart) {
                //This is the case where we have a variant that is on the - strand
                //We need to check if the next variant is consecutive
                if(nextStart.position == (currentEnd.position - 1)) {
                    currentEnd = nextEnd
                }
                else {
                    mergedConsecutiveVariants.add(Pair(currentStart,currentEnd))
                    currentStart = nextStart
                    currentEnd = nextEnd
                }
            }
            else if(nextStart.position == (currentEnd.position + 1)) {
                currentEnd = nextEnd
            }
            else {
                mergedConsecutiveVariants.add(Pair(currentStart,currentEnd))
                currentStart = nextStart
                currentEnd = nextEnd
            }
        }
        mergedConsecutiveVariants.add(Pair(currentStart,currentEnd))
        return mergedConsecutiveVariants
    }

    /**
     * Function to convert a VariantContext into a Pair<Position,Position> which will be the assembly starts and ends of the variantContext
     */
    fun convertVariantContextToPositionRange(variant: VariantContext) : Pair<Position,Position> {
        //get out the assembly coords
        val contig = variant.getAttributeAsString("ASM_Chr","")
        val start = variant.getAttributeAsInt("ASM_Start",variant.start)
        val end = variant.getAttributeAsInt("ASM_End",variant.end)
        return Pair(Position(contig, start), Position(contig, end))
    }

    /**
     * Function to resize the Position range based on the new position.  If isFirst is true then it will resize the start position, otherwise it will resize the end position
     */
    fun resizePositionRange(positionRange: Pair<Position,Position>, newPosition : Int, isFirst: Boolean) : Pair<Position,Position> {
        return if(isFirst) {
            //Slide the start position to the new position
            Pair(Position(positionRange.first.contig,newPosition),positionRange.second)
        } else {
            //Slide the end position to the new position
            Pair(positionRange.first,Position(positionRange.second.contig,newPosition))
        }
    }


    /**
     * This function will bulk load sequences in from the AGC record and then will associate the returned sequences
     * with the metadata record which contains the coordiantes for the query and will add in the asmSeq.
     */
    fun addSequencesToMetaData(dbPath: String, metadata: List<HVCFRecordMetadata>) : List<HVCFRecordMetadata> {
        //get out the assembly coordinates and build them into the regions
        val metaDataToRangeLookup = metadata.map {
            val queries = mutableListOf<String>()
            val displayNames = mutableListOf<String>()

            for(range in it.asmRegions) {
                if(range.first.position-1 > range.second.position-1) {
                    queries.add("${range.first.contig}@${it.sampleName}:${range.second.position-1}-${range.first.position-1}")
                    displayNames.add("${range.first.contig}:${range.second.position-1}-${range.first.position-1}")
                }
                else {
                    queries.add("${range.first.contig}@${it.sampleName}:${range.first.position - 1}-${range.second.position - 1}")
                    displayNames.add("${range.first.contig}:${range.first.position-1}-${range.second.position-1}")
                }
            }

            Triple(it,queries, displayNames)
        }

        val ranges = metaDataToRangeLookup.flatMap { it.second }
        println("LCJ addSequencesToMetaData: calling retrieveAgcContigs with ranges.size = ${ranges.size}")
        val seqs = retrieveAgcContigs(dbPath,ranges)

        return metaDataToRangeLookup.map { it.first.copy(asmSeq = buildSeq(seqs,it.third,it.first)) } //This is a useful way to keep things immutable
    }


    /**
     * Function to build the haplotype sequence based on the list of display regions and the given haplotype sequence object.
     * The sequence is already extracted out of AGC and stored in the seqs map.
     * The seqs map is keyed by a Pair of (sampleName, displayRegion) and the value is the NucSeq object.
     * The Pair may look something like ("B97", "1:1-1000") or ("B97", "chr1")
     */
    fun buildSeq(seqs: Map<Pair<String,String>,NucSeq> ,displayRegions : List<String>, hvcfRecordMetadata: HVCFRecordMetadata) : String {
        val hapSeqRegions = hvcfRecordMetadata.asmRegions

        return displayRegions.mapIndexed{ idx, currentDisplayRegion ->
            val currentHapSeqRegion = hapSeqRegions[idx]

            val seq = seqs[Pair(hvcfRecordMetadata.sampleName,currentDisplayRegion)]!!

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
     * Simple function to convert the all the HVCFRecordMetadata records into VariantContext records
     */
    fun convertMetaDataToHVCFContexts(metaData: List<HVCFRecordMetadata>, asmHeaders: MutableMap<String,VCFHeaderLine>, dbPath:String): List<VariantContext> {
        return metaData.map { convertMetaDataRecordToHVCF(it, asmHeaders, dbPath) }
    }

    /**
     * Simple function to convert a single metadata record into a VariantContext.
     * This will also create the ALT tag and add it to the asmHeaders object for use during export.
     */
    fun convertMetaDataRecordToHVCF(metaDataRecord: HVCFRecordMetadata, asmHeaders: MutableMap<String, VCFHeaderLine>, dbPath: String): VariantContext {
        val assemblyHaplotypeSeq:String = metaDataRecord.asmSeq
        //md5 hash the assembly sequence
        val assemblyHaplotypeHash = getChecksumForString(assemblyHaplotypeSeq)
        check(metaDataRecord.refSeq.isNotEmpty()) { "Reference sequence is empty" }
        //md5 has the refSequence
        val refSeqHash = getChecksumForString(metaDataRecord.refSeq)

        //create the asmHeader lines
        if(!asmHeaders.containsKey(assemblyHaplotypeHash)) {
            asmHeaders[assemblyHaplotypeHash] =
            VCFAltHeaderLine(
                "<ID=${assemblyHaplotypeHash}, Description=\"haplotype data for line: ${metaDataRecord.sampleName}\">," +
                        "Source=\"${dbPath}/assemblies.agc\",SampleName=\"${metaDataRecord.sampleName}\"," +
                        "Regions=\"${metaDataRecord.asmRegions.map { "${it.first.contig}:${it.first.position}-${it.second.position}" }.joinToString(",")}\"," +
                        "Checksum=\"Md5\",RefRange=\"${refSeqHash}\">",
                VCFHeaderVersion.VCF4_2
            )
        }


        //build a variant context of the HVCF with the hashes
        return createHVCFRecord(metaDataRecord.sampleName, Position(metaDataRecord.refContig,metaDataRecord.refStart),
            Position(metaDataRecord.refContig, metaDataRecord.refEnd ),
            Pair(metaDataRecord.refSeq[0].toString(), assemblyHaplotypeHash))
    }


    /**
     * Function will return -1 if unable to resize the variantcontext due to its type(mostly an INDEL)
     * If the requested position is outside of the current variants coordinates it will return the ASM_Start for + strand and ASM_End for - strand
     */
    fun resizeVariantContext(variant: VariantContext, position: Int, strand : String) : Int {
        //check to see if the variant is either a RefBlock or is a SNP with equal lengths
        return if(isVariantResizable(variant)) {
            when {
                position < variant.start -> variant.getAttributeAsInt("ASM_Start",variant.start)
                position > variant.end -> variant.getAttributeAsInt("ASM_End",variant.end)
                strand == "+" -> {
                    val offset = position - variant.start
                    variant.getAttributeAsInt("ASM_Start",variant.start) + offset
                }
                strand == "-" -> {
                    val offset = position - variant.start
                    variant.getAttributeAsInt("ASM_Start",variant.end) - offset
                }
                else -> -1
            }
        } else {
            -1
        }
    }

    /**
     * Function to check if a variant is resizable.  Only RefBlocks and SNPs are resizable
     */
    fun isVariantResizable(variant: VariantContext) : Boolean {
        return when {
            variant.getReference().baseString.length == 1 && variant.end - variant.start > 0 && variant.type == VariantContext.Type.SYMBOLIC -> true //refBlock
            variant.reference.baseString.length == variant.getAlternateAllele(0).baseString.length -> true //This covers both SNPs and multiallelic polymorphisms
            else -> false
        }
    }

    override fun run() {
        createASMHvcfs(dbPath, bed, referenceFile, mafDir, outputDir)
    }

}