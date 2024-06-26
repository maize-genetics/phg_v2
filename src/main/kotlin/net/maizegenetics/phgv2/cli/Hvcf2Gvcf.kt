package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.VariantContextComparator
import htsjdk.variant.vcf.VCFEncoder
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeader
import htsjdk.variant.vcf.VCFHeaderLine
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.exportVariantContext
import net.maizegenetics.phgv2.utils.parseALTHeader
import net.maizegenetics.phgv2.utils.verifyURI
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * This class takes a phgv2 hvcf file that was created via path-finding, and
 * turns it into a gvcf file.  It assumes the gvcf file for all samples
 * reflected in the hvcf file live in the tiledb database.  These will be
 * exported for use.
 *
 * The file name for the newly created gvcf file will be the same as the hvcf file
 * (minus extension).
 *
 * Steps:
 * 1.  verify initial data
 * 2.  read hvcf file, get sample names represented
 *    - sample names come from the ALT header lines - grab all the SampleName fields
 *    and add to a set.
 * 3.  export the gvcf files from the tiledb database
 * 3.  bgzip and index the gvcf files
 * 4.  ALso from the hvcf files, create a list of ReferenceRange objects for each sample
 * 5.  Read the gvcf files, a sample at a time, and for each sample, pull the gvcf records
 * that overlap the ReferenceRanges for that sample.  This will be a list of VariantContext records.
 * 6.  Merge the list of vcRecords for each sample together, sorted by reference ranges, and write
 * to a new gvcf file in the output directory.
 * 7.  The header for the gvcf file will be the same as the header from one of the gvcf files.
 *     With the exception that sampleName will be changed to the sample name from the hvcf file.
 * 8.  The gvcf file will be written to the output directory.
 *
 *
 */
class Hvcf2Gvcf: CliktCommand(help = "Create  h.vcf files from existing PHG created g.vcf files")  {
    private val myLogger = LogManager.getLogger(Hvcf2Gvcf::class.java)

    val hvcfDir by option("--hvcf-dir", help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--hvcf-dir must not be blank"
            }
        }

    val condaEnvPrefix by option (help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
        .default("")

    val outputDir by option (help = "Output directory for the gVCF files.  If not provided, the current working directory is used.")
        .default("")

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        .default("")
    val referenceFile by option(help = "Path to local Reference FASTA file needed for sequence dictionary")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--reference-file must not be blank"
            }
        }

    override fun run() {
        //TODO("Not yet implemented")
        val dbPath = if (dbPath.isBlank()) {
            System.getProperty("user.dir")
        } else {
            dbPath
        }

        println("Hvcf2Gvcf: dbPath = $dbPath, hvcfDir = $hvcfDir, outputDir = $outputDir, condaEnvPrefix = $condaEnvPrefix")

        // Verify the tiledbURI - verifyURI will throw an exception if the URI is not valid
        val validDB = verifyURI(dbPath,"hvcf_dataset",condaEnvPrefix)

        buildGvcfFromHvcf(dbPath, referenceFile, outputDir, hvcfDir, condaEnvPrefix)
    }

    fun buildGvcfFromHvcf(dbPath: String, referenceFile: String, outputDir: String, hvcfDir: String, condaEnvPrefix: String) {
        // load the reference file
        println("LCJ - in buildGvcfFromHvcf")
        val refSeq = CreateMafVcf().buildRefGenomeSeq(referenceFile)
        // get list of hvcf files
        // walk the gvcf directory process files with g.vcf.gz extension
        File(hvcfDir).walk().filter { !it.isHidden && !it.isDirectory }
            .filter { it.name.endsWith("h.vcf.gz")  || it.name.endsWith("h.vcf")  }.toList()
            .forEach { hvcfFile ->
                println("buildGvcfFromHvcf: Processing hvcf file: ${hvcfFile.name}")
                val records = processSingleHVCF(outputDir, hvcfFile, dbPath,condaEnvPrefix)

                val sample = hvcfFile.toString().substringAfterLast("/").substringBefore(".")
                val gvcfFile = "$outputDir/${sample}.g.vcf"
                exportVariantContext(sample,records,gvcfFile, refSeq,setOf())
            }

    }

    fun processSingleHVCF(outputDir:String, hvcfFile: File, dbPath: String, condaEnvPrefix: String): List<VariantContext> {

        val reader = VCFFileReader(hvcfFile,false)
        // We also need to print to the new gvcf all the headers from 1 of the accessed gvcf files
        val header = reader.fileHeader
        val altHeaders = parseALTHeader(header = header)
        val sampleNames = altHeaders.values.map { it.sampleName() }.toSet()
        val contigNames = header.contigLines.map { it.id } // contig names needed for sorting

        // This checks for existing gvcf files, and if they don't exist, exports them
        val exportSuccess = exportGvcfFiles(sampleNames, outputDir, dbPath, condaEnvPrefix)

        // Using the hvcfdFileReader, walk the hvcf file and for each entry create
        // a ReferenceRange object and add that to a list of Reference Range objects
        // for the sample to which it applies.  You will get the sample by indexing the
        // altHeaders map with the key of the current record's ID.
        val sampleToRefRanges = mutableMapOf<String, MutableList<ReferenceRange>>()

        // Get the sampleName from the hvcf file.  THis is the name
        // without the path, and without the extension of h.vcf, hvcf, h.vcf.gz or hvcf.gz
        val pathSample = hvcfFile.toString().substringAfterLast("/").substringBefore(".")

        // process the hvcf records into the sampleToRefRanges map
        reader.forEach { context ->
            val id = context.getAlternateAllele(0).displayString.removeSurrounding("<", ">")
            val altHeader = altHeaders[id]
            val sampleName = altHeader?.sampleName() ?: throw IllegalStateException("No ALT header found for record: ${context.id}")
            val refRange = ReferenceRange(context.contig, context.start, context.end)
            val sampleList = sampleToRefRanges.getOrDefault(sampleName, mutableListOf())
            sampleList.add(refRange)
            sampleToRefRanges[sampleName] = sampleList

        }
        reader.close()

        // Now we can read the gvcfs, a sample at a time, into the gvcfRecords list and process.
        // or we could process these in parallel - how much memory do we need to hold all?  it isn't
        // the genome, just the gvcf but that could be large.

        // For each sample, read the gvcf file and pull the gvcf records that overlap the ReferenceRanges
        // for that sample.  This will be a list of VariantContext records.  After we have the list of vcRecords
        // for each sample, we will merge them together, sorted by reference ranges, and write to a new gvcf file
        // in the output directory.

        // THis may not be the correct version, either.  Do I need a map of REferenceRange to VariantContext records?
        // in which case, do I lose the sample names?

        //val refRangeToVariantContext = mutableMapOf<ReferenceRange, List<VariantContext>>()
        val refRangeToVariantContext = mutableMapOf<ReferenceRange, MutableList<VariantContext>>()

        sampleToRefRanges.forEach { sample, ranges ->
            val gvcfFile = "$outputDir/${sample}.vcf" // tiledb wrote with extension .vcf

            val gvcfReader = VCFFileReader(File(gvcfFile),false)

            // This needs to loop through both the List of ReferenceRanges and the gvcfRecords
            // It should find entries in the gvcf file whose positions overlap those of the reference ranges
            // and add them to the gvcfRecords list.
            val gvcfVariants  = mutableListOf<VariantContext>()
            gvcfReader.use { reader ->
                for (vc in reader) {
                    gvcfVariants.add(vc)
                }
            }
            gvcfReader.close()

            //val rangeToGvcfRecords = findOverlappingRecords(ranges, gvcfVariants)
            val rangeToGvcfRecords = findOverlappingRecordsForSample(ranges, gvcfVariants)

            //Add the rangeToGvcfRecords to the refRangeToVariantContext map
            // we can use "plus" but it creates a new map containing the combined entries
            // and assigns it back to refRangeToVariantContext
            // chatGPT says that below is more efficient than using "plus"
            // refRangeToVariantContext.plus(rangeToGvcfRecords)
            // refRangeToVariantContext.addAll(rangeToGvcfRecords) fails as "addAll" is not an option of a mutableMap

            for ((range, variantList) in rangeToGvcfRecords) {
                val existingList = refRangeToVariantContext.getOrPut(range) { mutableListOf() }
                existingList.addAll(variantList)
            }
        }
        // Put all the VariantContext records from the refRangeToVariantContext map
        // onto a list, then sort that list.
        val allRecords = mutableListOf<VariantContext>()
        refRangeToVariantContext.values.forEach { allRecords.addAll(it) }
        val variants = allRecords.sortedWith(VariantContextComparator(contigNames))

        return variants
    }

    // function to export the gvcf file if they don't exist
    fun exportGvcfFiles(sampleNames:Set<String>, outputDir:String, dbPath:String, condaEnvPrefix:String):Boolean {
        val success = true
        // Check if the files listed in gvcfFIles already exist
        val gvcfFiles = mutableListOf<String>()
        sampleNames.forEach { sampleName ->
            val gvcfFile = "$outputDir/${sampleName}.vcf"
            gvcfFiles.add(gvcfFile)
        }

        println("exportGvcfFiles: gvcfFiles= ${gvcfFiles}")
        val missingFiles = gvcfFiles.filter { !File(it).exists() }
        // For the entries in the missingFiles list, create a list of sampleNames.
        // The sampleNames are the entry in the missingFiles list, remove up to and
        // including the first "/" and remove the ".vcf" extension
        val missingSampleNames = missingFiles.map { it.substringAfterLast("/").substringBeforeLast(".") }

        // Setup the conda enviroment portion of the command
        var command = if (condaEnvPrefix.isNotBlank()) mutableListOf("conda","run","-p",condaEnvPrefix) else mutableListOf("conda","run","-n","phgv2-conda")
        // Call ExportVcf with the outputDir and the missingSamleNames list
        // to create the gvcf files
        var dataCommand = mutableListOf(
            "conda",
            "run",
            "-n",
            "phgv2-conda",
            "tiledbvcf",
            "export",
            "--uri",
            "$dbPath/gvcf_dataset",
            "-O",
            "v",
            "--sample-names",
            missingSampleNames.joinToString(","),
            "--output-dir",
            outputDir
        )

        // join the conda and data portion of the commands
        command.addAll(dataCommand)
        val builder = ProcessBuilder(command)

        val redirectError = "$outputDir/export_gvcf_error.log"
        val redirectOutput = "$outputDir/export_gvcf_output.log"
        builder.redirectOutput(File(redirectOutput))
        builder.redirectError(File(redirectError))

        myLogger.info("ExportVcf Command: " + builder.command().joinToString(" "))
        println("TILEDB ExportVcf Command: " + builder.command().joinToString(" "))
        val process = builder.start()
        val error = process.waitFor()
        if (error != 0) {
            myLogger.error("tiledbvcf export for: $missingSampleNames run via ProcessBuilder returned error code $error")
            throw IllegalStateException("Error running tiledbvcf export of dataset $dbPath/gvcf_dataset for: $missingSampleNames. error: $error")
        }

        // For all files in missingSampleNames, rename them to g.vcf
        // tiledb writes the exported files as .vcf
//        missingSampleNames.forEach { sample ->
//            File("$outputDir/$sample.vcf").renameTo(File("$outputDir/${sample}.g.vcf"))
//        }

        return success
    }

    fun findOverlappingRecordsForSample(ranges:List<ReferenceRange>, variantContexts:List<VariantContext>):MutableMap<ReferenceRange,MutableList<VariantContext>> {
        // split ranges and variantContexts by chromosome
        val overlappingRecords = mutableMapOf<ReferenceRange, MutableList<VariantContext>>()
        val rangesByChrom = ranges.groupBy { it.contig }
        val variantContextsByChrom = variantContexts.groupBy { it.contig }
        // call findOverlappingRecords for each chromosome, add the results to overlappingRecords
        rangesByChrom.keys.forEach { chrom ->
            val chromRanges = rangesByChrom[chrom] ?: emptyList()
            val chromVariants = variantContextsByChrom[chrom] ?: emptyList()
            val chromOverlappingRecords = findOverlappingRecords(chromRanges, chromVariants)
            overlappingRecords.putAll(chromOverlappingRecords)
        }
        return overlappingRecords
    }

    // This is based on CreateMafVCF:convertGVCFToHVCFForChrom() that determines if a variant or part
    // of a variant is in a reference range.  It differs in that we are not creating hvcf meta data,
    // but rather, at this stage, are merely finding the overlapping variants
    // Another function will be called to split the gvcf records if they extend  beyond the reference range
    // either at the beginning or the end.
    fun findOverlappingRecords(ranges:List<ReferenceRange>, variantContexts:List<VariantContext>):MutableMap<ReferenceRange,MutableList<VariantContext>> {
        val refRangeToVariantContext = mutableMapOf<ReferenceRange, MutableList<VariantContext>>() // this will be returned
        var currentVariantIdx = 0
        println("LCJ - findOverlappingRecords, size of ranges: ${ranges.size}, size of variantContexts: ${variantContexts.size}")
        for (range in ranges) {
            val regionStart = range.start
            val regionEnd = range.end
            val regionChrom = range.contig
            val tempVariants = mutableListOf<VariantContext>()

            while (currentVariantIdx < variantContexts.size) {
                val currentVariant = variantContexts[currentVariantIdx]

                // check different cases for the variant
                //If variant is fully contained in Bed region add to temp list and increment currentVariantIdx
                //If variant is partially contained in Bed region add to temp list do not increment as we need to see if the next bed also overlaps
                //If variant is not contained in Bed region, skip and do not increment as we need to see if the next bed overlaps
                if(CreateMafVcf().bedRegionContainedInVariant(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)), currentVariant)) {
                    // THis is the case where the region is completely contained within the variant,
                    // meaning the variant may overlap the region.  We need to adjust the asm positions
                    val fixedVariants = fixPositions(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)), listOf(currentVariant))
                    // Add the fixedVariants to the refRangeToVariantContext map, for the current range
                    val currentVariants = refRangeToVariantContext.getOrPut(range, { mutableListOf() })
                    currentVariants.addAll(fixedVariants)
                    refRangeToVariantContext[range] = currentVariants
                    tempVariants.clear()
                    break
                }
                if(CreateMafVcf().variantFullyContained(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)), currentVariant)) {
                    //This is the case where the variant is completely contained within the region
                    tempVariants.add(currentVariant)
                    currentVariantIdx++
                }
                else if(CreateMafVcf().variantPartiallyContainedStart(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)),currentVariant)) {
                    tempVariants.add(currentVariant)
                    break
                }
                else if(CreateMafVcf().variantPartiallyContainedEnd(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)), currentVariant)) {
                    tempVariants.add(currentVariant)
                    currentVariantIdx++
                }
                else if(CreateMafVcf().variantAfterRegion(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)), currentVariant)) {
                    //write out what is in tempVariants
                    if(tempVariants.isNotEmpty()) {
                        val fixedVariants = fixPositions(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)), tempVariants)
                        val currentVariants = refRangeToVariantContext.getOrPut(range, { mutableListOf() })
                        currentVariants.addAll(fixedVariants)
                        refRangeToVariantContext[range] = currentVariants
                        tempVariants.clear()
                    }
                    //move up Bed region
                    break
                }
                else { //this is the case where the Variant is behind the BED region
                    //move up Variant
                    currentVariantIdx++
                }
            }
            // Process the last variants in the list
            if(tempVariants.isNotEmpty()) {
                val fixedVariants = fixPositions(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)), tempVariants)
                val currentVariants = refRangeToVariantContext.getOrPut(range, { mutableListOf() })
                currentVariants.addAll(fixedVariants)
                refRangeToVariantContext[range] = currentVariants
                tempVariants.clear()
            }
        }

        return refRangeToVariantContext
    }

    // This fixes both the ref start/end and the asm positions in the tempVariants list.  It should return an
    // ammended version of the list. Code is based on CreateMafVcf:convertGVCFRecordsToHVCFMetaData()
    // but only deals with the ASM positions portion of this code.  In addition, it adds code
    // to adjust the variant's start/end positions. Indels are not currently handled.
    fun fixPositions(region: Pair<Position,Position>, variants: List<VariantContext> ): List<VariantContext> {
        val fixedVariants = mutableListOf<VariantContext>()

        // TODO ADD INDEL SUPPORT
        //Take the first and the last variantContext
        val firstVariant = variants.first()
        val lastVariant = variants.last()

        //val check strandedness of the variants
        val firstStrand = firstVariant.getAttributeAsString("ASM_Strand","+")

        val lastStrand = lastVariant.getAttributeAsString("ASM_Strand","+")
        //Resize the first and last variantContext ref/ASM start and end based on the regions
        var newStartPositions = resizeVCandASMpositions(firstVariant, region.first.position, firstStrand)
        if(newStartPositions == Pair(-1,-1)) {
            newStartPositions = if(firstStrand == "+") Pair(firstVariant.start,firstVariant.getAttributeAsInt("ASM_Start",region.first.position))
            else Pair(firstVariant.end,firstVariant.getAttributeAsInt("ASM_End",region.first.position))
        }

        var newEndPositions = resizeVCandASMpositions(lastVariant, region.second.position, lastStrand)
        if(newEndPositions == Pair(-1,-1)) {
            newEndPositions = if(lastStrand == "+") Pair(lastVariant.end,lastVariant.getAttributeAsInt("ASM_End",region.second.position))
            else Pair(lastVariant.start,lastVariant.getAttributeAsInt("ASM_Start",region.second.position))
        }

        // At this point, we have changes for the first and last regions.  If the list size
        // is only 1, we change it based on newStartPositions and newEndPositions
        // If there are multiple, we change the first and last entries.  The first entry gets its start changed,
        // THe last entry gets  end values changed.  The middle entries are not changed.
        if (variants.size == 1) {
            val updatedFirstVariant = VariantContextBuilder(firstVariant)
                .start(newStartPositions.first.toLong())
                .stop(newEndPositions.first.toLong())
                .attribute("ASM_Start", newStartPositions.second)
                .attribute("ASM_End", newEndPositions.second)
                .make()
            fixedVariants.add(updatedFirstVariant)
        }
        else {
            // update the first and last variants, leaving those
            // in the middle unchanged
            val updatedFirstVariant = VariantContextBuilder(firstVariant)
                .start(newStartPositions.first.toLong())
                .attribute("ASM_Start", newStartPositions.second)
                .make()
            fixedVariants.add(updatedFirstVariant)
            val updatedLastVariant = VariantContextBuilder(lastVariant)
                .stop(newEndPositions.first.toLong())
                .attribute("ASM_End", newEndPositions.second)
                .make()

            if (variants.size > 2) {
                fixedVariants.addAll(variants.subList(1,variants.size-1))
            }
            fixedVariants.add(updatedLastVariant)
        }

        return fixedVariants

    }


    // The position is the position in the reference range.  The strand is the strand of the variant
    fun resizeVCandASMpositions(variant: VariantContext, position: Int, strand : String) : Pair<Int,Int> {
        //check to see if the variant is either a RefBlock or is a SNP with equal lengths
        var refAsmPos = Pair<Int,Int>(-1,-1) // this is what is returned
        return if (CreateMafVcf().isVariantResizable(variant)) {
            // if the position is < the start of the variant, then we return <position,asm_start>
            // if the postiion is > variant end, then we return <position, asm_end>
            // But the "position" must be modified in the above with an offset.
            // these 2 checks verify the position is within the variant range
            // if the strand is +, then we return position + offset
            // if the strand is -, then we return position - offset
            when {
                position < variant.start -> {
                    // The reference position starts before the variant, so we should
                    // keep it as the variant start
                    Pair(variant.start,variant.getAttributeAsInt("ASM_Start",variant.start))
                }
                position > variant.end -> {
                    // The reference position is after the variant, so we should
                    // keep it as the variant end
                    Pair(variant.end,variant.getAttributeAsInt("ASM_End",variant.end))
                }
                strand == "+" -> {
                    val offset = position - variant.start
                    Pair(position+offset,variant.getAttributeAsInt("ASM_Start",variant.start) + offset)
                }
                strand == "-" -> {
                    val offset = position - variant.start
                    Pair(position-offset,variant.getAttributeAsInt("ASM_Start",variant.end) - offset)
                }
                else -> refAsmPos
            }

        }
        else {
            refAsmPos
        }

    }

}