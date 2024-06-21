package net.maizegenetics.phgv2.cli

import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.parseALTHeader
import net.maizegenetics.phgv2.utils.verifyURI
import org.apache.logging.log4j.LogManager
import org.jetbrains.kotlinx.dataframe.io.SupportedFormatSample
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

    override fun run() {
        //TODO("Not yet implemented")
        val dbPath = if (dbPath.isBlank()) {
            System.getProperty("user.dir")
        } else {
            dbPath
        }

        // Verify the tiledbURI
        // If it doesn't an exception will be thrown
        val validDB = verifyURI(dbPath,"hvcf_dataset",condaEnvPrefix)

        buildGvcfFromHvcf(dbPath, outputDir, hvcfDir, condaEnvPrefix)
        //buildFastaFromHVCF(dbPath, output, fastaType, hvcfDir, hvcfFile,condaEnvPrefix)
    }

    fun buildGvcfFromHvcf(dbPath: String, outputDir: String, hvcfDir: String, condaEnvPrefix: String) {
        // get list of hvcf files
        // walk the gvcf directory process files with g.vcf.gz extension
        File(hvcfDir).walk().filter { !it.isHidden && !it.isDirectory }
            .filter { it.name.endsWith("g.vcf.gz")  || it.name.endsWith("g.vcf")  }.toList()
            .forEach { hvcfFile ->
                val hvcfFileReader = VCFFileReader(hvcfFile,false)
                val records = processSingleHVCF(hvcfFile, dbPath,condaEnvPrefix)

            }

    }

    fun processSingleHVCF(hvcfFile: File, dbPath: String, condaEnvPrefix: String) {
        val reader = VCFFileReader(hvcfFile,false)
        val header = reader.fileHeader
        val altHeaders = parseALTHeader(header)
        val sampleNames = altHeaders.values.map { it.sampleName() }.toSet()

        // THis checks for existing gvcf files, and if they don't exist, exports them
        val exportSuccess = exportGvcfFiles(sampleNames, outputDir, dbPath, condaEnvPrefix)


        // Using the hvcfdFileReader, walk the hvcf file and for each entry create
        // a ReferenceRange object and add that to a list of Reference Range objects
        // for the sample to which it applies.  You will get the sample by indexing the
        // altHeaders map with the key of the current record's ID.
        val sampleToRefRanges = mutableMapOf<String, MutableList<ReferenceRange>>()

        // process the hvcf records into the sampleToRefRanges map
        reader.forEach { context ->
            val altHeader = altHeaders[context.id]
            val sampleName = altHeader?.sampleName() ?: throw IllegalStateException("No ALT header found for record: ${context.id}")
            val refRange = ReferenceRange(context.contig, context.start, context.end)
            val sampleList = sampleToRefRanges.getOrDefault(sampleName, mutableListOf())
            sampleList.add(refRange)
            sampleToRefRanges[sampleName] = sampleList

        }

        // Now we can read the gvcfs, a sample at a time, into the gvcfRecords list and process.
        // or we could process these in parallel - how much memory do we need to hold all?  it isn't
        // the genome, just the gvcf but that could be large.

        // For each sample, read the gvcf file and pull the gvcf records that overlap the ReferenceRanges
        // for that sample.  This will be a list of VariantContext records.  After we have the list of vcRecords
        // for each sample, we will merge them together, sorted by reference ranges, and write to a new gvcf file
        // in the output directory.

        // THis may not be the correct version, either.  Do I need a map of REferenceRange to VariantContext records?
        // in which case, do I lose the sample names?

        val refRangeToVariantContext = mutableMapOf<ReferenceRange, MutableList<VariantContext>>()
        // This is WRONG, but might give me ideas on what to do
        sampleToRefRanges.forEach { sample, ranges ->
            val gvcfFile = "$outputDir/${sample}.g.vcf"

            val gvcfReader = VCFFileReader(File(gvcfFile),false)

            // THis needs to loop through both the List of ReferenceRanges and the gvcfRecords
            // It should find entries in the gvcf file whose positions overlap those of the reference ranges
            // and add them to the gvcfRecords list.
            val rangeToGvcfRecords = findOverlappingRecords(ranges, gvcfReader);

            // now need to split the gvcf records that extend beyond the reference range
            // Once this is done we can put this set of records into the refRangeToVariantContext map


            // We will not yet be writing anything to a file

        }


        return
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

        val missingFiles = gvcfFiles.filter { !File(it).exists() }
        //For the entries in the missingFiles list, create a list of sampleNames.
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

        // join the 2 conda and data portion of the commands
        command.addAll(dataCommand)
        val builder = ProcessBuilder(command)

        val redirectError = "$outputDir/export_gvcf_error.log"
        val redirectOutput = "$outputDir/export_gvcf_output.log"
        builder.redirectOutput(File(redirectOutput))
        builder.redirectError(File(redirectError))

        myLogger.info("ExportVcf Command: " + builder.command().joinToString(" "))
        val process = builder.start()
        val error = process.waitFor()
        if (error != 0) {
            myLogger.error("tiledbvcf export for: $missingSampleNames run via ProcessBuilder returned error code $error")
            throw IllegalStateException("Error running tiledbvcf export of dataset $dbPath/gvcf_dataset for: $missingSampleNames. error: $error")
        }

        // For all files in missingSampleNames, rename them to g.vcf
        missingSampleNames.forEach { sample ->
            File("$outputDir/$sample.vcf").renameTo(File("$outputDir/${sample}.g.vcf"))
        }

        return success
    }

    // This is based on CreateMafVCF:convertGVCFToHVCFForChrom() that determines if a variant or part
    // of a variant is in a reference range.  It differs in that we are not creating hvcf meta data,
    // but rather, at this stage, or merely finding the overlapping variants
    // Another function will be called to split the gvcf records if they extend  beyond the reference range
    // either at the beginning or the end.
    fun findOverlappingRecords(ranges:List<ReferenceRange>, reader:VCFFileReader):Map<ReferenceRange,List<VariantContext>> {
        val refRangeToVariantContext = mutableMapOf<ReferenceRange, MutableList<VariantContext>>() // this will be returned
        var currentVariant = reader.iterator().next()
        for (range in ranges) {
            val regionStart = range.start
            val regionEnd = range.end
            val regionChrom = range.contig
            val tempVariants = mutableListOf<VariantContext>()

            while (currentVariant != null) {


                //check different cases for the variant
                //If variant is fully contained in Bed region add to temp list and increment currentVariantIdx
                //If variant is partially contained in Bed region add to temp list do not increment as we need to see if the next bed also overlaps
                //If variant is not contained in Bed region, skip and do not increment as we need to see if the next bed overlaps
                if(CreateMafVcf().bedRegionContainedInVariant(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)), currentVariant)) {
                    // THis is the case where the region is completely contained within the variant,
                    // meaning the variant may overlap the region.  We need to adjust the asm positions
                    val fixedVariants = fixASMPositions(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)), listOf(currentVariant))
                    refRangeToVariantContext.getOrPut(range, { mutableListOf() }).addAll(fixedVariants)
                    tempVariants.clear()
                    break
                }
                if(CreateMafVcf().variantFullyContained(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)), currentVariant)) {
                    //This is the case where the variant is completely contained within the region
                    tempVariants.add(currentVariant)
                    currentVariant = reader.iterator().next()
                }
                else if(CreateMafVcf().variantPartiallyContainedStart(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)),currentVariant)) {
                    tempVariants.add(currentVariant)
                    break
                }
                else if(CreateMafVcf().variantPartiallyContainedEnd(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)), currentVariant)) {
                    tempVariants.add(currentVariant)
                    currentVariant = reader.iterator().next()
                }
                else if(CreateMafVcf().variantAfterRegion(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)), currentVariant)) {
                    //write out what is in tempVariants
                    if(tempVariants.isNotEmpty()) {
                        val fixedVariants = fixASMPositions(Pair(Position(regionChrom,regionStart),Position(regionChrom,regionEnd)), listOf(currentVariant))
                        refRangeToVariantContext.getOrPut(range, { mutableListOf() }).addAll(fixedVariants)
                        tempVariants.clear()
                    }
                    //move up Bed region
                    break
                }
                else { //this is the case if the Variant is behind the BED region
                    //move up Variant
                    currentVariant = reader.iterator().next()
                }
            }
        }

        return refRangeToVariantContext
    }

    // This fixes the asm position in the tempVariants list.  It should return an
    // ammended version of the list. COde is based on CreateMafVcf:convertGVCFRecordsToHVCFMetaData()
    // but only deals with the ASM positions portion of this code.
    fun fixASMPositions( region: Pair<Position,Position>, variants: List<VariantContext> ): List<VariantContext> {
        val fixedVariants = mutableListOf<VariantContext>()

        //Take the first and the last variantContext
        val firstVariant = variants.first()
        val lastVariant = variants.last()

        //val check strandedness of the variants
        val firstStrand = firstVariant.getAttributeAsString("ASM_Strand","+")

        val lastStrand = lastVariant.getAttributeAsString("ASM_Strand","+")
        //Resize the first and last variantContext ASM start and end based on the regions
        var newASMStart = CreateMafVcf().resizeVariantContext(firstVariant, region.first.position, firstStrand)
        if(newASMStart == -1) {
            newASMStart = if(firstStrand == "+") firstVariant.getAttributeAsInt("ASM_Start",region.first.position)
            else firstVariant.getAttributeAsInt("ASM_End",region.first.position)
        }

        var newASMEnd = CreateMafVcf().resizeVariantContext(lastVariant, region.second.position, lastStrand)
        if(newASMEnd == -1) {
            newASMEnd = if(lastStrand == "+") lastVariant.getAttributeAsInt("ASM_End",region.second.position)
            else lastVariant.getAttributeAsInt("ASM_Start",region.second.position)
        }

        val regions = CreateMafVcf().buildNewAssemblyRegions(newASMStart,newASMEnd,variants)

        //TODO - need to add the new regions to the variantContexts, replacing the old ASM_* values
        return fixedVariants

    }
    // this function takes a list of REferenceRange objects and a VarientContext Reader, and finds
    // positions in the gvcf file that overlap the ReferenceRange objects.  It returns a list of
    // VariantContext records that overlap the ReferenceRange objects.
    // Hmmm . this may not be correct.  I believe I need a map of ReferenceRange to VariantContext records

    fun findOverlappingRecords1(ranges:List<ReferenceRange>, reader:VCFFileReader):Map<ReferenceRange,List<VariantContext>> {
        // Assume the List<ReferenceRanges> is sorted by CHromosome and start position
        // Read the gvcf file, and for each VariantContext record, check if it overlaps.  Keep reading until we get to a variantContext
        // that overlaps with our current reference range.  When we find a context that overlaps, add it to the list of VariantContext records
        // for that range.  There may be more than 1 record that overlaps with a range.  Do not move to the next ReferenceRange until
        // the variant context record is beyond the end of the current ReferenceRange..  And do not move to the next variant Context record
        // if the current one extends beyond the current reference range.  It may get added to a list for both reference ranges.
        // This is a bit tricky, but I think it can be done.
        val refRangeToVariantContext = mutableMapOf<ReferenceRange, MutableList<VariantContext>>()
        var currentRange = ranges[0]
        var currentRangeIndex = 0
        var currentVariantContext = reader.iterator().next()
        //This is not totally correct but is a good start.
        // the chromosomes will change, and that is fine, we may need to move to the next refRange
        while (currentVariantContext != null) {
            if (currentVariantContext.contig != currentRange.contig) {
                //  At some point either the ReferenceRange
                // entry and/or the VariantContext entry is going to move up to the next
                // chromosome.  That means we need to move to the next ReferenceRange
                // and check if the VariantContext is on that chromosome.  If it is, we
                // need to check if it overlaps the ReferenceRange.  If it does, we add it
                // to the list of VariantContext records for that ReferenceRange.  If it doesn't
                // we move to the next VariantContext record.

                // how do I know if it is the VariantContext or the ReferenceRange that needs to move up?
                // I think it is the ReferenceRange that needs to move up.  If the VariantContext is on the
                // next chromosome, it will not overlap the current ReferenceRange.  If the ReferenceRange
                // is on the next chromosome, it will not overlap the current VariantContext.  So, I think
                // I need to move the ReferenceRange up to the next chromosome.
                currentRangeIndex++
                if (currentRangeIndex >= ranges.size) {
                    // We are done
                    break
                }
                currentRange = ranges[currentRangeIndex]
                // need to check again if the contigs match - not sure this is right, check Zack's code

            }
            if (currentVariantContext.start > currentRange.end) {
                // This variantContext is beyond the current ReferenceRange
                // Move to the next ReferenceRange
                currentRangeIndex++
                if (currentRangeIndex >= ranges.size) {
                    // We are done
                    break
                }
                currentRange = ranges[currentRangeIndex]
            } else if (currentVariantContext.end < currentRange.start) {
                // This variantContext is before the current ReferenceRange
                // Move to the next variantContext
                currentVariantContext = reader.iterator().next()
            } else {
                // This variantContext overlaps the current ReferenceRange
                // Add it to the list of VariantContext records for the current ReferenceRange
                val vcList = refRangeToVariantContext.getOrDefault(currentRange, mutableListOf())
                vcList.add(currentVariantContext)
                refRangeToVariantContext[currentRange] = vcList
                currentVariantContext = reader.iterator().next()
            }
        }

        return refRangeToVariantContext
    }
}