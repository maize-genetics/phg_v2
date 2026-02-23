package net.maizegenetics.phgv2.cli

import biokotlin.genome.AssemblyVariantInfo
import biokotlin.genome.MAFToGVCF
import biokotlin.genome.refDepth
import biokotlin.seq.NucSeq
import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextComparator
import htsjdk.variant.vcf.VCFFileReader
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.*
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
 * 3.  export the gvcf files from the tiledb database
 * 4.  If the reference sample is in the sample names list, and the reference gvcf file does not exist,
 *     create a gvcf of ref blocks for all reference ranges
 * 5.  Also from the hvcf files, create a list of ReferenceRange objects for each sample
 * 6.  Read the gvcf files, a sample at a time, and for each sample, pull the gvcf records
 *     that overlap the ReferenceRanges for that sample.  This will be a list of VariantContext records.
 * 7.  Merge the list of vcfRecords for each sample together, sorted by reference ranges, and write
 *     to a new gvcf file in the output directory.
 * 8.  The gvcf file will be written to the output directory.
 *
 */


data class HvcfVariant(val refRange: ReferenceRange,val sampleName: String, val hapId: String)


class Hvcf2Gvcf :
    CliktCommand(help = "Create g.vcf file for a PHG pathing h.vcf using data from existing PHG created g.vcf files") {
    private val myLogger = LogManager.getLogger(Hvcf2Gvcf::class.java)
    // These values come from BioKotlin:MAFToGVCF - they should remain consistent - what is a better
    // way to handle this?  Perhaps in Biokotlin they should be moved outside the class and thus
    // accessible as is done with MAFToGVCF.kt:refDepth, which we import to this class?

    val hvcfDir by option(
        "--hvcf-dir",
        help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB"
    )
        .required()

    val condaEnvPrefix by option(help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
        .default("")

    val outputDir by option(help = "Output directory for the gVCF files.  If not provided, the current working directory is used.")
        .default("")
    val batchSize by option(help = "Number of samples to export at a time.")
        .int()
        .default(1)

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        .default("")
    val referenceFile by option(help = "Path to local Reference FASTA file needed for sequence dictionary")
        .required()

    override fun run() {
        logCommand(this)

        val dbPath = dbPath.ifBlank {
            System.getProperty("user.dir")
        }

        // Verify the tiledbURI - verifyURI will throw an exception if the URI is not valid
        verifyURI(dbPath, "hvcf_dataset", condaEnvPrefix)

        // read and store ref file data for later use when creating the vcf sequence dictionary
        val time = System.nanoTime()
        val refSeq = CreateMafVcf().buildRefGenomeSeq(referenceFile)
        myLogger.info("Time to build refSeq: ${(System.nanoTime() - time) / 1e9} seconds")

        val timeHapIdMap = System.nanoTime()
        val hapIdAndRangeToSampleMap = createASMHapIdMap(dbPath)
        myLogger.info("Time to build hapId+Range -> Sample Map: ${(System.nanoTime() - timeHapIdMap)/1e9} seconds")


        buildGvcfFromHvcf(dbPath, refSeq, hapIdAndRangeToSampleMap,outputDir, hvcfDir, condaEnvPrefix, batchSize)
    }

    /**
     * Function to build gvcfs from an hvcf.
     *
     * If the hvcf is diploid it will output 2 separate gvcfs
     */
    private fun buildGvcfFromHvcf(
        dbPath: String,
        refSeq: Map<String, NucSeq>,
        hapIdAndRangeToSampleMap: Map<Pair<ReferenceRange,String>, List<HvcfVariant>>,
        outputDir: String,
        hvcfDir: String,
        condaEnvPrefix: String,
        batchSize: Int
    ) {
        // get list of hvcf files
        // walk the gvcf directory process files with g.vcf.gz extension

        File(hvcfDir).walk().filter { !it.isHidden && !it.isDirectory }
            .filter {
                it.name.endsWith(".h.vcf.gz") || it.name.endsWith(".h.vcf") ||
                        it.name.endsWith(".hvcf.gz") || it.name.endsWith(".hvcf")
            }
            .forEach { hvcfFile ->
                myLogger.info("buildGvcfFromHvcf: Processing hvcf file: ${hvcfFile.name}")
                var time = System.nanoTime()
                val sample = when {
                    (hvcfFile.name.endsWith(".h.vcf") || hvcfFile.name.endsWith(".h.vcf.gz")) ->
                        hvcfFile.name.substringBeforeLast(".h.vcf")

                    (hvcfFile.name.endsWith(".hvcf") || hvcfFile.name.endsWith(".hvcf.gz")) ->
                        hvcfFile.name.substringBeforeLast(".hvcf")

                    else -> error("Unexpected file extension")
                }

                val (gamete1Records, gamete2Records) = processHVCFtoVariantContext(sample, refSeq, hapIdAndRangeToSampleMap,outputDir, hvcfFile, dbPath, condaEnvPrefix, batchSize)

                myLogger.info("Time to processHVCFtoVariantContext: ${(System.nanoTime() - time) / 1e9} seconds")

                val (outputFileName1, outputFileName2) = buildOutputFileNames(sample, outputDir, gamete2Records)

                myLogger.info("buildGvcfFromHvcf: exporting VariantContexts to gvcf file: $outputFileName1")
                time = System.nanoTime()
                exportVariantContext(sample, gamete1Records, outputFileName1, refSeq, setOf())
                myLogger.info("Time to exportVariantContext for Gamete1: ${(System.nanoTime() - time) / 1e9} seconds")

                //Handle 2nd gamete if diploid
                if(gamete2Records.isNotEmpty()) {
                    myLogger.info("buildGvcfFromHvcf: exporting VariantContexts to gvcf file: $outputFileName2")
                    time = System.nanoTime()
                    exportVariantContext(sample, gamete2Records, outputFileName2, refSeq, setOf())
                    myLogger.info("Time to exportVariantContext: ${(System.nanoTime() - time) / 1e9} seconds")
                }
            }

    }

    /**
     * Function to build the output file names.  This works with both haploid and diploid input files.
     */
    private fun buildOutputFileNames(
        sampleName: String,
        outputDir: String,
        gamete2Records: List<VariantContext>
    ): Pair<String, String> {
        val outputFileName1 = if (gamete2Records.isEmpty()) {
            "$outputDir/${sampleName}.g.vcf"
        } else {
            "$outputDir/${sampleName}_1.g.vcf"
        }

        val outputFileName2 = if (gamete2Records.isEmpty()) {
            ""
        } else {
            "$outputDir/${sampleName}_2.g.vcf"
        }
        return Pair(outputFileName1, outputFileName2)
    }

    /**
     * Function to build the GVCF output sample name.
     *
     * This works with both haploid and diploid hvcfs
     */
    private fun buildOutputSampleNames(
        isHaploid: Boolean,
        outputSampleName: String
    ): Pair<String, String> {
        val outputSampleName1 = if (isHaploid) {
            outputSampleName
        } else {
            "${outputSampleName}_1"
        }
        val outputSampleName2 = if (isHaploid) {
            ""
        } else {
            "${outputSampleName}_2"
        }
        return Pair(outputSampleName1, outputSampleName2)
    }

    //The Map<Pair<RefRange,HapId>,List<HvcfVariant>>
    // the list does not matter as the sequence is the same but we should collect it in case
    //We retain all of them as we will need those boundaries later and it is easier to just keep them together
    fun createASMHapIdMap(dbPath: String): Map<Pair<ReferenceRange,String>, List<HvcfVariant>> {
        //check to see if there are hvcfFiles in the dbPath

        //If so we need to loop through them and extract out the HapIds and pair them with the refRanges.
        //TODO Eventually we can replace this with direct tileDB calls
        return File("${dbPath}/hvcf_files").walk().filter { !it.isHidden && !it.isDirectory }
            .filter {
                it.name.endsWith(".h.vcf.gz") || it.name.endsWith(".h.vcf") ||
                        it.name.endsWith(".hvcf.gz") || it.name.endsWith(".hvcf")
            }
            .flatMap { hvcfFile ->
                extractASMHapIds(hvcfFile)
            }
            .groupBy { Pair(it.refRange, it.hapId) } // We then can collect across all samples
    }


    /**
     * Function to pull out the assembly hvcf information so we can look it up later.
     * Here we can use the first alt allele as it is guaranteed to be haploid because it is coming from tiledb.
     */
    fun extractASMHapIds(file: File): List<HvcfVariant> {
        val hvcfReader = VCFFileReader(file, false)
        val sampleName = hvcfReader.header.sampleNamesInOrder.first()
        val hvcfVariants = hvcfReader.mapNotNull { context ->
            if (context.alternateAlleles.isEmpty() || context.alternateAlleles.any { it.isNoCall }) {
                null
            } else {
                //Here we can get the first alt allele as the ASM hvcfs are guaranteed to be single sample haploid
                val hapId = context.getAlternateAllele(0).displayString.removeSurrounding("<", ">")
                val refRange = ReferenceRange(context.contig, context.start, context.end)
                HvcfVariant(refRange, sampleName,hapId)
            }
        }

        return hvcfVariants
    }


    /**
     * Function to process the HVCF file into a set of Variant contexts
     */
    private fun processHVCFtoVariantContext(
        outputSampleName: String,
        refSeq: Map<String, NucSeq>,
        hapIdAndRangeToSampleMap: Map<Pair<ReferenceRange,String>, List<HvcfVariant>> ,outputDir: String, hvcfFile: File,
        dbPath: String,
        condaEnvPrefix: String,
        batchSize: Int
    ): Pair<List<VariantContext>,List<VariantContext>> {

        val reader = VCFFileReader(hvcfFile, false)

        val header = reader.fileHeader
        val altHeaders = parseALTHeader(header = header)
        val sampleNames = altHeaders.values.map { it.sampleName() }.toSet()
        val contigNames = header.contigLines.map { it.id } // contig names needed for sorting
        //This will also export all the gvcfs that are not already exported
        val sampleNameToVCFMap = buildSampleNameToGVCFMap(sampleNames, outputDir, dbPath, condaEnvPrefix, refSeq, batchSize)


        // Using the hvcfFileReader, walk the hvcf file and for each entry create
        // a ReferenceRange object and add that to a list of Reference Range objects
        // for the sample to which it applies.  You will get the sample by indexing the
        // altHeaders map with the key of the current record's ID.
        val (isHaploid, hvcfVariantsList) = extractHVCFVariantsFromSample(reader, hapIdAndRangeToSampleMap)


        //Instead of making a big map, make 2 maps if diploid
        val gvcfSampleNameToRefRangesGamete1 = hvcfVariantsList.map { it.first }.filterNotNull()
            .groupBy({ it.sampleName },{ it.refRange })

        val gvcfSampleNameToRefRangesGamete2 = if(isHaploid) {
            emptyMap()
        }
        else {
            hvcfVariantsList.map { it.second }.filterNotNull()
                .groupBy ({ it.sampleName },{ it.refRange })
        }

        val gvcfSampleNames = gvcfSampleNameToRefRangesGamete1.keys + gvcfSampleNameToRefRangesGamete2.keys

        return extractGVCFVariantsForRefRanges(
            gvcfSampleNames,
            sampleNameToVCFMap,
            isHaploid,
            outputSampleName,
            refSeq,
            contigNames,
            gvcfSampleNameToRefRangesGamete1,
            gvcfSampleNameToRefRangesGamete2
        )
    }

    /**
     * Function to extract the gvcf variants for the reference ranges that overlap the hvcf.
     */
    private fun extractGVCFVariantsForRefRanges(
        gvcfSampleNames: Set<String>,
        sampleNameToVCFMap: MutableMap<String, String>,
        isHaploid: Boolean,
        outputSampleName: String,
        refSeq: Map<String, NucSeq>,
        contigNames: List<String>,
        gvcfSampleNameToRefRangesGamete1: Map<String, List<ReferenceRange>>,
        gvcfSampleNameToRefRangesGamete2: Map<String, List<ReferenceRange>>
    ): Pair<List<VariantContext>, List<VariantContext>> {
        val gamete1VariantContexts = mutableListOf<VariantContext>()
        val gamete2VariantContexts = mutableListOf<VariantContext>()

        for (gvcfSampleName in gvcfSampleNames) {
            myLogger.info("processHVCFtoVariantContext: processing sample: $gvcfSampleName")
            val time = System.nanoTime()
            check(sampleNameToVCFMap.containsKey(gvcfSampleName)) { "Sample $gvcfSampleName is not in the sampleNameToVCFMap.  This should not happen as we should have exported all gvcfs at the beginning of this function." }
            val gvcfFile = sampleNameToVCFMap[gvcfSampleName]!! // tiledb wrote with extension .vcf
            val gvcfReader = VCFFileReader(File(gvcfFile), false)

            val gvcfVariants = mutableListOf<VariantContext>()
            gvcfReader.use { reader ->
                for (vc in reader) {
                    gvcfVariants.add(vc)
                }
            }
            gvcfReader.close()

            val (outputSampleName1, outputSampleName2) = buildOutputSampleNames(isHaploid, outputSampleName)


            //Now find overlapping snps for gamete 1
            if (gvcfSampleNameToRefRangesGamete1.containsKey(gvcfSampleName)) {
                val ranges = gvcfSampleNameToRefRangesGamete1[gvcfSampleName]!!

                val variants =
                    findOverlappingRecordsForSample(outputSampleName1, ranges, gvcfVariants, refSeq).values.flatten()
                gamete1VariantContexts.addAll(variants)
            }

            if (gvcfSampleNameToRefRangesGamete2.isNotEmpty() && gvcfSampleNameToRefRangesGamete2.containsKey(
                    gvcfSampleName
                )
            ) {
                val ranges = gvcfSampleNameToRefRangesGamete2[gvcfSampleName]!!

                val variants =
                    findOverlappingRecordsForSample(outputSampleName2, ranges, gvcfVariants, refSeq).values.flatten()
                gamete2VariantContexts.addAll(variants)
            }

            myLogger.info("Time to process sampleToRefRanges for sample $gvcfSampleName : ${(System.nanoTime() - time) / 1e9} seconds")
        }

        //Make sure the outputs are sorted correctly
        val outputGamete1VariantContexts = gamete1VariantContexts.sortedWith(VariantContextComparator(contigNames))
        val outputGamete2VariantContexts = gamete2VariantContexts.sortedWith(VariantContextComparator(contigNames))

        return Pair(outputGamete1VariantContexts, outputGamete2VariantContexts)
    }


    /**
     * Function to extract the HVCF variants for each sample.  This will give us a pair of variants if diploid.
     */
    fun extractHVCFVariantsFromSample(
        reader: VCFFileReader,
        hapIdAndRangeToSampleMap: Map<Pair<ReferenceRange, String>, List<HvcfVariant>>
    ): Pair<Boolean,List<Pair<HvcfVariant?, HvcfVariant?>>> {
        var rangesSkipped = 0
        var totalHvcfVariants = 0
        var isHaploid = true
        val hvcfRecordsList = reader.mapNotNull { context ->
            totalHvcfVariants++
            //This might need to be relaxed as it might not handle if one of the haps is missing
            if (context.alternateAlleles.isEmpty() || context.alternateAlleles.any { it.isNoCall }) {
                rangesSkipped++
                null
            } else {
                val (isVariantHaploid, variants) = processCurrentHVCFVariant(context, hapIdAndRangeToSampleMap)
                if(!isVariantHaploid) {
                    isHaploid = false
                }
                variants
            }
        }

        myLogger.info("Hvcf2Gvcf:processHVCFtoVariantContext: rangesSkipped = $rangesSkipped, totalHvcfVariants = $totalHvcfVariants")

        return Pair(isHaploid, hvcfRecordsList)
    }

    /**
     * Function to process the HVCF record and extract out the variants
     */
    fun processCurrentHVCFVariant(
        context: VariantContext,
        hapIdAndRangeToSampleMap: Map<Pair<ReferenceRange, String>, List<HvcfVariant>>
    ): Pair<Boolean,Pair<HvcfVariant?, HvcfVariant?>> {
        var isHaploid = true
        //Need to extract out both the first genotype and the second's hapIds
        //We assume that we only have one sample in this hvcf
        val refRange = ReferenceRange(context.contig, context.start, context.end)

        val genotypes = context.genotypes.first()
        val hapId1 = genotypes.getAllele(0).displayString.removeSurrounding("<", ">")
        val hapId2 = if ((genotypes?.alleles?.size ?: 0) > 1) {
            genotypes.getAllele(1).displayString.removeSurrounding("<", ">")
            isHaploid = false
        } else {
            hapId1
        }

        //Get the first hvcfVariant for each gamete.
        // If There are more than one its fine as they have the same sequence at this reference range so they
        // should have the same variants
        //If they are null it means that it is a missing call...
        //By going to the hapIdAndRangeMap from the assembly hvcfs we make sure we are not pulling out incorrect records.
        val hvcfVariants1 = hapIdAndRangeToSampleMap[Pair(refRange, hapId1)]?.first()
        val hvcfVariants2 = hapIdAndRangeToSampleMap[Pair(refRange, hapId2)]?.first()

        return Pair(isHaploid,Pair(hvcfVariants1, hvcfVariants2))
    }

    /**
     * This function will build the sampleNameToGVCF map and
     * will export any missing gvcf files and make a reference gvcf if missing
     */
    private fun buildSampleNameToGVCFMap(
        sampleNames: Set<String>,
        outputDir: String,
        dbPath: String,
        condaEnvPrefix: String,
        refSeq: Map<String, NucSeq>,
        batchSize: Int
    ): MutableMap<String, String> {
        // This checks for existing gvcf files, and if they don't exist, exports them
        // Get list of sample names and the gvcf file names for each sample
        val sampleNameToVCFMap = getSampleNameToVCFMap(sampleNames.toList(), outputDir).toMutableMap()

        // There is no gvcf for the reference, so we need to create one if
        // the reference is in the sampleNames list.  The query below gets the
        // ref sample name from the agc file (which is <dbPath>/assemblies.agc)
        val refSampleName = retrieveRefSampleName(dbPath, condaEnvPrefix)

        val refFileExists = listOf(
            File(outputDir, "$refSampleName.vcf"),
            File(outputDir, "$refSampleName.g.vcf"),
            File(outputDir, "$refSampleName.g.vcf.gz")
        ).any { it.exists() }

        // If the refSampleName is in the sampleNames list, and the refGvcfFile does not exist,
        // create a gvcf of ref blocks for all reference ranges. Because there may be multiple
        // hvcf files, and each one may have different reference ranges, we need to create a
        // gvcf file for the reference sample that contains all the reference ranges from the bed file
        if (sampleNames.contains(refSampleName) && !refFileExists) {
            // find the bed file.  When the reference was loaded, the bed file was copied
            // to the dbPath/reference folder.  Look for any file with extension .bed in
            // folder.
            val bedFile = File("$dbPath/reference").walk().filter { it.name.endsWith(".bed") }.toList()
                .firstOrNull()?.absolutePath
            check(bedFile != null) { "No bed file found in $dbPath/reference, cannot pull ref haplotypes" }

            // create a gvcf of ref blocks for all reference ranges
            createRefGvcf(refSampleName, bedFile, refSeq, outputDir)
            sampleNameToVCFMap[refSampleName] = "$outputDir/$refSampleName.vcf"
        }
        val missingSamples = sampleNames - sampleNameToVCFMap.keys

        if (missingSamples.isNotEmpty()) {
            myLogger.info("calling exportGvcfFiles with ${missingSamples.size} missing samples")

            runBlocking {
                val exportSuccess = exportSamplesGvcfFiles(missingSamples, outputDir, dbPath, condaEnvPrefix, batchSize)
                if (!exportSuccess) {
                    myLogger.error("Failed to export gvcf files for missing samples")
                    throw IllegalStateException("Failed to export gvcf files for missing samples")
                }
            }

            // Add missing samples to the map - should be MissingSample -> <outputDir>/<MissingSample>.vcf
            missingSamples.forEach { sampleName ->
                sampleNameToVCFMap[sampleName] = "$outputDir/$sampleName.vcf"
            }
        }
        return sampleNameToVCFMap
    }

    // get the vcf file name for all samples in the list.  Should be a map<SampleName, fileName>
    private fun getSampleNameToVCFMap(sampleNames: List<String>, outputDir: String): Map<String, String> {
        val directory = File(outputDir)

        // Return a map of sampleName to the actual VCF file found
        return sampleNames.mapNotNull { sampleName ->
            // Check if any of the expected files exist
            val vcfFile = listOf(
                File(directory, "$sampleName.vcf"),
                File(directory, "$sampleName.g.vcf"),
                File(directory, "$sampleName.g.vcf.gz")
            ).firstOrNull { it.exists() } // Get the first file that exists

            // If a file exists, return the pair (sampleName, vcfFile name), otherwise return null
            if (vcfFile != null) {
                sampleName to vcfFile.absolutePath
            } else {
                null
            }
        }.toMap() // Convert the list of pairs to a map
    }

    // This function creates a gvcf of refBlocks based on the ReferenceRanges for the ref sample
    fun createRefGvcf(refSampleName: String, bedFile: String, refSeq: Map<String, NucSeq>, outputDir: String) {

        val totalTime = System.nanoTime()
        val mafToGvcf = MAFToGVCF()
        val variantList = mutableListOf<AssemblyVariantInfo>()
        // Read the bedZFile, use it to create a list of VariantInfo records
        bufferedReader(bedFile).use { reader ->
            // buildRefBLockVariantInfo for each range in the bed file
            reader.forEachLine { line ->
                val fields = line.split("\t")
                val contig = fields[0]
                val start = fields[1].toInt()
                val end = fields[2].toInt()
                variantList += buildRefBlockVariantInfo(
                    refSeq,
                    contig,
                    Pair<Int, Int>(start + 1, end),
                    contig, // asm data is same as ref data
                    Pair<Int, Int>(start + 1, end),
                    "+" // strand is always positive for ref
                )
            }
        }

        val vcs = mafToGvcf.createVariantContextsFromInfo(refSampleName, variantList, false, false, 0)

        val refGvcfFile = "$outputDir/${refSampleName}.vcf"
        // Once created, we export the ref gvcf file.  The extension is ".vcf" as that
        // is consistent with the extension used by tiledbvcf when exporting the gvcf files
        exportVariantContext(refSampleName, vcs, refGvcfFile, refSeq, setOf())
        myLogger.info("Time to createRefGvcf: ${(System.nanoTime() - totalTime) / 1e9} seconds")
    }

    // This function exists in BioKotlin, but is private.  If it is moved to public,
    // we can remove this and access the BioKotlin version.
    fun buildRefBlockVariantInfo(
        refSequence: Map<String, NucSeq>,
        chrom: String,
        currentRefBlockBoundaries: Pair<Int, Int>,
        assemblyChrom: String,
        currentAssemblyBoundaries: Pair<Int, Int>,
        assemblyStrand: String
    ): AssemblyVariantInfo {
        // -1, NucSeq is 0-based
        return AssemblyVariantInfo(
            chrom, currentRefBlockBoundaries.first, currentRefBlockBoundaries.second, "REF",
            refSequence[chrom]!!.get(currentRefBlockBoundaries.first - 1).toString(), ".", false,
            refDepth, assemblyChrom, currentAssemblyBoundaries.first, currentAssemblyBoundaries.second, assemblyStrand
        )
    }

    suspend fun exportSamplesGvcfFiles(
        sampleNames: Set<String>,
        outputDir: String,
        dbPath: String,
        condaEnvPrefix: String,
        batchSize: Int
    ): Boolean {

        val processingChannel = Channel<Deferred<Boolean>>(25)

        CoroutineScope(Dispatchers.IO).launch {

            sampleNames.chunked(batchSize).forEach { sampleList ->
                processingChannel.send(async {
                    exportGvcfFiles(
                        sampleList,
                        outputDir,
                        dbPath,
                        condaEnvPrefix
                    )
                })
            }

            processingChannel.close()

        }

        for (deferred in processingChannel) {
            try {
                deferred.await()
            } catch (e: Exception) {
                myLogger.error("Error exporting gvcf files: ${e.message}")
                return false
            }
        }

        return true

    }

    // function to export from tiledb the gvcf files if they don't exist
    // If we get the tiledb-java API working for MAC, we can hold these in memory
    // while processing and skip the export.
    private fun exportGvcfFiles(
        sampleNames: List<String>,
        outputDir: String,
        dbPath: String,
        condaEnvPrefix: String
    ): Boolean {

        // Set up the conda environment portion of the command
        val condaCommandPrefix = if (condaEnvPrefix.isNotBlank()) {
            mutableListOf("conda", "run", "-p", condaEnvPrefix)
        } else {
            mutableListOf("conda", "run", "-n", "phgv2-tiledb")
        }

        val dataCommand = mutableListOf(
            "tiledbvcf",
            "export",
            "--uri",
            "$dbPath/gvcf_dataset",
            "--output-format",
            "v",
            "--sample-names",
            sampleNames.joinToString(","),
            "--output-dir",
            outputDir
        )

        // Join the conda and data portions of the command
        val command = condaCommandPrefix + dataCommand
        val builder = ProcessBuilder(command)

        val redirectError = "$outputDir/export_gvcf_error.log"
        val redirectOutput = "$outputDir/export_gvcf_output.log"
        builder.redirectOutput(File(redirectOutput))
        builder.redirectError(File(redirectError))

        // Log the command
        myLogger.info("Command: " + builder.command().joinToString(" "))

        // Start the process
        val process = builder.start()
        val error = process.waitFor()

        // Handle error if the process fails
        if (error != 0) {
            myLogger.error("tiledbvcf export returned error code $error")
            throw IllegalStateException(
                "Error running tiledbvcf export of dataset $dbPath/gvcf_dataset: $error"
            )
        }

        return true

    }


    // Process the gvcf records for a sample, and return a map of ReferenceRange to VariantContext
    // These are processed per-sample, then per-chromosome, to reduce memory usage
    private fun findOverlappingRecordsForSample(
        outputSampleName: String,
        ranges: List<ReferenceRange>,
        variantContexts: List<VariantContext>,
        refSeq: Map<String, NucSeq>
    ): MutableMap<ReferenceRange, MutableList<VariantContext>> {
        // Query AGC to get sequence for the sample
        println("findOverlappingRecordsForSample: outputSampleName = $outputSampleName")
        // split ranges and variantContexts by chromosome
        val time = System.nanoTime()
        val overlappingRecords = mutableMapOf<ReferenceRange, MutableList<VariantContext>>()
        val rangesByChrom = ranges.groupBy { it.contig }
        val variantContextsByChrom = variantContexts.groupBy { it.contig }

        // call findOverlappingRecords for each chromosome, add the results to overlappingRecords
        rangesByChrom.keys.forEach { chrom ->
            val chromRanges = rangesByChrom[chrom] ?: emptyList()
            val chromVariants = variantContextsByChrom[chrom] ?: emptyList()
            myLogger.info("findOverlappingRecordsForSample: Processing ${chromRanges.size} ranges and ${chromVariants.size} variants for chromosome $chrom")
            val chromOverlappingRecords = findOverlappingRecords(outputSampleName, chromRanges, chromVariants, refSeq)
            overlappingRecords.putAll(chromOverlappingRecords)
        }
        myLogger.info("Time to run findOverlappingRecordsForSample: ${(System.nanoTime() - time) / 1e9} seconds")
        return overlappingRecords
    }


    // This is based on CreateMafVCF:convertGVCFToHVCFForChrom() which determines if a variant or part
    // of a variant is in a reference range.  It differs in that we are not creating hvcf meta data,
    // but rather, at this stage, are merely finding the overlapping variants
    private fun findOverlappingRecords(
        sample: String,
        ranges: List<ReferenceRange>,
        variantContexts: List<VariantContext>,
        refSeq: Map<String, NucSeq>
    ): MutableMap<ReferenceRange, MutableList<VariantContext>> {

        val refRangeToVariantContext =
            mutableMapOf<ReferenceRange, MutableList<VariantContext>>() // this will be returned
        var currentVariantIdx = 0

        for (range in ranges) {
            val regionStart = range.start
            val regionEnd = range.end
            val regionChrom = range.contig
            val tempVariants = mutableListOf<VariantContext>()

            while (currentVariantIdx < variantContexts.size) {
                val currentVariant = variantContexts[currentVariantIdx]

                val region = Pair(
                    Position(regionChrom, regionStart),
                    Position(regionChrom, regionEnd)
                )

                when {
                    VariantContextUtils.bedRegionContainedInVariant(region, currentVariant) -> {
                        // This is the case where the region is completely contained within the variant,
                        // meaning the variant may overlap the region.  We need to adjust the asm positions
                        val fixedVariants = fixPositions(
                            sample,
                            Pair(Position(regionChrom, regionStart), Position(regionChrom, regionEnd)),
                            listOf(currentVariant),
                            refSeq
                        )
                        val currentVariants = refRangeToVariantContext.getOrPut(range) { mutableListOf() }
                        currentVariants.addAll(fixedVariants)
                        tempVariants.clear()
                        break
                    }

                    VariantContextUtils.variantFullyContained(region, currentVariant) -> {
                        // This is the case where the variant is completely contained within the region
                        tempVariants.add(currentVariant)
                        currentVariantIdx++
                    }

                    VariantContextUtils.variantPartiallyContainedStart(region, currentVariant) -> {
                        tempVariants.add(currentVariant)
                        break
                    }

                    VariantContextUtils.variantPartiallyContainedEnd(region, currentVariant) -> {
                        tempVariants.add(currentVariant)
                        currentVariantIdx++
                    }

                    VariantContextUtils.variantAfterRegion(region, currentVariant) -> {
                        // write data from tempVariants
                        if (tempVariants.isNotEmpty()) {
                            val fixedVariants = fixPositions(
                                sample,
                                Pair(Position(regionChrom, regionStart), Position(regionChrom, regionEnd)),
                                tempVariants,
                                refSeq
                            )
                            val currentVariants = refRangeToVariantContext.getOrPut(range) { mutableListOf() }
                            currentVariants.addAll(fixedVariants)
                            tempVariants.clear()
                        }
                        // move up bed region
                        break
                    }

                    else -> {
                        // move up variant
                        currentVariantIdx++
                    }
                }
            }

            // Process the last variants in the list
            if (tempVariants.isNotEmpty()) {
                val fixedVariants = fixPositions(
                    sample,
                    Pair(Position(regionChrom, regionStart), Position(regionChrom, regionEnd)),
                    tempVariants,
                    refSeq
                )
                val currentVariants = refRangeToVariantContext.getOrPut(range) { mutableListOf() }
                currentVariants.addAll(fixedVariants)
                tempVariants.clear()
            }
        }

        return refRangeToVariantContext
    }

    // This fixes both the ref start/end and the asm positions in the tempVariants list.  It should return an
    // ammended version of the list. Code is based on CreateMafVcf:convertGVCFRecordsToHVCFMetaData()
    // but only deals with the ASM positions portion of this code.  In addition, it adds code
    // to adjust the variant's start/end positions. Indels are not currently handled.
    // Is the variant is a SNP, we would not get into this function, as the SNP would be fully contained
    // in the region.  This function is only called when the variant is partially contained in the region.

    private fun fixPositions(
        sampleName: String,
        region: Pair<Position, Position>,
        variants: List<VariantContext>,
        refSeq: Map<String, NucSeq>
    ): List<VariantContext> {
        val fixedVariants = mutableListOf<VariantContext>()

        // TODO ADD INDEL SUPPORT ??
        // Take the first and the last variantContext
        val firstVariant = variants.first()
        val lastVariant = variants.last()

        // val check strandedness of the variants
        val firstStrand = firstVariant.getAttributeAsString("ASM_Strand", "+")

        val lastStrand = lastVariant.getAttributeAsString("ASM_Strand", "+")

        // Resize the first and last variantContext ref/ASM start and end based on the regions
        // This returns a list of Pairs of Ints.  The first Int is the reference position
        // and the second Int is the ASM position.  The first pair is for the start of the first variant,
        // the second pair is for the end of the last variant..  If there is only 1 variant, the 2 should be the same.

        val newGvcfPositions = resizeVCandASMpositions(
            Pair(firstVariant, lastVariant),
            Pair(region.first.position, region.second.position),
            Pair(firstStrand, lastStrand)
        )

        // At this point, we have changes for the first and last regions.  If the list size
        // is only 1, we change it based on newStartPositions and newEndPositions
        // If there are multiple, we change the first and last entries.  The first entry gets its start changed,
        // The last entry gets end values changed.
        if (variants.size == 1) {

            try {
                // createRefRangeVC will update the positions of the variant as well as the ref allele
                // value based on the adjusted ref position.
                val vcb = createRefRangeVC(
                    refSeq,
                    sampleName,
                    Position(firstVariant.contig, newGvcfPositions[0].first),
                    Position(firstVariant.contig, newGvcfPositions[1].first),
                    Position(firstVariant.contig, newGvcfPositions[0].first),
                    Position(firstVariant.contig, firstVariant.end),
                    firstVariant.getAttributeAsString("ASM_Strand", "+")
                )

                fixedVariants.add(vcb)
            } catch (exc: Exception) {
                myLogger.error("Error: Size=1 updating first variant with newGvcfPositions: ${newGvcfPositions}, created from region: ${region}")
                myLogger.error("  the firstVariant that FAILS: ${firstVariant}")
                throw exc
            }

        } else {
            // update the first and last variants, leaving those
            // in the middle unchanged.  CreateRefRangeVC will update the positions of the variant
            // as well as the ref allele
            try {
                // update the first variant in the list
                val vcb = createRefRangeVC(
                    refSeq,
                    sampleName,
                    Position(firstVariant.contig, newGvcfPositions[0].first),
                    Position(firstVariant.contig, firstVariant.end),
                    Position(firstVariant.contig, newGvcfPositions[0].second),
                    Position(firstVariant.contig, firstVariant.end),
                    firstVariant.getAttributeAsString("ASM_Strand", "+")
                )
                fixedVariants.add(vcb)
            } catch (exc: Exception) {
                myLogger.error("Error: Size > 1, updating first variant with newGvcfPositions: ${newGvcfPositions}, created from region: ${region}")
                myLogger.error("  the firstVariant that FAILS: ${firstVariant}")
                throw exc
            }

            try {
                // Update the last variant in the list
                if (variants.size > 2) {
                    fixedVariants.addAll(variants.subList(1, variants.size - 1))
                }
                val vcb = createRefRangeVC(
                    refSeq,
                    sampleName,
                    Position(lastVariant.contig, lastVariant.start),
                    Position(lastVariant.contig, newGvcfPositions[1].first),
                    Position(lastVariant.contig, lastVariant.getAttributeAsInt("ASM_Start", lastVariant.start)),
                    Position(lastVariant.contig, newGvcfPositions[1].second),
                    lastVariant.getAttributeAsString("ASM_Strand", "+")
                )
                fixedVariants.add(vcb)
            } catch (exc: Exception) {
                myLogger.error("Error: Size > 1, updating last variant with newGvcfPositions: ${newGvcfPositions}, created from region: ${region}")
                myLogger.error("  the lastVariant that FAILS: ${lastVariant}")
                throw exc
            }

        }

        return fixedVariants
    }


    //  Based on CreateMafVcf:resizeVariantContext() - but this version deals with both
    // the reference start/end as well as the ASM_* positions
    // This new version: returns a List of Pairs of Ints.  The first Int is the reference position
    // and the second Int is the ASM position.  This first pair is for the start of the first variant,
    // the second pair is for the end of the last variant.  If there is only 1 variant, the 2 should be
    // the same.
    // NOTE: reverse-strand processing is only relevant to the asm coordinates.
    // The reference coordinates will always be based on the forward strand.
    fun resizeVCandASMpositions(
        variants: Pair<VariantContext, VariantContext>,
        positions: Pair<Int, Int>,
        strands: Pair<String, String>
    ): List<Pair<Int, Int>> {

        val updatedPositions = mutableListOf<Pair<Int, Int>>()
        var refAsmPos_first = Pair<Int, Int>(-1, -1)
        var refAsmPos_last = Pair<Int, Int>(-1, -1)
        val firstVariant = variants.first
        val lastVariant = variants.second

        //check to see if the variant is either a RefBlock or is a SNP with equal lengths
        if (VariantContextUtils.isVariantResizable(firstVariant)) {
            // if the position is < the start of the variant, then we return <variant.start,asm_start>

            //  The "position" must be modified in the above with an offset.
            // these 2 checks verify the position is within the variant range
            // if the strand is +, then we return position + offset
            // if the strand is -, then we return position - offset
            refAsmPos_first = when {
                // Both postions.first and firstVariant.start are reference positions
                // We are adjusting the reference position and the asm position based on how far
                // away the reference position is from the variant start.
                positions.first < firstVariant.start -> {
                    // The reference position starts before the variant, so we keep the new gvcf
                    // entry start equal to the current variant start
                    Pair(firstVariant.start, firstVariant.getAttributeAsInt("ASM_Start", firstVariant.start))
                }

                strands.first == "+" -> {
                    // This and the case below are hit when the ref start position is within the variant
                    val offset = positions.first - firstVariant.start

                    // We need to offset the ASM_Start, by the difference between the position and the variant start
                    // However, the ref position should be the same as the ref range start as it is equal to or
                    // greater than the variant start
                    Pair(positions.first, firstVariant.getAttributeAsInt("ASM_Start", firstVariant.start) + offset)
                }

                strands.first == "-" -> {
                    val offset = positions.first - firstVariant.start

                    // offset for reverse strand
                    Pair(positions.first, firstVariant.getAttributeAsInt("ASM_Start", firstVariant.end) - offset)
                }

                else -> {
                    Pair(-1, -1) // should never get here
                }
            }
        } else {
            // not resizable, so only change the ASM_* values
            val newASMStart = if (strands.first == "+") firstVariant.getAttributeAsInt("ASM_Start", positions.first)
            else firstVariant.getAttributeAsInt("ASM_End", positions.first)
            refAsmPos_first = Pair(firstVariant.start, newASMStart)
        }
        updatedPositions.add(refAsmPos_first)

        // Processing the last variant that overlaps the reference range.
        // This could be the same as the first variant if the list only has 1 variant
        // The updated values depend on whether the refRange overlaps the beginning of the
        // variant, the end, or is completely contained within the variant.
        if (VariantContextUtils.isVariantResizable(lastVariant)) {

            refAsmPos_last = when {
                positions.second >= lastVariant.end -> {
                    // The end of the ref range is beyond the end of the variant, so
                    // keep the variant.end and ASM_End the same
                    Pair(lastVariant.end, lastVariant.getAttributeAsInt("ASM_End", lastVariant.end))
                }

                strands.first == "+" -> {
                    // ref range ends before the lastVariant.end and is forward strand
                    val offset = positions.second - lastVariant.start

                    // RefRanges ends is before the lastVariant.end
                    // We need to offset the ASM_End, by the difference between the position and the variant start
                    Pair(positions.second, lastVariant.getAttributeAsInt("ASM_Start", lastVariant.start) + offset)
                }

                strands.first == "-" -> {
                    val offset = positions.second - lastVariant.start

                    // variant end at or beyond the ref range end.  Move the ASM_Start by the offset
                    Pair(positions.second, lastVariant.getAttributeAsInt("ASM_Start", firstVariant.end) - offset)
                }

                else -> {
                    Pair(-1, -1) // should never get here
                }
            }

        } else {
            // not resizable so only change the ASM_* values
            val newASMEnd = if (strands.second == "+") lastVariant.getAttributeAsInt("ASM_End", positions.second)
            else lastVariant.getAttributeAsInt("ASM_Start", positions.second)
            refAsmPos_last = Pair(lastVariant.start, newASMEnd)
        }
        updatedPositions.add(refAsmPos_last)
        return updatedPositions
    }
}
