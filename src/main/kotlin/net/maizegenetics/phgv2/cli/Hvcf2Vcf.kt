package net.maizegenetics.phgv2.cli

import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriter
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.utils.HvcfVariant
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.VCFConversionUtils
import net.maizegenetics.phgv2.utils.addSequenceDictionary
import net.maizegenetics.phgv2.utils.createGenericHeader
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.TreeMap

data class HvcfRangeHapIdSampleGamete(val refRange: ReferenceRange, val hapId: String, val sampleGametes: List<SampleGamete>)


class Hvcf2Vcf:
    CliktCommand(help = "Create vcf file for a PHG pathing h.vcf using data from existing PHG created vcf files") {

    private val myLogger = LogManager.getLogger(Hvcf2Vcf::class.java)

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        .default("")


    val hvcfDir by option(help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        .required()

    val pangenomeVcfFile by option(help = "Path to the VCF file containing all the PHG SNPs.  This is typically created by running merge-gvcf on the Assembly Gvcf files.")
        .required()

    val outputFile by option(help = "Output file.")
        .required()

    //This is needed to create the refSeqDictionary
    val referenceFile by option(help = "Path to local Reference FASTA file needed for sequence dictionary")
        .required()


    override fun run() {
        processHVCFAndBuildVCF(dbPath, hvcfDir, pangenomeVcfFile, outputFile, referenceFile)
    }

    /**
     * Function to process all the hvcfs in hvcfDir and extract out the donorVCF variants and write to the outputfile
     */
    fun processHVCFAndBuildVCF(dbPath: String, hvcfDir: String, donorVcfFile: String, outputFile: String, referenceFile: String) {
        //load in the refSeq
        myLogger.info("Building RefGenome Sequence")
        val refSeq = CreateMafVcf().buildRefGenomeSeq(referenceFile)

        myLogger.info("Creating the ASM HapId map from the Assembly hVCF files.")
        //Load in the ASM hapId map so we can make sure we pull out the right SNPs.
        //This is in form (RefRange,ASMName) -> hapID to allow for easy lookup
        val asmHapIdMap = VCFConversionUtils.createASMNameAndRefRangeMap(dbPath)

        myLogger.info("Creating map keeping track of the refRanges and HapIds.")
        //We then need to load in the hvcfFiles and create a map of refRange+SampleName -> HapId1 + HapId2
        //Or maybe it needs to be refRange+HapId -> List<Pair<SampleName,Gamete>>
        val refRangeAndHapIdMap = createRangeHapMapToSampleGamete(hvcfDir)

        val ranges = refRangeAndHapIdMap.map {it.key.first}.toSet()

        myLogger.info("Building Position to refRange Lookup")
        //We also need to build a fast Position -> refRange lookup
        val positionToRangeMap = buildPositionToRefRangeMap(ranges)

        myLogger.info("Walking through the Merged VCF file and exporting out the new VCF file.")
        extractVcfAndExport(asmHapIdMap, refRangeAndHapIdMap, positionToRangeMap, donorVcfFile, outputFile, refSeq)
    }


    /**
     * Function to create a RangeMap+HapId -> List<SampleGamete> for the hvcfs to extract out the SNPs from
     */
    fun createRangeHapMapToSampleGamete(hvcfDir: String): Map<Pair<ReferenceRange, String>, List<SampleGamete>> {
        //walk through the hvcfDir and process each hvcf file
        return File(hvcfDir).walk().filter { !it.isHidden && !it.isDirectory }
            .filter {
                it.name.endsWith(".h.vcf.gz") || it.name.endsWith(".h.vcf") ||
                        it.name.endsWith(".hvcf.gz") || it.name.endsWith(".hvcf")
            }
            .flatMap { hvcfFile ->
                processSingleHvcfFile(hvcfFile)
            }.groupBy({Pair(it.refRange, it.hapId)}, {
                it.sampleGametes
            }).mapValues { entry ->
                entry.value.flatten()
            }
    }

    /**
     * Function to process a single HVCF File to determine which sample gametes are in each haplotype for each reference range.
     * This is needed to pull out the correct SNPs from the donor VCF file.
     */
    fun processSingleHvcfFile(hvcfFile: File): List<HvcfRangeHapIdSampleGamete> {
        val hvcfReader = VCFFileReader(hvcfFile, false)

        //This could be a multisample hvcf file so we need to collect all of the gametes
        return hvcfReader.mapNotNull { context ->
            processSingleHvcfVariant(context)
        }.flatten()
    }

    /**
     * Function to process a single Hvcf Variant extracting out which hapIds match each sampleGamete
     */
    fun processSingleHvcfVariant(context: VariantContext): List<HvcfRangeHapIdSampleGamete> {
        //Pull out the reference range info
        val refRange = ReferenceRange(context.contig, context.start, context.end)

        //Get the genotypes for each sample
        val genotypes = context.genotypes.flatMap { genotype ->
            val baseSampleName = genotype.sampleName

            //walk through each sample and convert to a pair Allele, sampleGamete
            genotype.alleles.mapIndexed { index, allele ->
                val cleanedAllele = allele.displayString.removeSurrounding("<", ">")
                Pair(cleanedAllele, SampleGamete(baseSampleName, index))
            }
        }

        //Now we need to group the cleaned alleles and put the sample gametes into a list
        return genotypes.groupBy({ it.first }, { it.second }).map { (cleanedAllele, sampleGametes) ->
            HvcfRangeHapIdSampleGamete(refRange, cleanedAllele, sampleGametes)
        }
    }


    /**
     * Build a TreeMap for both the start and the end point for each reference range.
     * This should allow us to easily hit the correct reference ranges.
     *
     * TODO Check to see if we can just use the start...We should be able to just makes it a bit trickier logic
     */
    fun buildPositionToRefRangeMap(ranges: Set<ReferenceRange>): TreeMap<Position, ReferenceRange> {
        //Use a TreeMap to allow for fast lookup of the position  A treeRangeMap would also work, but they tend to be slower than simple TreeMap
        //TODO implement primitive map

        val positionMap = TreeMap<Position, ReferenceRange>()
        for (range in ranges) {
            positionMap[Position(range.contig, range.start)] = range
            positionMap[Position(range.contig, range.end)] = range
        }
        return positionMap
    }

    /**
     * Function to extract out the VCF variants and export to a merged VCF file
     *
     * asmHapIdMap is a Map of RefRange + SampleId -> HapId
     * refRangeAndHapIdMap is Map of RefRange + HapId -> list<SampleGamete>
     * positionToRangeMap is a TreeMap of Position -> RefRange to allow for fast lookup of the RefRange for each position in the donor VCF file
     */
    fun extractVcfAndExport( asmHapIdMap: Map<Pair<ReferenceRange,String>, List<HvcfVariant>>,
                             refRangeAndHapIdMap: Map<Pair<ReferenceRange, String>, List<SampleGamete>>,
                             positionToRangeMap: TreeMap<Position, ReferenceRange>,
                             donorVcfFile: String,
                             outputFileName: String,
                             refSeq: Map<String, NucSeq>) {

        //Get out the sampleNames
        val sampleNamesSorted = refRangeAndHapIdMap.values
            .flatten()
            .map { it.name }
            .toSet()
            .sorted()

        //Build an output VCF With the header and sampleNames coming from the outputSampleNames
        val writer = VariantContextWriterBuilder()
            .unsetOption(Options.INDEX_ON_THE_FLY)
            .setOutputFile(File(outputFileName))
            .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
            .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
            .build()

        val header = createGenericHeader(sampleNamesSorted,setOf())
        addSequenceDictionary(header, refSeq)
        writer.writeHeader(header)

        extractVariantsAndWrite(asmHapIdMap, refRangeAndHapIdMap, positionToRangeMap, donorVcfFile, writer)

        writer.close()
    }

    /**
     * function to extract the variants out and write to the outputWriter
     */
    fun extractVariantsAndWrite(asmHapIdMap: Map<Pair<ReferenceRange,String>, List<HvcfVariant>>,
                                refRangeAndHapIdMap: Map<Pair<ReferenceRange, String>, List<SampleGamete>>,
                                positionToRangeMap: TreeMap<Position, ReferenceRange>,
                                donorVcfFile: String,
                                outputWriter: VariantContextWriter) {
        //Then we can walk through the mergedVCF file
        //For each position, lookup the RefRange,
        //Then lookup RefRange+ASMName in asmHapIdMap to get out the HapId for each sampleName in the VCF
        //Then for each RefRange+ASMNames hapIds we lookup the SampleName and Gamete and write it to that genotype

        val donorVCFReader = VCFFileReader(File(donorVcfFile), false)

        donorVCFReader.forEach { context ->
            val newContext = buildVariantContext(context, positionToRangeMap, asmHapIdMap, refRangeAndHapIdMap) ?: return@forEach
            outputWriter.add(newContext)
        }
    }

    /**
     *  Function to build a new variant context for the output VCF file.
     *  This is where we pull out the correct genotypes for each sample based on the haplotype information from the hvcf files and the donor VCF file.
     */
    fun buildVariantContext(
        vcfContext: VariantContext,
        positionToRangeMap: TreeMap<Position, ReferenceRange>,
        asmHapIdMap: Map<Pair<ReferenceRange, String>, List<HvcfVariant>>,
        refRangeAndHapIdMap: Map<Pair<ReferenceRange, String>, List<SampleGamete>>
    ): VariantContext? {
        val position = Position(vcfContext.contig, vcfContext.start)
        val refRange = positionToRangeMap.floorEntry(position)?.value ?: return null

        if(refRange.contig != vcfContext.contig || refRange.start > vcfContext.start || refRange.end < vcfContext.end) {
            return null
        }

        //Walk through the gvcfGenotypeAlleles
        val newContextBuilder = VariantContextBuilder()
            .chr(vcfContext.contig)
            .start(vcfContext.start.toLong())
            .stop(vcfContext.end.toLong())
            .alleles(vcfContext.alleles)
            .id(vcfContext.id)
            .attributes(vcfContext.attributes)

        //he go through the genotypes those sample names should match the sample names needed from the asmHapIdMap
        val sampleGameteAndAllelePairs = extractAllelesForEachSampleGamete(vcfContext, asmHapIdMap, refRange, refRangeAndHapIdMap)

        //Now we need to group them by the sample names
        val outputGenotypes = buildOutputGenotypes(sampleGameteAndAllelePairs)

        //Add the genotypes to the variant
        newContextBuilder.genotypes(outputGenotypes)
        return newContextBuilder.make()
    }

    /**
     * Function to build a list of output Genotypes based on the sampleGameteAndAllelePairs
     */
    fun buildOutputGenotypes(sampleGameteAndAllelePairs: List<Pair<SampleGamete, Allele>>): List<Genotype> {
        return sampleGameteAndAllelePairs.groupBy { it.first.name }
            .map { (sampleName, listOfAlleles) ->
                val alleles = listOfAlleles.sortedBy { it.first.gameteId }
                    .map { sampleGameteAndAllele -> sampleGameteAndAllele.second }

                val genotype = GenotypeBuilder()
                    .name(sampleName)
                    .alleles(alleles)
                    .make()

                genotype
        }
    }

    /**
     * Function for extracting out the alleles for each SampleGamete and output the list for output
     * Context here is a vcf variant
     */
    fun extractAllelesForEachSampleGamete(
        vcfContext: VariantContext,
        asmHapIdMap: Map<Pair<ReferenceRange, String>, List<HvcfVariant>>,
        refRange: ReferenceRange,
        refRangeAndHapIdMap: Map<Pair<ReferenceRange, String>, List<SampleGamete>>
    ): List<Pair<SampleGamete, Allele>> {

        return vcfContext.genotypes.flatMap { genotype ->
            val sampleName = genotype.sampleName
            check(
                asmHapIdMap.containsKey(
                    Pair(
                        refRange,
                        sampleName
                    )
                )
            ) { "Unable to find donor hapId for $refRange, sampleName: $sampleName" }
            val donorHapId = asmHapIdMap[Pair(refRange, sampleName)]!!.first().hapId

            // We actually do not need to a normal check here.
            // IF a haplotype is not hit in the imputed it will not be in the refRangeAndHapIdMap so we can skip
            if(!refRangeAndHapIdMap.containsKey(Pair(refRange, donorHapId))) {
                return@flatMap emptyList()
            }

            val sampleGametesForThisHapId = refRangeAndHapIdMap[Pair(refRange, donorHapId)]!!
            sampleGametesForThisHapId.map { sampleGamete ->
                Pair(sampleGamete, genotype.getAllele(0)) // This should work correctly as the donor VCF is haploid
            }
        }
    }

}