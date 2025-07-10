package net.maizegenetics.phgv2.cli

import biokotlin.util.*
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFHeader
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import net.maizegenetics.phgv2.utils.MASTER_GVCF_NAME
import net.maizegenetics.phgv2.utils.exportGvcfFilesAndMerge
import org.apache.logging.log4j.LogManager
import java.io.File

class CreateMasterGvcf : CliktCommand(help = "Create master GVCF file") {

    private val myLogger = LogManager.getLogger(CreateMasterGvcf::class.java)

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        .default("")

    val hvcfFile by option()
        .default("")

    val outputDir by option(help = "Output directory for the gVCF files.  If not provided, the current working directory is used.")
        .default("")

    override fun run() {
        logCommand(this)
        exportGvcfFilesAndMerge(dbPath)
        val gvcfFile = "$outputDir/${File(hvcfFile).name.substringBeforeLast(".h.vcf")}.g.vcf"
        convertHvcfToGvcf(dbPath, hvcfFile, gvcfFile)
    }

    fun convertHvcfToGvcf(dbPath: String, hvcfFile: String, gvcfFile: String) {

        runBlocking {

            val masterGVCF = "$dbPath/$MASTER_GVCF_NAME"
            val masterGvcfReader = vcfReader(masterGVCF, debug = false)

            val hvcfReader = vcfReader(hvcfFile, debug = false)
            val hvcfSamples = hvcfReader.samples
            val hvcfAltHeaders: Map<String, AltHeaderMetaData> = hvcfReader.altHeaders

            val variantContextChannel = Channel<Deferred<List<VariantContext>>>(100)

            launch(Dispatchers.IO) {

                var numHvcfRecords = 0
                while (true) {

                    val hvcfVariant = hvcfReader.variant() ?: break

                    numHvcfRecords++

                    val masterGvcfVariants = overlappingVariants(hvcfVariant, masterGvcfReader)
                    variantContextChannel.send(async {
                        convertHvcfVariantToGvcfVariants(
                            hvcfSampleToGvcfSample(hvcfAltHeaders, hvcfVariant),
                            hvcfVariant,
                            masterGvcfVariants
                        )
                    })

                    hvcfReader.advanceVariant()

                }

                println("Processed $numHvcfRecords h.vcf records")

                variantContextChannel.close()

            }

            writeOutputVCF(gvcfFile, hvcfSamples, variantContextChannel)

        }


    }

    /**
     * Convert the hvcf variant to gvcf variants.
     */
    private fun convertHvcfVariantToGvcfVariants(
        hvcfSampleToGvcfSample: Map<String, String?>,
        hvcfVariant: SimpleVariant,
        masterGvcfVariants: List<SimpleVariant>
    ): List<VariantContext> {
        return masterGvcfVariants
            .mapNotNull { variant ->
                if (refRangeContainsVariant(hvcfVariant, variant)) {
                    convertHvcfVariantToGvcfVariant(hvcfSampleToGvcfSample, hvcfVariant, variant)
                } else {
                    // TODO: Handle the case where the hvcf variant is not in the range of the master gvcf variant
                    null
                }
            }
    }

    private fun refRangeContainsVariant(refRangeVariant: SimpleVariant, variant: SimpleVariant): Boolean {
        return refRangeVariant.contig == variant.contig &&
                refRangeVariant.start <= variant.start &&
                refRangeVariant.end >= variant.end
    }

    private fun convertHvcfVariantToGvcfVariant(
        hvcfSampleToGvcfSample: Map<String, String?>,
        hvcfVariant: SimpleVariant,
        masterGvcfVariant: SimpleVariant
    ): VariantContext? {

        if (masterGvcfVariant.refAllele == "N") return null

        val refAllele = alleleRef(masterGvcfVariant.refAllele)

        val alleleMap = mutableMapOf<String, Allele>()
        masterGvcfVariant.altAlleles.forEach { alleleMap[it] = alleleAlt(it) }
        // reference added last to overwrite a matching alt allele
        alleleMap[masterGvcfVariant.refAllele] = refAllele

        val samples = hvcfVariant.samples

        val genotypes = samples
            .map { sample ->

                val gvcfSample = hvcfSampleToGvcfSample[sample]
                val alleles = if (gvcfSample != null) {
                    masterGvcfVariant.genotypeStrs(gvcfSample)
                } else {
                    listOf(".")
                }

                val alleleObjs = alleles
                    .map { allele ->
                        when (allele) {
                            "." -> Allele.NO_CALL

                            "<INS>" -> Allele.SV_SIMPLE_INS

                            "<DEL>" -> Allele.SV_SIMPLE_DEL

                            "REF" -> refAllele

                            else -> alleleMap[allele]
                                ?: throw IllegalArgumentException("Allele not found: $allele at contig: ${masterGvcfVariant.contig} position: ${masterGvcfVariant.positionRange} sample: $sample alleleMap: $alleleMap")
                        }
                    }

                val isPhased = if (gvcfSample != null) {
                    masterGvcfVariant.isPhased(gvcfSample)
                } else {
                    true
                }

                GenotypeBuilder(sample, alleleObjs)
                    .phased(isPhased)
                    .make()

            }

        val builder = VariantContextBuilder()
            .source(".")
            .alleles(alleleMap.values)
            .chr(masterGvcfVariant.contig)
            .start(masterGvcfVariant.start.toLong())
            .stop(masterGvcfVariant.end.toLong())
            .genotypes(genotypes)

        if (masterGvcfVariant.end != masterGvcfVariant.start) {
            builder.attribute("END", masterGvcfVariant.end)
        }

        return builder.make()

    }

    /**
     * Creates a HTSJDK reference allele.
     */
    private fun alleleRef(allele: String): Allele {
        return Allele.create(allele, true)
    }

    /**
     * Creates a HTSJDK alternate allele.
     */
    private fun alleleAlt(allele: String): Allele {
        return Allele.create(allele, false)
    }

    private suspend fun overlappingVariants(
        hvcfVariant: SimpleVariant,
        masterGvcfReader: VCFReader
    ): List<SimpleVariant> {

        val contig = hvcfVariant.contig
        val start = hvcfVariant.start
        val end = hvcfVariant.end

        var variant = masterGvcfReader.variant()

        // Advance to the first variant that is on the same contig
        while (variant != null && variant.contig != contig) {
            masterGvcfReader.advanceVariant()
            variant = masterGvcfReader.variant()
        }

        // Advance to the first variant that is within the range of the hvcf variant
        while (variant != null && variant.contig == contig && variant.end < start) {
            masterGvcfReader.advanceVariant()
            variant = masterGvcfReader.variant()
        }

        val overlappingVariants = mutableListOf<SimpleVariant>()
        while (variant != null && variant.contig == contig && variant.start <= end) {

            overlappingVariants.add(variant)

            if (variant.end <= end) {
                masterGvcfReader.advanceVariant()
                variant = masterGvcfReader.variant()
            } else {
                break
            }

        }

        return overlappingVariants

    }

    private fun hvcfSampleToGvcfSample(
        altHeaders: Map<String, AltHeaderMetaData>,
        variant: SimpleVariant
    ): Map<String, String?> {
        return variant.samples
            .associateWith { sample ->
                val hapid = variant.genotypeStrs(sample).first().substringAfter("<").substringBeforeLast(">")
                if (hapid != ".") {
                    altHeaders[hapid]?.sampleName()
                        ?: throw IllegalStateException("Hapid: $hapid for sample: $sample not found in altHeaders")
                } else {
                    null
                }

            }
    }

    private suspend fun writeOutputVCF(
        outputFile: String,
        samples: List<String>,
        variantContextChannel: Channel<Deferred<List<VariantContext>>>
    ) {

        // Write the output GVCF file, using the HTSJDK VariantContextWriterBuilder
        VariantContextWriterBuilder()
            .unsetOption(Options.INDEX_ON_THE_FLY)
            .setOutputFile(File(outputFile))
            .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
            .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
            .build()
            .use { writer ->

                val header = VCFHeader(createGenericVCFHeaders(samples))
                writer.writeHeader(header)

                var numVariantsWritten = 0
                var numDeferredProcessed = 0
                for (deferred in variantContextChannel) {
                    val variantContexts = deferred.await()
                    numDeferredProcessed++
                    variantContexts.forEach { variantContext ->
                        numVariantsWritten++
                        writer.add(variantContext)
                    }
                }
                println("Processed $numDeferredProcessed deferreds")
                println("Wrote $numVariantsWritten variants to $outputFile")

            }

    }

}