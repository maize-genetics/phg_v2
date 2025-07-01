package net.maizegenetics.phgv2.brapi.service

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.exportMultiSampleHVCF
import net.maizegenetics.phgv2.brapi.model.DataFormatEnum
import net.maizegenetics.phgv2.brapi.model.FileFormatEnum
import net.maizegenetics.phgv2.brapi.model.VariantSet
import net.maizegenetics.phgv2.brapi.model.VariantSetAvailableFormats
import net.maizegenetics.phgv2.brapi.utilities.BrAPIConfig
import net.maizegenetics.phgv2.cli.ExportVcf
import net.maizegenetics.phgv2.cli.StartServer
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * This class is used to generate the VariantSets object that is returned by the "variantsets" endpoint.
 * It looks for a vcf file called allSamplesMerged.vcf.gz in the tiledb_uri directory.
 * If it doesn't find it , one is created using the ExportVCF plugin.
 * A URL pointing to the file is returned in the "availableFormats" section of the VariantSets object.
 */
object VariantSetsService {

    private val myLogger = LogManager.getLogger(VariantSetsService::class.java)

    private val individualSamplesDir = "${BrAPIConfig.tiledbURI}/individualSamples/"
    private val variantSetsDir = "${BrAPIConfig.tiledbURI}/variantsets/"
    const val allSamplesFileName = "allSamplesMerged.h.vcf"
    val allSamplesHvcf = "${variantSetsDir}$allSamplesFileName"

    init {
        File(individualSamplesDir).mkdirs()
        File(variantSetsDir).mkdirs()
    }

    suspend fun getVariantSet(): VariantSet {

        val availableFormats = VariantSetAvailableFormats(DataFormatEnum.VCF, fileFormat = FileFormatEnum.TEXT_TSV, fileURL = "${StartServer.serverURL()}$allSamplesFileName")

        // These values may not be correct, but the availableFormat.fileURI should be correct
        // What is our studDbId now?
        // Is the callSetCount the number of samples?  We can probably get that by querying the VCF file, or the
        // AGC or tiledb file.  This assumes we made the allSamplesMerged.vcf.gz file from what is currently in the db.
        return VariantSet(callSetCount = 0, referenceSetDbId = "all", studyDbId = "all", variantCount = 0, variantSetDbId = "all", variantSetName = "all", availableFormats = listOf(availableFormats))

    }

    /**
     * This function creates a single hvcf file for each sample in the database.
     */
    private fun createSingleSampleHVCFs() {

        val sampleNames = SamplesService.allTaxaNames().joinToString(",") { it.sampleName }

        val exportVcfCommand = "--db-path ${BrAPIConfig.tiledbURI} --sample-names $sampleNames -o $individualSamplesDir"
        myLogger.info("createSingleSampleHVCFs: ExportVcf Command: $exportVcfCommand")
        ExportVcf().test(exportVcfCommand)

    }

    /**
     * This function creates a single hvcf file from all the individual hvcf files.
     * Returns the name of the file created.
     */
    fun createAllSamplesHVCF(): String {

        createSingleSampleHVCFs()

        val graph = HaplotypeGraph(individualSamplesDir)

        exportMultiSampleHVCF(graph, allSamplesHvcf)

        return allSamplesHvcf

    }

}