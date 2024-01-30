package net.maizegenetics.phgv2.brapi.service

import net.maizegenetics.phgv2.brapi.model.DataFormatEnum
import net.maizegenetics.phgv2.brapi.model.FileFormatEnum
import net.maizegenetics.phgv2.brapi.model.VariantSet
import net.maizegenetics.phgv2.brapi.model.VariantSetAvailableFormats
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * This class is used to generate the VariantSets object that is returned by the "variantsets" endpoint.
 * It looks for a vcf file called allSamplesMerged.vcf.gz in the tiledb_uri directory.
 * If it doesn't find it , one is created using the ExportVCF plugin.
 * A URL pointing to the file is returned in the "availbleFormats" section of the VariantSets object.
 */
class VariantSetsService {

    private val myLogger = LogManager.getLogger(VariantSetsService::class.java)

    init{
        myLogger.info("VariantSetsService initialized.")
        // anything to go in here?
    }

    // Should this return a VarianatSet?
    fun generateVariantSets(tildbUri:String): VariantSet {
        // LCJ - I made up this name.  What should it be called?
        // will there be more than 1 variantset?  I named it "all"
        // as that is what we called the reference cache.  What is a better name?
        val variantSetURI = "${tildbUri}/variantsets/allSamplesMerged.h.vcf.gz"
        val availableFormats = VariantSetAvailableFormats(DataFormatEnum.VCF, fileFormat = FileFormatEnum.TEXT_TSV, fileURL=variantSetURI)
        if (!File(variantSetURI).exists()) {
            // TODO create and export the file using Terry's code.
            // First create the graph, then call ExportVCF
            // Are the hvcf files currently written to tiledb_uri/hvcf_files when they are created?
            // Or do we need to export the hvcf files to tiledb_uri/hvcf_files and then create the graph
            // from them?
        }
        // These values may not be correct, but the availableFormat.fileURI should be correct
        // What is our studDbId now?
        // Is the callSetCount the number of samples?  We can probably get that by querying the VCF file, or the
        // AGC or tiledb file.  This assumes we made the allSamplesMerged.vcf.gz file from what is currently in the db.
        return VariantSet(callSetCount = 0, referenceSetDbId = "all", studyDbId = "all", variantCount = 0, variantSetDbId = "all", variantSetName = "all", availableFormats = listOf(availableFormats))
    }
}