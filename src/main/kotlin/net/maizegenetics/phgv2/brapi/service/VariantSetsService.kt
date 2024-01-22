package net.maizegenetics.phgv2.brapi.service

import net.maizegenetics.phgv2.brapi.model.DataFormatEnum
import net.maizegenetics.phgv2.brapi.model.FileFormatEnum
import net.maizegenetics.phgv2.brapi.model.VariantSet
import net.maizegenetics.phgv2.brapi.model.VariantSetAvailableFormats
import org.apache.logging.log4j.LogManager
import java.io.File

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
        val variantSetURI = "${tildbUri}/variantsets/allSamplesMerged.vcf.gz"
        val availableFormats = VariantSetAvailableFormats(DataFormatEnum.VCF, fileFormat = FileFormatEnum.APPLICATION_ZIP, fileURL=variantSetURI)
        if (!File(variantSetURI).exists()) {
            // create the file using Terry's code.
        }
        // These values may not be correct, but the availableFormat.fileURI should be correct
        // What is our studDbId now?
        // Is the callSetCount the number of samples?  We can probably get that by querying the VCF file, or the
        // AGC or tiledb file.  This assumes we made the allSamplesMerged.vcf.gz file from what is currently in the db.
        return VariantSet(callSetCount = 0, referenceSetDbId = "all", studyDbId = "all", variantCount = 0, variantSetDbId = "all", variantSetName = "all", availableFormats = listOf(availableFormats))
    }
}