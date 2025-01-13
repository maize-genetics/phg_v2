package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.utils.*
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * Class to load the HVCF file into a core tiledb array
 *
 * This will load to the core tiledb array the HVCF files that are in the directory
 * This array will live in the same folder that contains the vcf dataset .
 *
 * Both hvcfs from aligned assemblies and from imputation may be loaded from the same directory.
 * Imputed hvcfs will not have alt headers but that is handled in the code below
 *
 * In the future: Initdb should create both the tiledbvcf array for the gvcf files, and the tiledbcore array for the hvcf files
 * For now, we call InitHvcfArray to create the core array for the hvcf files
 * before calling this function.
 */

class LoadHvcf: CliktCommand(help = "Load  h.vcf files into TileDB core datasets") {
    private val myLogger = LogManager.getLogger(LoadHvcf::class.java)

    val hvcfDir by option(help = "Full path to HVCF file directory")
        .required()

    val dbPath by option(help = "Folder holding TileDB datasets")
        .default("")

    // Pre-compile the Regex pattern - used when creating the output fasta file names
    val HVCF_PATTERN = Regex("""(\.hvcf|\.h\.vcf|\.hvcf\.gz|\.h\.vcf\.gz)$""")

    override fun run() {
        logCommand(this)

        val dbPath = if (dbPath.isBlank()) {
            System.getProperty("user.dir")
        } else {
            dbPath
        }
        // Verify the tiledbURI - an exception is thrown from verifyURI if the URI is not valid
        myLogger.info("LoadHvcf: verifying array")
        val goodArray = TileDBCoreHvcfUtils.verifyHvcfArray(dbPath)
        myLogger.info("LoadHvcf: goodArray: $goodArray")

        processHvcfsToTiledbArrays(hvcfDir, dbPath)

    }

    // function to write the altheader  and context data to the tiledb
    // alt_header_array and hvcf_variants_array respectively.
    fun processHvcfsToTiledbArrays(hvcfDir:String, dbPath:String) {
        // Load the arrays
        // Must first parse the hvcfs into the proper format
        // Then call function to load.

        // do I want to parse and load in pairs, or parse all to a single map and then load all at once?
        // Could be memory issue.  But ... I'm not pulling sequence here - not like when creating the hvcfs
        // when we had to pull sequence to create the hapid/sequence-hash.

        val combinedHvcfHeaderData = mutableListOf<Map<String, String>>() // list of maps holding combined parsed hvcf header data
        val combinedHvcfVariantData = mutableListOf<Map<String, String>>() // list of maps holding combined parsed hvcf context data

        myLogger.info("processHvcfsToTIledbArrays: hvcfDir: $hvcfDir, dbPath: $dbPath printing file list:")
        val hvcfFiles = File(hvcfDir).listFiles { file ->
            myLogger.info("file: ${file.name}")
            HVCF_PATTERN.containsMatchIn(file.name)
        }
        if (hvcfFiles == null || hvcfFiles.isEmpty()) {
            myLogger.error("No hvcf files found in $hvcfDir")
            return
        }
        myLogger.info("hvcfFiles size: ${hvcfFiles.size}, begin parsing hvcf files.")
        val altArrayName = "${dbPath}/alt_header_array"
        val variantArrayName = "${dbPath}/hvcf_variants_array"
        hvcfFiles.forEach{ file ->
            // data is written for each file as it is processed.
            // This is to ensure batching for tiledb arrays.
            val vcfReader = VCFFileReader(file, false)
            val hvcfData = TileDBCoreHvcfUtils.parseTiledbAltHeaders(vcfReader)
            myLogger.info("hvcfData size: ${hvcfData.size} for file ${file.name}")
            if (hvcfData.isNotEmpty()) {
                // hvcf files from imputation may not have alt headers
                TileDBCoreHvcfUtils.writeAltDataToTileDB(altArrayName, hvcfData)
                // combinedHvcfHeaderData.addAll(hvcfData)
                myLogger.info("Finished writing ALT data for file ${file.name}")
            }

            // Then parse the body of this file
            // We want to store the SampleName, plus ID1 and ID2 to the tiledbArray
            // ID1 and ID2 will be the same if this is haploid data

            // TODO:  need to handle multipsample hvcf files
            // WIll we have them from imputation or are these always single sample unless someone merges them?
            val variantData = TileDBCoreHvcfUtils.parseTiledbVariantData(vcfReader)
            if (variantData.isEmpty()) {
                myLogger.warn("No variant data found in ${file.name}")
                return
            }
            TileDBCoreHvcfUtils.writeVariantsDataToTileDB(variantArrayName, variantData)
            myLogger.info("Finished writing Variant data for file ${file.name}")
            //combinedHvcfVariantData.addAll(variantData)

        }

    }
}