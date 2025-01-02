package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import htsjdk.variant.vcf.VCFFileReader
import io.tiledb.java.api.Array
import io.tiledb.java.api.QueryType
import io.tiledb.java.api.*
import net.maizegenetics.phgv2.utils.*
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * Class to load the HVCF file into a core tiledb array
 *
 * This will load to the core tiledb array the HVCF files that are in the directory
 * This array will live in the same folder that contains the vcf dataset .
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

    // threads is not required.  Defaults to 1
    // While technically, this should be an "int", it is passed to ProcessBuilder as a string
    // so we're leaving it defined as a string here.
    val threads by option(help = "Number of threads for use by tiledb")
        .default("1")

    val type by option(help = "Type of hvcf file  to process: choices are assembly or imputed")
        .default("assembly")

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
        println("LoadHvcf: verifying array")
        val goodArray = verifyHvcfArray(dbPath)
        println("LoadHvcf: goodArray: $goodArray")

        processHvcfsToTiledbArrays(hvcfDir, dbPath)

    }
    // Hvcf array must exist in the dbPath folder
    fun verifyHvcfArray(dbPath:String): Boolean {
        // Verify the folder exists and the hvcf_array exists in that folder
        // as a tiledb array
        val arrayName = "${dbPath}/alt_header_array"
        if (!File(dbPath).exists() || !File(dbPath).isDirectory()) {
            myLogger.warn("Folder $dbPath does not exist - creating.")
            File(dbPath).mkdirs()
            createTileDBCoreArrays(arrayName)
        } else {
            myLogger.info("Folder $dbPath exists, trying to open array name $arrayName")
            try {
                val array = Array(Context(), arrayName, QueryType.TILEDB_READ)
            } catch (exc: TileDBError) {
                myLogger.error("Array $arrayName does not exist or exists but is not a tiledb array.\n" +
                        "Please check your arrayName for accuracy and run  Initdb to set up the arrays if necessary.")
                throw exc
            }
        }
        return true

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
        hvcfFiles.forEach{ file ->
            val vcfReader = VCFFileReader(file, false)
            val hvcfData = parseTiledbAltHeaders(vcfReader)
            combinedHvcfHeaderData.addAll(hvcfData)

            // Then parse the body of this file
            // We want to store the sampleGamete, the ID and the GT value to the tiledbArray
            // so the datatype should be another List<Map<String,String>> where the map is
            // sampleGamete, ID
            //val haplotypeVariants = vcfReader.iterator().asSequence().toList()

            // TODO:  need to handle multipsample hvcf files
            // WIll we have them from imputation or are these always single sample unless someone merges them?
            val sampleName = vcfReader.header.genotypeSamples.firstOrNull()
            // Iterate through VariantContext objects in the VCF
            vcfReader.iterator().forEach { variantContext ->
                val genotype = variantContext.getGenotype(sampleName)

                if (genotype.isNoCall) {
                    println("No GT call for variant at ${variantContext.contig}:${variantContext.start}")
                    return@forEach
                }

                val altAlleles = variantContext.alternateAlleles
                val gtField = genotype.genotypeString // need this to determine how many alleles there are
                val alleles = gtField.split("[/|]".toRegex()) // Split on "/" or "|" for diploid

                // create the refRange mapping
                val chr = variantContext.contig
                val start = variantContext.start
                val end = variantContext.end
                val refRange = "${chr}:${start}-${end}"

                // Create the entry.  Only 1 if is haploid, 2 if is diploid
                val entry = mutableMapOf<String, String>()
                entry["RefRange"] = refRange
                entry["ID1"] = altAlleles[0].displayString
                entry["SampleName"] = "${sampleName}"
                entry["ID2"] = altAlleles[0].displayString

                // DO we want this to be 0, or do we want this to be
                // the same value as the ID1 field?
                // If there is only 1 allele, then the ID2 field should be the same as the ID1 field
                // when it is diploid.  But what about haploid - how do we know it is hapoid if we
                // fill in the ID2 field?
                if (alleles.size == 2 && altAlleles.size == 2) {
                    entry["ID2"] = altAlleles[1].displayString // If there are 2 alt alleles, set ID2 to the second allele
                }
                combinedHvcfVariantData.add(entry)
            }
        }

        // Now load the data
        myLogger.info("Writing to tiledb: size of combinedHvcfALTHeaderData: ${combinedHvcfHeaderData.size}")
        val arrayName = "${dbPath}/alt_header_array"
        writeAltDataToTileDB(arrayName, combinedHvcfHeaderData)
        myLogger.info("Finished writing ALT data")
        myLogger.info("\nWriting to tiledb: size of combinedHvcfVariantData: ${combinedHvcfVariantData.size}")
        val variantArrayName = "${dbPath}/hvcf_variants_array"
        writeVariantsDataToTileDB(variantArrayName, combinedHvcfVariantData)
        myLogger.info("Finished writing Variant data")
    }
}