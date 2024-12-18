package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import io.tiledb.java.api.Array
import io.tiledb.java.api.QueryType
import io.tiledb.java.api.*
import net.maizegenetics.phgv2.utils.createTileDBCoreArrays
import net.maizegenetics.phgv2.utils.parseTiledbAltHeaders
import net.maizegenetics.phgv2.utils.writeAltDataToTileDB
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

        processAltHeaderToArray(hvcfDir, dbPath)

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

    // function to write the altheader data to the tiledb array
    fun processAltHeaderToArray(hvcfDir:String, dbPath:String) {
        // Load the arrays
        // Must first parse the hvcfs into the proper format
        // Then call function to load.

        // do I want to parse and load in pairs, or parse all to a single map and then load all at once?
        // Could be memory issue.  But ... I'm not pulling sequence here - not like when creating the hvcfs
        // when we had to pull sequence to create the hapid/sequence-hash.

        val combinedHvcfData = mutableListOf<Map<String, String>>() // list of maps holding combined parsed hvcf data
        //val hvcfFiles = File(hvcfDir).listFiles().filter{it.endsWith("h.vcf")}.toList()

        val hvcfFiles = File(hvcfDir).listFiles { file ->
            HVCF_PATTERN.containsMatchIn(file.name)
        }
        myLogger.info("hvcfFiles size: ${hvcfFiles.size}, begin parsing hvcf files.")
        hvcfFiles.forEach{ file ->
            val hvcfData = parseTiledbAltHeaders(file.toString())
            combinedHvcfData.addAll(hvcfData)
        }

        // Now load the data
        myLogger.info("Before writing to tiledb: size if combinedHvcfData: ${combinedHvcfData.size}")
        val arrayName = "${dbPath}/alt_header_array"
        writeAltDataToTileDB(arrayName, combinedHvcfData)
    }
}