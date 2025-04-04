package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import net.maizegenetics.phgv2.utils.TileDBCoreHvcfUtils
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * This function creates the tiledb array for the hVCF files
 * It uses the Java API to create the array
 * The schema used is the original schema:  RefRanges and SampleNames
 * as dimensions, RefChecksum, Regions and IDs as the attributes.
 * While the schema may change, examples of how to use the Java API for
 * tiledb core arrays should remain constant.
 */
class InitHvcfArray: CliktCommand(help = "Create TileDB core array  h.vcf files") {
    private val myLogger = LogManager.getLogger(InitHvcfArray::class.java)

    val dbPath by option(help = "Folder name under which TileDB datasets will be created. If this folder does not exist, it will be created.")
        .default("")

    override fun run() {
        logCommand(this)
        // Set the dbPath to the current directory if it is not set
        val dbPath = if (dbPath.isBlank()) {
            System.getProperty("user.dir")
        } else {
            dbPath
        }

        // if dbPath does not exist as a directly, throw an exception
        val dbfolder = File(dbPath)
        if (!dbfolder.exists()) {
            myLogger.warn("Folder $dbPath does not exist .")
            dbfolder.mkdirs()
        }

        // Create the tiledb core array
        TileDBCoreHvcfUtils.createTileDBCoreArrays(dbPath)
        myLogger.info("TileDB core array created for h.vcf files in $dbPath")
    }
}