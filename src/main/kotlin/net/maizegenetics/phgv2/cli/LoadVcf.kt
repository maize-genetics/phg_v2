package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.utils.inputStreamProcessing
import net.maizegenetics.phgv2.utils.verifyURI
import org.apache.logging.log4j.LogManager
import java.io.BufferedInputStream
import java.io.File
import java.lang.Exception
import java.nio.file.Files
import java.nio.file.Paths
import java.util.stream.Collectors

/**
 * This class will load gvcf and hvcf files into TileDB datasets.
 * It has these steps:
 *  1. create list of gvcf vs hvcf files for loading
 *      a. if no hvcf.gz, h.vcf.gz, gvcf.gz, g.vcf.gz files are found, return error
 *  2. verify the dbpath exists and TileDB datasets exists in that folder
 *      a. if no datasets exist, create them
 *  3. load the gvcf and hvcf files
 *  4. when loading the hvcf files, make a copy of the file in the db-path/hvcf_files folder
 *     This will be used by the brAPI server to create the initial graph used to service
 *     the brAPI calls.
 */
class LoadVcf : CliktCommand(help = "Load g.vcf and h.vcf files into TileDB datasets") {

    private val myLogger = LogManager.getLogger(LoadVcf::class.java)

    val vcfDir by option(help = "Full path to VCF file directory")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--vcf-dir must not be blank"
            }
        }

    val dbPath by option(help = "Folder holding TileDB datasets")
        .default("")

    // threads is not required.  Defaults to 1
    // While technically, this should be an "int", it is passed to ProcessBuilder as a string
    // so we're leaving it defined as a string here.
    val threads by option(help = "Number of threads for use by tiledb")
        .default("1")

    val condaEnvPrefix by option (help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
        .default("")

    override fun run() {

        val dbPath = if (dbPath.isBlank()) {
            System.getProperty("user.dir")
        } else {
            dbPath
        }
        // Verify the tiledbURI - an exception is thrown from verifyURI if the URI is not valid
        val validDB = verifyURI(dbPath,"hvcf_dataset")

        loadVcfFiles(vcfDir,dbPath,condaEnvPrefix,threads)
    }

    fun loadVcfFiles(vcfDir:String,dbPath:String,condaEnvPrefix:String,threads:String="1") {

        // Check the type of files in the vcfDir
        // anything with h.vcf.gz or hvcf.gz is a hvcf file
        // anything with g.vcf.gz or gvcf.gz is a gvcf file
        // These are compressed files as tiledb needs the bgzipped and indexed files
        val fileLists = getFileLists(vcfDir)

        if (fileLists.first.isEmpty() && fileLists.second.isEmpty()) {
            // no gvcf or hvcf files found in the user supplied vcfDir
            myLogger.warn("No files ending in g.vcf.gz or h.vcf.gz found in $vcfDir.  Note that both the bgzipped and indexed files must exist in the specified folder. \nPlease check the folder and try again.")
            throw IllegalArgumentException("LoadVCF: No files ending in g.vcf.gz or h.vcf.gz found in $vcfDir.  Note that both the bgzipped and indexed files must exist in the specified folder. \nPlease check the folder and try again.")
        }

        // verify the URI is correct,
        if (fileLists.first.isNotEmpty()) {
            // look for duplicate sample names - tiledb uri was verified at start of run()
            // Get samples from tiledb
            val uri = dbPath + "/gvcf_dataset"
            val tiledbSampleList = getTileDBSampleLists(uri)
            // verify the samples in the vcf files are not already in the tiledb dataset
            val vcfSampleList = getVcfSampleLists(fileLists.first)
            // verify there are no overlaps between the tiledbSampleList and the vcfSampleList
            val overlap = tiledbSampleList.intersect(vcfSampleList.keys)
            // verify there are no overlaps between the tiledbSampleList and the vcfSampleList
            if (overlap.isNotEmpty()) {
                myLogger.warn("The following samples are already in the gvcf tiledb dataset: ${overlap.joinToString(",")}")
                myLogger.warn("Please remove the vcf files containing these samples and try again.")
            }
            // Load the gvcf files!
            myLogger.info("No overlap between tiledb and vcf files.  Loading gvcf files.")
            val datasetURI = "${dbPath}/gvcf_dataset"
            loadVCFToTiledb(fileLists.first, dbPath,"gvcf_dataset", threads,condaEnvPrefix)
        }
        if (fileLists.second.isNotEmpty()) {
            // look for duplicate sample names - tiledb uri was verified at start of run()
            // Get samples from tiledb
            val uri = dbPath + "/hvcf_dataset"
            val tiledbSampleList = getTileDBSampleLists(uri)
            // verify the samples in the vcf files are not already in the tiledb dataset
            val vcfSampleList = getVcfSampleLists(fileLists.second)
            // verify there are no overlaps between the tiledbSampleList and the vcfSampleList
            val overlap = tiledbSampleList.intersect(vcfSampleList.keys)
            // verify there are no overlaps between the tiledbSampleList and the vcfSampleList
            if (overlap.isNotEmpty()) {
                // tiledb dataset was created with the -n/--no-duplicates option.  It will not load
                // duplicate positions, but will quietly skip them.
                // Warning message is printed here to alert the user of the duplicate sample names.
                myLogger.warn(
                    "The following samples are already in the tiledb hvcf dataset: ${
                        overlap.joinToString(
                            ","
                        )
                    }"
                )
                myLogger.warn("Tiledb will not load duplicate positions for these samples.")
            }
            // Load the hvcf files!
            loadVCFToTiledb(fileLists.second, dbPath, "hvcf_dataset", threads,condaEnvPrefix)

            // Copy the hvcf files to the hvcf_files folder, omitting any that are duplicates
            // of those already loaded
            val tiledbHvcfDir = "${dbPath}/hvcf_files"
            Files.createDirectories(Paths.get(tiledbHvcfDir)) // skips folders that already exist
            for (sample in vcfSampleList.keys) {
                val sourceFile = vcfSampleList[sample]
                val destFile = "${tiledbHvcfDir}/${sourceFile?.split("/")?.last()}"
                if (!File(destFile).exists()) {
                    myLogger.info("Copying ${sourceFile} to ${destFile}")
                    File(sourceFile).copyTo(File(destFile))
                }
            }
        }
    }
    // This function walks the files in the vcf folder and returns a pair of lists,
    // one for gvcf files and one for hvcf files.
    fun getFileLists(vcfDir:String):Pair<List<String>,List<String>> {
        // get a list of files in the vcfDir
        // return a pair of lists, one for gvcf files and one for hvcf files
        val vcfFiles = File(vcfDir).listFiles()
        val gvcfFiles = mutableListOf<String>()
        val hvcfFiles = mutableListOf<String>()
        for (file in vcfFiles) {
            if (file.name.endsWith("g.vcf.gz") || file.name.endsWith("gvcf.gz")) {
                gvcfFiles.add(file.toString())
            } else if (file.name.endsWith("h.vcf.gz") || file.name.endsWith("hvcf.gz")) {
                hvcfFiles.add(file.toString())
            }
        }
        return Pair(gvcfFiles,hvcfFiles)
    }


    // This function runs tiledbvcf list --uri <uri> and returns a list of sample names
    fun getTileDBSampleLists(uri:String ): List<String> {
        try {
            var builder = ProcessBuilder("conda","run","-n","phgv2-conda","tiledbvcf","list","--uri",uri)
                .start()
            val sampleListOut = BufferedInputStream(builder.inputStream, 5000000)
            var samples = inputStreamProcessing(sampleListOut)
            return samples

        } catch (exc: Exception) {
            myLogger.error("Error reading tiledb list output on ${uri}.")
            throw IllegalArgumentException("Error running ProcessBuilder to list samples from ${uri}: ${exc}")
        }
    }

    fun getVcfSampleLists(fileLists:List<String>):Map<String,String> {
        // get the sample names from the vcf files
        val vcfSampleToFileName = mutableMapOf<String,String>()
        for (file in fileLists) {
            // use htsjdk to read the vcf file and parse the sample names
            val vcfReader = VCFFileReader(File(file),false)
            val vcfHeader = vcfReader.fileHeader
            val samples = vcfHeader.sampleNamesInOrder
            for (sample in samples) {
                vcfSampleToFileName.put(sample,file)
            }
        }
        return vcfSampleToFileName
    }

    // Because this is called from both the LoadVcf class and the CreateRefVcf class,
    // we must explicitly pass dbPath as otherwise the parameter will not be known
    // when coming from CreateRefVcf, resulting in the error:
    //    Cannot read from option delegate before parsing command line
    fun loadVCFToTiledb(vcfList:List<String>, dbPath:String, dataSet:String, threads:String,condaEnvPrefix:String) {

        // declare tne temp folder
        // This will be used to write output files from ProcessBuilder commands
        // called elsewhere in this class.  It should have been created in the
        // verifyURI function,
        val tempDir = "${dbPath}/temp"
        val uri = "${dbPath}/${dataSet}"

        // get just last part of uri string, ie just the last name in this folder
        val uriName = uri.split("/").last()
        val vcfListFile = "${tempDir}/${uriName}_vcfList.txt"
        File(vcfListFile).writeText(vcfList.joinToString("\n"))
        // print the vcfListFile contents
        myLogger.info("vcfListFile contents: ${File(vcfListFile).readText()}")

        val command = if (condaEnvPrefix.isNotBlank()) mutableListOf("conda","run","-p",condaEnvPrefix,"tiledbvcf","store","--uri",uri,"-t",threads,"-f",vcfListFile,"--remove-sample-file")
        else mutableListOf("conda","run","-n","phgv2-conda","tiledbvcf","store","--uri",uri,"-t",threads,"-f",vcfListFile,"--remove-sample-file")
        // Store the files to tiledb
        var builder = ProcessBuilder(command)
        var redirectOutput = tempDir + "/tiledb_store_output.log"
        var redirectError = tempDir + "/tiledb_store_error.log"
        builder.redirectOutput( File(redirectOutput))
        builder.redirectError( File(redirectError))

        myLogger.info("begin Command:" + builder.command().joinToString(" "))
        println("begin Command:" + builder.command().joinToString(" "))

        try {
            var process = builder.start()
            var error = process.waitFor()
            if (error != 0) {
                myLogger.error("loadVCFToTiledb run via ProcessBuilder returned error code ${error}. See ${redirectError} and ${redirectOutput} for more details.")
                throw IllegalArgumentException("Error running ProcessBuilder to store vcfs to tiledb ${uri} array: ${error}")
            }

        } catch (exc: Exception) {
            myLogger.error("Error: loadVCFToTiledb: could not load tiledb array ${uri}.  See ${redirectError} and ${redirectOutput} for more details.")
            throw IllegalArgumentException("Error running ProcessBuilder to store vcfs to tiledb array: ${exc.message}")
        }
        myLogger.info("Finished! vcfs stored to tiledb dataset ${uri}")
    }

}