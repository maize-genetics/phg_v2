package net.maizegenetics.phgv2.cli

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.launch
import kotlinx.coroutines.runBlocking
import kotlinx.coroutines.withContext
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * This class takes a keyfile containing a list of fasta files and their sample names.
 * The fasta files should be the full path name to the fasta file.
 * The sample name is the name that identifies the genome to be loaded.  This name should
 * be consistent with the sample names used in the GVCF files for these genomes.
 *
 * Fasta files may be compressed or uncompressed.  If compressed, the extension should be .gz.
 * The output directory is the directory where the updated fasta files will be written.
 *
 * THis class will take the fasta files, read them into memory, update the idlines with "sampleName=<sampleName>"
 * The updated fasta files are written to the output directory with <sampleName>.fa as the file name.  Bgzip is no
 * longer performed on any of the files as this step precedes the align-assemblies step, and anchorwave does not
 * support compressed files.  The user can compress the files after this step if desired.
 *
 * "threads" is an optional parameter.  This value determines how many fasta files will be processed concurrently.
 * When determining the number of threads to use, keep in mind the memory capacity of the machine,
 * and the size of the fasta files.  Each fasta is read into memory, idlines updated, then re-written
 * to a new file of the same name in the user specified output directory.
 *
 */
class PrepareAssemblies : CliktCommand(help = "Annotate FASTA file Id lines with sample names, write new FASTAs with name <sampleName>.fa") {
    private val myLogger = LogManager.getLogger(PrepareAssemblies::class.java)

    val keyfile by option(help = "Tab-delimited file containing 2 columns name Fasta and SampleName.  Fasta column contains full path name for the fasta files.  SampleName contains the sample name for that assembly, e.g. B73 or CML247. ")
        .required()

    // This is the folder where the output files are written.  They will have the same
    // name as the original fasta files.
    val outputDir by option(help = "Folder where updated fastas will be written.")
        .required()

    val threads by option(help = "number of threads for running concurrently")
        .int()
        .default(1)


    data class InputChannelData(val fastaFile:String, val sampleName:String, val compressed:Boolean, val outputDir:String)

    override fun run() {
        logCommand(this)

        // create list of assemblies to align from the assemblies file")
        myLogger.info("creating assembliesList, calling createParallelAnnotatedFastas")
        // Read the keyfile, parse the fasta file names and the sampleName
        // Create a list of pairs of fasta file name and sampleName.  Add checking to validate
        // the format of the keyfile.

        val assemblies = File(keyfile).bufferedReader().readLines()
            .filter { !it.startsWith("#") }
            .map { line ->
                val parts = line.split("\t")
                if (parts.size != 2) {
                    throw IllegalArgumentException("Invalid format in line: \"$line\". Each line must contain exactly one tab between filename and samplename columns.")
                }
                Pair(parts[0], parts[1])
            }

        createParallelAnnotatedFastas(assemblies, outputDir)
    }


    // This method sets up parallel processing of the fasta files.
    // The files are placed in an input queue, worker threads are launched
    // to process the queue entries.  The number of worker threads is based on the
    // number of threads specified by the user.
    private fun createParallelAnnotatedFastas(assemblies:List<Pair<String,String>>, outputDir:String) {

        runBlocking {
            // Setup
            val inputChannel = Channel<InputChannelData>(100)
            launch {
                myLogger.info("Adding entries to the inputChannel:")
                assemblies.forEach { entry ->
                    val asmFile = entry.first.trim()
                    // Allow for compressed or non-compressed, e.g. .fa or .fa.gz as extensions
                    val sampleName = entry.second.trim()
                    val compressed = if (asmFile.endsWith(".gz")) true else false

                    myLogger.info("adding ${sampleName} to the inputChannel")
                    inputChannel.send(InputChannelData(asmFile, sampleName,  compressed, outputDir))
                }
                myLogger.info("Done adding data to the inputChannel")
                inputChannel.close()

            }

            val workerThreads = (1..threads).map {run ->
                launch { annotateFasta(inputChannel) }
            }
        }
    }

    // This function takes a list of fasta files and to each idline in each fasta file,
    // appends " sampleName=<sampleName>" where sampleName is the fasta file name minus the extension
    // If the idline already contains " sampleName=<sampleName>" then it is not appended.
    private suspend fun annotateFasta(inputChannel: Channel<InputChannelData>) =
        withContext (
            Dispatchers.Default
        ) {
            for (entry in inputChannel) {
                myLogger.info("annotateFasta: entry = ${entry.sampleName}")
                val fastaFile = entry.fastaFile
                val sampleName = entry.sampleName
                val outputDir = entry.outputDir
                // Write all lines in the fasta file to the new file.
                // the sample name is the filename minus the extension
                // the new file name is the same as the original but it is
                // written to the outputDir
                // if  the line starts with >, append " sampleName=${sampleName}" to the line

                // Previously if the original file was compressed, we compressed the new file.
                // We no longer do this, as the align-assemblies step requires uncompressed fasta files.
                // This function will only create un-compressed files.
                // The file names will be <sampleName>.fa
                val justName =  "${sampleName}.fa"
                val newFilename = "${outputDir}/${justName}"
                File(newFilename).bufferedWriter().use { writer ->
                    bufferedReader(fastaFile).forEachLine { line ->
                        if (line.startsWith(">")) {
                            // Only append if the line does not already contain "sampleName=${sampleName}"
                            if (!line.contains("sampleName=${sampleName}") )
                                writer.write("${line} sampleName=${sampleName}\n")
                            else
                                writer.write("${line}\n")
                        } else {
                            writer.write("${line}\n")
                        }
                    }
                }
            }
        }

}