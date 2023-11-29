package net.maizegenetics.phgv2.agc

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.int
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.launch
import kotlinx.coroutines.runBlocking
import kotlinx.coroutines.withContext
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * This class takes a list of fasta files and updates the idline with the sample name.
 * The sample name is derived from the fasta file name, minus the extension.
 * Fasta files should be named with just the sample name, extension can be either .fa or .fasta.
 *
 * This is consistent naming with the AGG compress code..  AGC loads the fasta name minus extension as
 * the sample name.
 *
 * This just a helper class.  Users may already have updated their fasta files.
 */
class AnnotateFasta : CliktCommand() {
    private val myLogger = LogManager.getLogger(AnnotateFasta::class.java)

    val fastaList by option(help = "File containing full path name for the fasta files, one per line, that will be updated. ")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--fasta-list must not be blank"
            }
        }


    // This is the folder where the output files are written.  They will have the same
    // name as the original fasta files.
    val outputDir by option(help = "Folder where updated fastas will be written.")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--output-dir must not be blank"
            }
        }

    val threads by option(help = "number of threads for running concurrently")
        .int()
        .default(1)


    data class InputChannelData(val fastaFile:String, val sampleName:String, val outputDir:String)

    override fun run() {
        // create list of assemblies to align from the assemblies file")
        myLogger.info("creating assembliesList, calling createParallelANnotatedFastas")
        val assembliesList = File(fastaList).readLines().filter { it.isNotBlank() }
        createParallelAnnotatedFastas(assembliesList, outputDir)
    }

//    fun createAnnotatedFastas(fastaList:String, outputDir:String) {
//        val fastaFiles = File(fastaList).readLines()
//        fastaFiles.forEach { fastaFile ->
//            // Write all lines in the fasta file to the new file.
//            // the sample name is the filename minus the extension
//            // the new file name is the same as the original but it is
//            // written to the outputDir
//            // if  the line starts with >, append " sampleName=${sampleName}" to the line
//            val sampleName = File(fastaFile).nameWithoutExtension
//            val newFilename = "${outputDir}/${File(fastaFile).name}"
//            File(newFilename).bufferedWriter().use { writer ->
//                File(fastaFile).forEachLine { line ->
//                    if (line.startsWith(">")) {
//                        // Only append if the line does not already contain "sampleName=${sampleName}"
//                        if (!line.contains("sampleName=${sampleName}") )
//                             writer.write("${line} sampleName=${sampleName}\n")
//                        else
//                            writer.write("${line}\n")
//                    } else {
//                        writer.write("${line}\n")
//                    }
//                }
//            }
//        }
//    }

    private fun createParallelAnnotatedFastas(assemblies:List<String>, outputDir:String) {
        println("createParallelAnnotatedFastas: assemblies: ${assemblies.joinToString("\n")}")
        runBlocking {
            // Setup
            val inputChannel = Channel<InputChannelData>()
            launch {
                myLogger.info("Adding entries to the inputChannel:")
                println("adding entries to the inputChannel , size of assemblies = ${assemblies.size} ")
                assemblies.forEach { asmFile ->
                    val sampleName = File(asmFile).nameWithoutExtension
                    println("adding ${sampleName} to the inputChannel")
                    inputChannel.send(InputChannelData(asmFile, sampleName, outputDir))
                }
                println("Done adding entries to the inputChannel")
                myLogger.info("Done adding data to the inputChannel")
                inputChannel.close()

                val workerThreads = (1..threads).map {
                    launch { annotateFasta(inputChannel) }
                }
            }
        }
    }

    private suspend fun annotateFasta(inputChannel: Channel<InputChannelData>) =
        withContext (
            Dispatchers.Default
        ) {
            for (entry in inputChannel) {
                println("annotateFasta: entry = ${entry.sampleName}")
                val fastaFile = entry.fastaFile
                val sampleName = entry.sampleName
                val outputDir = entry.outputDir
                // Write all lines in the fasta file to the new file.
                // the sample name is the filename minus the extension
                // the new file name is the same as the original but it is
                // written to the outputDir
                // if  the line starts with >, append " sampleName=${sampleName}" to the line
                val newFilename = "${outputDir}/${File(fastaFile).name}"
                File(newFilename).bufferedWriter().use { writer ->
                    File(fastaFile).forEachLine { line ->
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