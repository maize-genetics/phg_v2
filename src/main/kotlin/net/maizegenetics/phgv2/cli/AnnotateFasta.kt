package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
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
                "--assembly-file must not be blank"
            }
        }

    override fun run() {
        createAnnotatedFastas(fastaList, outputDir)
    }

    fun createAnnotatedFastas(fastaList:String, outputDir:String) {
        val fastaFiles = File(fastaList).readLines()
        fastaFiles.forEach { fastaFile ->
            // Write all lines in the fasta file to the new file.
            // the sample name is the filename minus the extension
            // the new file name is the same as the original but it is
            // written to the outputDir
            // if  the line starts with >, append " sampleName=${sampleName}" to the line
            val sampleName = File(fastaFile).nameWithoutExtension
            val newFilename = "${outputDir}/${File(fastaFile).name}"
            File(newFilename).bufferedWriter().use { writer ->
                File(fastaFile).forEachLine { line ->
                    if (line.startsWith(">")) {
                        writer.write("${line} sampleName=${sampleName}\n")
                    } else {
                        writer.write("${line}\n")
                    }
                }
            }
        }
    }
}