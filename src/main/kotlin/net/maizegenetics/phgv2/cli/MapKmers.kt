package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.float
import net.maizegenetics.phgv2.utils.KeyFileData
import net.maizegenetics.phgv2.utils.alignReadsToHaplotypes
import java.io.File
import java.util.logging.Logger

/**
 * Sealed classes to handle either a keyFile or a list of readFiles.

 */
sealed class ReadInputFile {
    abstract fun getReadFiles() : List<KeyFileData>
    data class KeyFile(val keyFile: String): ReadInputFile() {
        @Override
        override fun getReadFiles(): List<KeyFileData> {
            check(File(keyFile).exists()) { "Key file $keyFile does not exist." }
            val linesWithHeader = File(keyFile).bufferedReader().readLines()
            val header = linesWithHeader.first().split("\t")
            //convert the header into a map of column name to column index
            val headerMap = header.mapIndexed { index, s -> s to index }.toMap()
            check(headerMap.containsKey("sampleName")) { "Key file $keyFile must have a column named sampleName." }
            check(headerMap.containsKey("filename")) { "Key file $keyFile must have a column named filename." }
            return linesWithHeader.drop(1).map{lines -> lines.split("\t")}.map { linesSplit ->
                KeyFileData(linesSplit[headerMap["sampleName"]!!], linesSplit[headerMap["filename"]!!], if(headerMap.containsKey("filename2") && linesSplit.indices.contains(headerMap["filename2"]!!)) linesSplit[headerMap["filename2"]!!] else "")
            }
        }
    }
    data class ReadFiles(val readFiles: String): ReadInputFile() {
        @Override
        override fun getReadFiles(): List<KeyFileData> {
            check(readFiles.isNotEmpty()) { "--read-files must have at least one file." }
            val fileNames = readFiles.split(",")
            check(fileNames.size <= 2) { "--read-files must have 1 or 2 files separated by commas.  You provided: ${fileNames.size}" }
            return listOf(KeyFileData("noSample",fileNames.first(), if(fileNames.size==1) "" else fileNames.last()))
        }

    }
}

class MapKmers : CliktCommand(help="Map Kmers to the pangenome reference") {

    private val myLogger = Logger.getLogger("net.maizegenetics.phgv2.cli.MapKmers")

    //./phg map-kmers \
    //    --kmer-index kmer_index.map \
    //    --reads my_reads.fastq \ // possibly thousands of samples being inputted
    //    --paired
    //    --output read_count_out.map \ // could we pipe this into impute method? // thousands of outputs
    //    // consider batch interface here ^^

    val hvcfDir by option(help = "Directory containing hvcf which build up the Haplotype Graph.")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--hvcf-dir must not be blank"
            }
        }

    val kmerIndex by option(help = "Kmer index file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--kmer-index must not be blank"
            }
        }



    val readInputFiles: ReadInputFile by mutuallyExclusiveOptions<ReadInputFile>(
        option("--key-file").convert{ ReadInputFile.KeyFile(it) },
        option("--read-files").convert{ ReadInputFile.ReadFiles(it) }
    ).single().required()



    val paired by option(help = "Flag to indicate if reads are paired end").flag(default=false)


    val outputDir by option("-o", "--output-dir", help = "Name for output ReadMapping file Directory")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--output-dir/-o must not be blank"
            }
        }



    override fun run() {
        myLogger.info("Begin mapping reads to the pangenome kmer index.")
        alignReadsToHaplotypes(hvcfDir, kmerIndex, readInputFiles.getReadFiles(), paired, outputDir)
    }
}