package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.*
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

    val hvcfDir by option(help = "Directory containing hvcf files used to build the HaplotypeGraph.")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--hvcf-dir must not be blank"
            }
        }

    val kmerIndex by option(help = "Kmer index file created by build-kmer-index.")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--kmer-index must not be blank"
            }
        }



    val readInputFiles: ReadInputFile by mutuallyExclusiveOptions<ReadInputFile>(
        option("--key-file", help = "Name of tab-delimited key file.  Columns for samplename and filename are required.  If using paired end fastqs, a filename2 column can be included").convert{ ReadInputFile.KeyFile(it) },
        option("--read-files", help = "Comma separated list of fastq files for a single sample.  Either 1(for single end) or 2(for paired end) files can be input at a time this way.  Any more and an error will be thrown.").convert{ ReadInputFile.ReadFiles(it) }
    ).single().required()


    val outputDir by option("-o", "--output-dir", help = "Name for output ReadMapping file Directory")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--output-dir/-o must not be blank"
            }
        }



    override fun run() {
        myLogger.info("Begin mapping reads to the pangenome kmer index.")
        alignReadsToHaplotypes(hvcfDir, kmerIndex, readInputFiles.getReadFiles(), outputDir)
    }
}