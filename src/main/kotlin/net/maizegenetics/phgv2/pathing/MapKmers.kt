package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.cli.logCommand
import org.apache.logging.log4j.LogManager
import java.io.File

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
            return linesWithHeader.drop(1).map{
                keyFileDataFromLine(it, headerMap)
            }
        }
    }
    data class ReadFiles(val readFiles: String): ReadInputFile() {
        @Override
        override fun getReadFiles(): List<KeyFileData> {
            check(readFiles.isNotEmpty()) { "--read-files must have at least one file." }
            val fileNames = readFiles.split(",")
            check(fileNames.size <= 2) { "--read-files must have 1 or 2 files separated by commas.  You provided: ${fileNames.size}" }
            check(fileNames[0].endsWith(".fq") || fileNames[0].endsWith(".fq.gz") || fileNames[0].endsWith(".fastq") || fileNames[0].endsWith(".fastq.gz")) {"file ${fileNames[0]} does not end in one of .fq, .fq.gz, .fastq, or .fastq.gz"}
            val fileBase = File(fileNames[0]).name.removeSuffix(".gz").removeSuffix(".fq").removeSuffix(".fastq")

            return listOf(KeyFileData(fileBase,fileNames.first(), if(fileNames.size==1) "" else fileNames.last()))
        }

    }

    fun keyFileDataFromLine(keyfileLine: String, headerMap: Map<String, Int>):KeyFileData {
        //possible headers are sampleName, filename, and filename2
        //check that filename has a fastq type extension
        //if it does assume filename2 does as well

        val splitLine = keyfileLine.split("\t")
        val file1 = splitLine[headerMap["filename"]!!]
        check(file1.endsWith(".fq") || file1.endsWith(".fq.gz") || file1.endsWith(".fastq") || file1.endsWith(".fastq.gz")) {"filename for ${splitLine[headerMap["sampleName"]!!]} does not end in one of .fq, .fq.gz, .fastq, or .fastq.gz"}
        val file2 = if (headerMap.contains("filename2") && splitLine.size >= headerMap["filename2"]!! + 1) splitLine[headerMap["filename2"]!!] else ""
        return KeyFileData(splitLine[headerMap["sampleName"]!!], file1, file2)
    }
}

class MapKmers : CliktCommand(help="Map Kmers to the pangenome reference") {

    private val myLogger = LogManager.getLogger(MapKmers::class.java)

    //./phg map-kmers \
    //    --kmer-index kmer_index.map \
    //    --reads my_reads.fastq \ // possibly thousands of samples being inputted
    //    --paired
    //    --output read_count_out.map \ // could we pipe this into impute method? // thousands of outputs
    //    // consider batch interface here ^^

    val hvcfDir by option(help = "Directory containing hvcf files used to build the HaplotypeGraph. (Required)")
        .required()

    val kmerIndex by option(help = "Kmer index file created by build-kmer-index. Default is <hvcfDir>/kmerIndex.txt.")
        .default("")

    val readInputFiles: ReadInputFile by mutuallyExclusiveOptions<ReadInputFile>(
        option("--key-file", help = "Name of tab-delimited key file.  Columns for samplename and filename are required.  If using paired end fastqs, a filename2 column can be included. A value must be entered for either --key-file or --read-files.").convert{ ReadInputFile.KeyFile(it) },
        option("--read-files", help = "Comma separated list of fastq files for a single sample.  Either 1(for single end) or 2(for paired end) files can be input at a time this way.  Any more and an error will be thrown.").convert{ ReadInputFile.ReadFiles(it) }
    ).single().required()


    val outputDir by option("-o", "--output-dir", help = "Name for output ReadMapping file Directory (Required)")
        .required()

    val threads by option(help = "Number of threads to use.")
        .int()
        .default(5)

    val minProportionOfMaxCount by option(help = "Minimum proportion of the maximum count for a read to be considered a match.")
        .double()
        .default(1.0)


    val minProportionSameReferenceRange by option(help = "Minimum proportion of the read that must align to the same reference range.")
        .double()
        .default(0.0) //Setting this to default of 0.0 as increasing it will filter out a lot of reads

    val diagnosticMode by option(help = "Run in diagnostic mode.  This will output additional information to a diagnostic log file for each read processed.  " +
            "The file name is automatically generated by appending a '_diagnostic.txt' to the end of the input fastq file name.  " +
            "This is intended to be used only during debugging as it will produce a lot of output.")
        .flag(default = false)

    override fun run() {

        logCommand(this)

        myLogger.info("Begin mapping reads to the pangenome kmer index.")
        // set the kmerIndex file name
        val kmerIndexFilename = kmerIndex.ifBlank { "${hvcfDir}/kmerIndex.txt" }

        // create a HaplotypeGraph from the list of hvcf files
        val graph = HaplotypeGraph(hvcfDir)
        AlignmentUtils.alignReadsToHaplotypes(graph, kmerIndexFilename, readInputFiles.getReadFiles(), outputDir, threads, minProportionOfMaxCount, minProportionSameReferenceRange, diagnosticMode)
    }
}