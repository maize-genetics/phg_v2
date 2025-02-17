package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.seqIO.NucSeqIO
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import org.apache.logging.log4j.LogManager
import java.io.BufferedWriter
import java.io.File

class RopeBwtChrIndex: CliktCommand( help = "Index a chromosome for RopeBwt") {
    private val myLogger = LogManager.getLogger(RopeBwtChrIndex::class.java)

    val keyfile by option(help = "Tab-delimited file containing 2 columns name Fasta and SampleName.  Fasta column contains full path name for the fasta files.  SampleName contains the sample name for that assembly, e.g. B73 or CML247. ")
        .required()

    val outputDir by option(help = "Output Directory")
        .required()

    val indexFilePrefix by option(help = "Prefix of the ropebwt3 index file.  This will be added to the output directory and ropeBWT3 will make a number of output files.")
        .required()

    val threads by option(help = "Number of threads to use for the index creation.")
        .int()
        .default(3)

    val deleteFmrIndex by option(help = "Delete the fmr index file after conversion to fmd.")
        .flag(default = true)

    val condaEnvPrefix by option (help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
        .default("")

    override fun run() {
        createChrIndex(keyfile, outputDir, indexFilePrefix, threads, deleteFmrIndex, condaEnvPrefix)
    }

    fun createChrIndex(keyfile: String, outputDir: String, indexFilePrefix: String, threads: Int, deleteFmrIndex: Boolean, condaEnvPrefix: String) {
        myLogger.info("Creating Rename Fasta directory")
        val renameFastaDir = "$outputDir/renamedFastas/"
        File("$outputDir/renamedFastas/").mkdirs()

        val allSeqLengths = mutableListOf<Pair<String,Int>>()

        myLogger.info("Parsing Key File")
        val keyFileParsed = parseKeyFile(keyfile)
        for(keyFileRecord in keyFileParsed) {
            myLogger.info("Indexing ${keyFileRecord.first} with sampleName ${keyFileRecord.second}")
            val (renamedFile, outputSeqLengths) = processKeyFileRecord(keyFileRecord.first, keyFileRecord.second, renameFastaDir)
            allSeqLengths.addAll(outputSeqLengths)

            myLogger.info("Indexing ${renamedFile}")
            addSeqToIndex(renamedFile, "$outputDir/$indexFilePrefix", threads, condaEnvPrefix)
        }

        RopeBWTUtils.convertBWTIndex(indexFilePrefix, condaEnvPrefix)

        //Delete the FMR if requested
        if (deleteFmrIndex) {
            RopeBWTUtils.deleteFMRIndex(indexFilePrefix)
        }

        //Build suffix array
        RopeBWTUtils.buildSuffixArray(indexFilePrefix, threads, condaEnvPrefix)

        //Write the contig lengths to a file
        buildChrLengthFile("$outputDir/$indexFilePrefix", allSeqLengths)

    }

    fun parseKeyFile(keyfile: String): List<Pair<String, String>> {
        return File(keyfile).bufferedReader().readLines()
            .filter { !it.startsWith("#") }
            .map { line ->
                val parts = line.split("\t")
                if (parts.size != 2) {
                    throw IllegalArgumentException("Invalid format in line: \"$line\". Each line must contain exactly one tab between filename and samplename columns.")
                }
                Pair(parts[0], parts[1])
            }
    }

    fun processKeyFileRecord(fastaFile: String, sampleName: String, renameFastaDir: String) : Pair<String,List<Pair<String,Int>>> {
        //Rename the contigs in the fasta file by the sample name
        val fastaFileName = File(fastaFile).nameWithoutExtension
        val renameFastaFile = "$renameFastaDir/${fastaFileName}_renamed.fa"
        val contigLengthPairs = mutableListOf<Pair<String,Int>>()

        bufferedWriter(renameFastaFile).use { writer ->
            renameFastaSeqs(fastaFile, sampleName, writer, contigLengthPairs)
        }

        return Pair(renameFastaFile,contigLengthPairs)
    }

    fun renameFastaSeqs(
        fastaFile: String,
        sampleName: String,
        writer: BufferedWriter,
        contigLengthPairs: MutableList<Pair<String, Int>>)
    {
        NucSeqIO(fastaFile).map { nucSeq ->
            val contigName = nucSeq.name
            val contigSequence = nucSeq.sequence
            val contigLength = contigSequence.size()
            //Contig needs to be chr_sampleGamete
            val outputName = "${contigName}_${sampleName}"
            writer.write(">$outputName\n$contigSequence\n")
            contigSequence.seq().chunked(80).forEach { chunk -> writer.write("$chunk\n") }
            contigLengthPairs.add(Pair(outputName, contigLength))
        }
    }

    fun addSeqToIndex(inputFasta: String, indexFilePrefix: String, threads: Int, condaEnvPrefix: String) {
        val isFirst = !File("$indexFilePrefix.fmd").exists()

        if(isFirst) {
           //Use the normal command
            RopeBWTUtils.runBuildStep(inputFasta, indexFilePrefix, threads, condaEnvPrefix)
        }
        else {
            //Use the update command
            RopeBWTUtils.runBuildUpdateStep(inputFasta, indexFilePrefix, threads, condaEnvPrefix)
        }
    }

    fun buildChrLengthFile(indexFilePrefix: String, contigLengthPairs: List<Pair<String,Int>>) {
        bufferedWriter("$indexFilePrefix.fmd.len.gz").use { writer ->
            contigLengthPairs.forEach { (contigName, contigLength) ->
                writer.write("$contigName\t$contigLength\n")
            }
        }
    }
}
