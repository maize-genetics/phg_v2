package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.float
import net.maizegenetics.phgv2.utils.alignReadsToHaplotypes
import java.util.logging.Logger

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

    val readFiles by option(help = "Fastq file(s) to be aligned.  Can be a comma separated list of files.")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--read-files must not be blank"
            }
        }

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
        alignReadsToHaplotypes(hvcfDir, kmerIndex, readFiles, paired, outputDir)
    }
}