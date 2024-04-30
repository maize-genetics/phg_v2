package net.maizegenetics.phgv2.simulation

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import net.maizegenetics.phgv2.utils.getBufferedWriter
import net.maizegenetics.phgv2.utils.retrieveAgcContigs
import org.apache.logging.log4j.LogManager
import java.io.File
import kotlin.random.Random

class GenerateReads: CliktCommand(help="Generate simulated reads") {
    private val myLogger = LogManager.getLogger(GenerateReads::class.java)

    val sampleNames by option(help = "sample names").default("")

    override fun run() {
        for (sampleName in sampleNames.split(",")) {
            singleReadsFromSample(sampleName, "/tiledb_maize")
        }
    }

    fun singleReadsFromSample(sampleName: String, dbPath: String) {
        //read header will be @<sampleName>:<chr>:<start>:<end>
        //then sequence, +, string of "E" same length as sequence
        val readLength = 100
        val outputFile = File("/workdir/output/$sampleName-single-reads.fq")

        myLogger.info("Generating single reads for $sampleName")
        getBufferedWriter(outputFile).use { fastqWriter ->
            for (chr in (1..10)) {
                val sequence = retrieveAgcContigs(dbPath, listOf("chr${chr.toString()}@$sampleName"))
                val nucseq = sequence.entries.first().value
                val chrlen = nucseq.size()
                val numberOfReads = chrlen / 1000
                val startLimit = chrlen - 100
                repeat(numberOfReads) {
                    val start = Random.nextInt(startLimit)
                    val end = start + readLength - 1
                    val readSeq = nucseq[start..end].seq()
                    fastqWriter.write(readSeq + "\n")
                    fastqWriter.write("+\n")
                    fastqWriter.write("${"E".repeat(readLength)}\n")
                }
            }
        }

    }


}

