package net.maizegenetics.phgv2.scaffolding

import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.TextCigarCodec
import java.io.BufferedWriter
import java.io.File

enum class MAPSTATUS {
    MAPPED,
    UNMAPPED,
    CLIPPED,
    MAP_OVER_2,
}

class ProcessPafFile : CliktCommand(help = "Process a PAF file to extract scaffold information") {

    val pafFile by option(help = "Path to the PAF file")
        .required()

    val outputDir by option(help = "Output Directory")
        .required()

    override fun run() {
        // Implement the logic to process the PAF file
        // This could involve reading the file, extracting relevant data,
        // and writing results to an output file or directory.

        processPAFFile(pafFile, outputDir)
    }

    fun processPAFFile(pafFile: String, outputDir: String) {
        val pafLinesSplit = bufferedReader(pafFile).readLines().map { it.split("\t") }
        File(outputDir).mkdirs()

        val outputFirstSNP = "$outputDir/firstSNP.txt"
        val outputEditDist = "$outputDir/editDistance.txt"
        val outputNumHits = "$outputDir/numHits.txt"

        val editDistances = mutableListOf<Int>()
        val firstSNPPos = mutableListOf<Int>()
        var numHitsAbove2 = 0
        var numReadsClipped = 0

        bufferedWriter(outputFirstSNP).use { firstSNPWriter ->
            bufferedWriter(outputEditDist).use { editDistancesWriter ->
                bufferedWriter(outputNumHits).use { numHitsWriter ->
                    for (line in pafLinesSplit) {
                        val mapStatus = processPafRecord(
                            line,
                            numHitsWriter,
                            firstSNPPos,
                            firstSNPWriter,
                            editDistances,
                            editDistancesWriter,
                        )
                        when (mapStatus) {
                            MAPSTATUS.MAPPED -> {
                                // Do nothing, we already processed the record
                            }
                            MAPSTATUS.UNMAPPED -> {
                                // Handle unmapped case if needed
                            }
                            MAPSTATUS.CLIPPED -> {
                                numReadsClipped++
                            }
                            MAPSTATUS.MAP_OVER_2 -> {
                                numHitsAbove2++
                            }
                        }
                    }
                }
            }
        }

        // Output summary statistics
        println("Total reads clipped: $numReadsClipped")
        println("Original NumReads: ${pafLinesSplit.size}")
        println("NumReads with hits > 2: ${numHitsAbove2}")
        println("Total reads: ${editDistances.size}")
        println("Total reads with edit distance == 0 : ${editDistances.count { it == 0 }}")
        println("Average edit distance: ${editDistances.average()}")
        println("Average first SNP position: ${firstSNPPos.average()}")


    }

    private fun processPafRecord(
        line: List<String>,
        numHitsWriter: BufferedWriter,
        firstSNPPos: MutableList<Int>,
        firstSNPWriter: BufferedWriter,
        editDistances: MutableList<Int>,
        editDistancesWriter: BufferedWriter
    ): MAPSTATUS {
        //index 15 is cigar string
        val cigarString = line[15].split(":").last()

        val numHits = line[14].split(":").last().toInt()
        numHitsWriter.write("$numHits\n")
        if (numHits <= 2) {
            //parse cigar string
            val cigarElements = TextCigarCodec.decode(cigarString).cigarElements

            var length = 0
            if (cigarElements.first().operator == CigarOperator.S || cigarElements.first().operator == CigarOperator.H
                || cigarElements.last().operator == CigarOperator.S || cigarElements.last().operator == CigarOperator.H
            ) {
                return MAPSTATUS.CLIPPED
            }
            var runningEditDist = 0
            //Loop through the cigar elements.  If we see an '=' operator we add to the length.  If anything else we add the current length to the firstSNPPos list
            for (element in cigarElements) {
                if (element.operator == CigarOperator.EQ) {
                    length += element.length
                } else if (element.operator == CigarOperator.X || element.operator == CigarOperator.I || element.operator == CigarOperator.D) {
                    //if we see an insertion or deletion, we add the current length to the firstSNPPos list
                    firstSNPPos.add(length)
                    firstSNPWriter.write("$length\n")
                    runningEditDist += element.length
                }
            }
            //add the running edit distance to the editDistances list
            editDistances.add(runningEditDist)
            editDistancesWriter.write("$runningEditDist\n")
            return MAPSTATUS.MAPPED
        } else {
            return MAPSTATUS.MAP_OVER_2
        }

    }
}