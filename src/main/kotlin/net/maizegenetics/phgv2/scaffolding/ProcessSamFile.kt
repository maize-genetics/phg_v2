package net.maizegenetics.phgv2.scaffolding

import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.ValidationStringency
import java.io.BufferedWriter
import java.io.File

class ProcessSamFile : CliktCommand(help = "Process a SAM file to extract scaffold information") {

    val samFile by option(help = "Path to the SAM file")
        .required()

    val outputDir by option(help = "Output Directory")
        .required()

    override fun run() {
        // Implement the logic to process the SAM file
        // This could involve reading the file, extracting relevant data,
        // and writing results to an output file or directory.

        processSAMFile(samFile, outputDir)
    }

    fun processSAMFile(samFile: String, outputDir: String) {
        File(outputDir).mkdirs() // Ensure the output directory exists
        // Logic to process the SAM file goes here
        SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(File(samFile)).use { reader ->
            var numSamRecords = 0

            val outputFirstSNP = "$outputDir/firstSNP.txt"
            val outputEditDist = "$outputDir/editDistance.txt"
            val outputNumHits = "$outputDir/numHits.txt"

            val editDistances = mutableListOf<Int>()
            val firstSNPPos = mutableListOf<Int>()
            var numHitsAbove2 = 0
            var numReadsClipped = 0

            // Create output files
            bufferedWriter(outputFirstSNP).use { firstSNPWriter ->
                bufferedWriter(outputEditDist).use { editDistancesWriter ->
                    bufferedWriter(outputNumHits).use { numHitsWriter ->
                        val iterator = reader.iterator()
                        var currentReadList = mutableListOf<SAMRecord>(iterator.next())
                        while (iterator.hasNext()) {
                            val record = iterator.next()

                            if(currentReadList.first().readName != record.readName) {
                                numSamRecords++
                                val mapStatus = processSamRecords(
                                    currentReadList,
                                    numHitsWriter,
                                    firstSNPPos,
                                    firstSNPWriter,
                                    editDistances,
                                    editDistancesWriter
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
                                currentReadList = mutableListOf(record)
                            }
                            else {
                                currentReadList.add(record)
                            }
                        }
                        numSamRecords++
                        // Process the last read group
                        val mapStatus = processSamRecords(
                            currentReadList,
                            numHitsWriter,
                            firstSNPPos,
                            firstSNPWriter,
                            editDistances,
                            editDistancesWriter
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

            // Output summary statistics
            println("Total reads clipped: $numReadsClipped")
            println("Original NumReads: ${numSamRecords}")
            println("NumReads with hits > 2: ${numHitsAbove2}")
            println("Total reads: ${editDistances.size}")
            println("Total reads with edit distance == 0 : ${editDistances.count { it == 0 }}")
            println("Average edit distance: ${editDistances.average()}")
            println("Average first SNP position: ${firstSNPPos.average()}")
        }
    }


    fun processSamRecords(
        currentReadList: MutableList<SAMRecord>,
        numHitsWriter: BufferedWriter,
        firstSNPPos: MutableList<Int>,
        firstSNPWriter: BufferedWriter,
        editDistances: MutableList<Int>,
        editDistancesWriter: BufferedWriter
    ): MAPSTATUS {
       //figure out optimal hits based on edit distance
        //remove reads that are clipped
        val filteredReads = currentReadList.filter { it.readUnmappedFlag.not() && it.cigarString.isNotEmpty() }
            .filter { it.cigar.firstOrNull()?.operator != CigarOperator.S && it.cigar.lastOrNull()?.operator != CigarOperator.S
                        && it.cigar.firstOrNull()?.operator != CigarOperator.H && it.cigar.lastOrNull()?.operator != CigarOperator.H }

        if( filteredReads.isEmpty()) {
            //means we have all clipped reads
            return MAPSTATUS.CLIPPED
        }

        val bestNM = filteredReads.minOfOrNull { it.getIntegerAttribute("NM") } ?: Int.MAX_VALUE
        val bestReads = filteredReads.filter { it.getIntegerAttribute("NM") == bestNM }


        //write the number of hits
        numHitsWriter.write("${bestReads.size}\n")

        if( bestReads.size > 2) {
            //if there are more than 2 hits, we skip processing
            numHitsWriter.write("${bestReads.size}\n")
            return MAPSTATUS.MAP_OVER_2
        }
        //Add the edit distance as well
        val editDistance = bestReads.first().getIntegerAttribute("NM")
        editDistances.add(editDistance)
        editDistancesWriter.write("$editDistance\n")

        if (editDistance == 0) {
            //if the edit distance is 0, we can skip processing
            firstSNPWriter.write("0\n")
        }
        else {
            for(alignment in bestReads) {
                //find the first SNP position
                var length = 0
                //Loop through the cigar elements.  If we see an '=' operator we add to the length.  If anything else we add the current length to the firstSNPPos list
                for (element in alignment.cigar) {
                    if (element.operator == CigarOperator.EQ) {
                        length += element.length
                    } else if (element.operator == CigarOperator.X || element.operator == CigarOperator.I || element.operator == CigarOperator.D) {
                        //if we see an insertion or deletion, we add the current length to the firstSNPPos list
                        firstSNPPos.add(length)
                        firstSNPWriter.write("$length\n")
                    }
                }
            }
        }
        return MAPSTATUS.MAPPED
    }
}