package net.maizegenetics.phgv2.utils

import biokotlin.util.bufferedReader
import java.io.BufferedWriter
import java.io.FileWriter

data class MAFRecord(val chr: String, val start: Int, val length: Int, val strand: String, val sequence: String)

fun testMergingMAF(inputFile: String, outputFile: String) {

    val reader = bufferedReader(inputFile)

    var currentLine = reader.readLine()
    val mafRecordList = mutableListOf<Pair<MAFRecord, MAFRecord>>()
    while (currentLine != null) {
        if (currentLine == "" || currentLine.startsWith("#")) {
            currentLine = reader.readLine()
            continue
        }

        if (currentLine.startsWith("a")) {
            val score = currentLine
            val firstRecord = reader.readLine()
            val secondRecord = reader.readLine()

            val firstRecordFields = firstRecord.split("\t")
            val secondRecordFields = secondRecord.split("\t")

            val firstRecordContig = firstRecordFields[1].trim()
            val firstRecordStart = firstRecordFields[2].trim().toInt()
            val firstRecordLength = firstRecordFields[3].trim().toInt()
            val seq = firstRecordFields[6].trim()

            val secondRecordContig = secondRecordFields[1].trim()
            val secondRecordStart = secondRecordFields[2].trim().toInt()
            val secondRecordLength = secondRecordFields[3].trim().toInt()
            val secondSeq = secondRecordFields[6].trim()

            mafRecordList.add(
                Pair(
                    MAFRecord(firstRecordContig, firstRecordStart, firstRecordLength, firstRecordFields[4].trim(), seq),
                    MAFRecord(
                        secondRecordContig,
                        secondRecordStart,
                        secondRecordLength,
                        secondRecordFields[4].trim(),
                        secondSeq
                    )
                )
            )

        }
        currentLine = reader.readLine()
    }

    // group by chromosome
    val groupedByChromosome = mafRecordList.groupBy { it.first.chr }
    // sort each group by start position
    val sortedByStart = groupedByChromosome.mapValues { (_, value) -> value.sortedBy { it.first.start } }

    // BufferedWriter(FileWriter("/Users/zrm22/Desktop/debug_mafToGVCF/testoutput_nonSplit.txt")).use { output ->
    BufferedWriter(FileWriter(outputFile)).use { output ->
        // for each group build the sequence
        for (chr in sortedByStart.keys) {
            val records = sortedByStart[chr]!!

            // loop through the records and concatenate the ref and asm seqs
            val refSeqBuilder = StringBuilder()
            val asmSeqBuilder = StringBuilder()

            var totalRefLength = 0
            var totalAsmLength = 0
            val startRef = records.first().first.start
            val startAsm = records.first().second.start
            for (record in records) {
                refSeqBuilder.append(record.first.sequence)
                asmSeqBuilder.append(record.second.sequence)
                totalRefLength += record.first.length
                totalAsmLength += record.second.length
            }
            output.write(">${chr}_ref $startRef ${totalRefLength}\n")
            output.write(refSeqBuilder.toString().chunked(80).joinToString("\n"))
            output.write("\n")
            output.write(">${chr}_asm $startAsm ${totalAsmLength}\n")
            output.write(asmSeqBuilder.toString().chunked(80).joinToString("\n"))
            output.write("\n")

        }
    }

}