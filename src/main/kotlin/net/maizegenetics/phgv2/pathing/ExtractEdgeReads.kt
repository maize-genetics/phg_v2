package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter

enum class AlignmentClass {
    PAIRUNIQUE, PAIRRARE, PAIRCOMMON, PAIRREADSPLIT, PAIRREADSPLITCONSEC, PAIRALIGNSPLIT, PAIRALIGNSPLITCONSEC, PAIROFFASM, UNALIGN,
    SINGLEUNIQUE, SINGLERARE, SINGLECOMMON, SINGLEALIGNSPLIT, SINGLEALIGNSPLITCONSEC,SINGLEOFFASM

}

/**
 * This class is to find edge case reads for use as integration tests for the read mapping pipeline
 *
 * IT IS CURRENTLY A WIP AND IS UNDER ACTIVE DEVELOPMENT
 */
class ExtractEdgeReads : CliktCommand( help = "Extract out Edge Case reads from SAMs/BAMs") {
    val bamDir by option(help = "Folder name where The BAM/SAM files are")
        .required()
    val hvcfDir by option(help = "Folder name where the HVCFS are")
        .required()
    val sampleName by option(help = "Sample name to use for the reads")
        .required()

    override fun run() {

        val hvcfFiles = File(hvcfDir).listFiles { file -> file.name.endsWith(".h.vcf") || file.name.endsWith(".h.vcf.gz") }.map { it.path }

        //load the graph in
        val graph = HaplotypeGraph(hvcfFiles)

        val bamFiles = File(bamDir).listFiles { file -> file.name.endsWith(".bam") || file.name.endsWith(".sam") }.map { it.path }

        //extract the edge reads
        for(bamFile in bamFiles) {
            extractReads(sampleName, bamFile, graph)
        }
        TODO("Not yet implemented")
    }

    fun extractReads(sampleName: String, bamFile: String, graph: HaplotypeGraph) {
        //load in the reads
        val samReader = SamReaderFactory.makeDefault().open(File(bamFile))
        val iterator = samReader.iterator()
        var currentReadId = ""
        var recordsForRead = mutableListOf<SAMRecord>()

        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()
        val hapIdsToSampleGametes = graph.hapIdsToSampleGametes()

        val numSampleGametes = graph.numSampleGametes()

        while(iterator.hasNext()) {
            val currentRecord = iterator.next()
            if(currentRecord.readName != currentReadId) {
                //process the reads
                processReads(sampleName, numSampleGametes, recordsForRead, hapIdToRefRangeMap, hapIdsToSampleGametes)
                //reset the records
                recordsForRead = mutableListOf()
                currentReadId = currentRecord.readName
            }
            recordsForRead.add(currentRecord)
        }
        //process the last read alignments
        processReads(sampleName, numSampleGametes, recordsForRead, hapIdToRefRangeMap, hapIdsToSampleGametes)
    }

    fun processReads(sampleName:String, numSampleGametes: Int, recordsForRead: List<SAMRecord>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, hapIdToSampleGamete: Map<String,List<SampleGamete>>) {
        //Pair off the reads by their alignment to haplotype ids
        val recordsGroupedByContig = recordsForRead.groupBy { record -> hapIdToSampleGamete[record.contig]!! }.map { filterAlignmentToPair(it.value) }
        //For the pair only keep track of the best ones based on edit distance
        //We are looking for various edge cases
        classifyAlignments(sampleName, numSampleGametes,recordsGroupedByContig, hapIdToRefRangeMap)
    }

    fun filterAlignmentToPair(records: List<SAMRecord>): Pair<SAMRecord?,SAMRecord?> {
        //these records are all hitting the same contig.  Need to split them by first in pair and second in pair
        val bestAlignments = records.groupBy { it.firstOfPairFlag }.map { keepBestAlignment(it.value) }
        return Pair(bestAlignments[0], bestAlignments[1])
    }

    fun keepBestAlignment(records : List<SAMRecord>) : SAMRecord? {
        return records.minByOrNull { it.getIntegerAttribute("NM") }
    }

    fun classifyAlignments(sampleName: String, numSampleGametes: Int ,records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>): AlignmentClass {
        return when {
            records.size == 1 -> classifyUniqueAlignments(records)
            isReadSplit(records) -> classifyReadSplitAlignments(records)
            isAlignSplit(records, hapIdToRefRangeMap) -> classifyAlignSplitAlignments(records)
            records.size in 2 .. numSampleGametes/2 -> classifyRareAlignments(sampleName, numSampleGametes ,records, hapIdToRefRangeMap )

            else -> AlignmentClass.UNALIGN
        }

    }

    fun isReadSplit(records: List<Pair<SAMRecord?,SAMRecord?>>): Boolean {
        //If each pair of records hit different haplotypes then it is a read split
        var isReadSplit = false
        for(recordPair in records) {
            //Skip if one of the records is null
            if(recordPair.first == null || recordPair.second == null) {
                continue
            }
            if(recordPair.first?.contig != recordPair.second?.contig) {
                isReadSplit = true
                break
            }
        }
        return isReadSplit
    }

    fun isAlignSplit(records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>): Boolean {
        //Check to see if all the records have a common refRange  if not its an align split
        var refRangeSet = setOf<ReferenceRange>()
        for(record in records) {
            val readOneHapId = record.first?.contig
            val readTwoHapId = record.second?.contig

            val currentRangeSet = mutableSetOf<ReferenceRange>()
            //Need to check nulls here
            if(readOneHapId != null) currentRangeSet.addAll(hapIdToRefRangeMap[readOneHapId]!!)
            if(readTwoHapId != null) currentRangeSet.addAll(hapIdToRefRangeMap[readTwoHapId]!!)

            if(currentRangeSet.isNotEmpty()) {
                refRangeSet = if(refRangeSet.isEmpty()) {
                    currentRangeSet
                } else {
                    refRangeSet.intersect(currentRangeSet)
                }
            }
            else {
                continue
            }

            if(refRangeSet.isEmpty()) {
                return true
            }
        }
        return false
    }


    fun classifyUniqueAlignments(records: List<Pair<SAMRecord?,SAMRecord?>>): AlignmentClass {
        return when {
            records.first().first == null && records.first().second == null -> AlignmentClass.UNALIGN
            (records.first().first == null || records.first().second == null) -> AlignmentClass.SINGLEUNIQUE
            else -> AlignmentClass.PAIRUNIQUE
        }
    }

    fun classifyRareAlignments(sampleName: String, numSampleGametes: Int, records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>): AlignmentClass {
        return when {
            //PAIRRARE
            //
            else -> AlignmentClass.UNALIGN
        }
    }
    fun classifyReadSplitAlignments(records: List<Pair<SAMRecord?,SAMRecord?>>): AlignmentClass {
        TODO()
    }
    fun classifyAlignSplitAlignments(records: List<Pair<SAMRecord?,SAMRecord?>>): AlignmentClass {
        TODO()
    }

    /**
     * Function that exports the filtered down fastq files and a table of [readID] -> [HapID hits] to outputFileName
     */
    fun outputReadsAndHapIdSets(
        outputFileName: String,
        // key = readId, values = alignments paired off by the correct strand that we want to export
        alignments: Map<String, List<Pair<SAMRecord, SAMRecord>>>
    ) {
        // Write the table to the file as tab-delimited text
        BufferedWriter(FileWriter(outputFileName)).use { writer ->
            writer.write("readID\tHapID hits")
            for ((readID, samlist) in alignments) {
                // add all HapIDs to hits
                val hits : MutableSet<String> = mutableSetOf<String>()
                for (pair in samlist) {
                    // SAMRecord.getContig() gets HapID
                    hits.add(pair.first.getContig())
                    hits.add(pair.second.getContig())
                }
                writer.write("$readID\t${hits.joinToString(separator = "\t")}\n")
            }
        }
        TODO("export filtered down fastq files")
    }
}