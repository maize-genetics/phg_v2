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
    PAIRUNIQUE, PAIRRARE, PAIRCOMMON, PAIRSPLIT, PAIROFFASM, PAIRUNALIGN,
    UNIQUE_RARE, UNIQUE_COMMON, UNIQUE_SPLIT, UNIQUE_OFFASM, UNIQUE_UNALIGN,
    RARE_COMMON, RARE_SPLIT, RARE_OFFASM, RARE_UNALIGN,
    COMMON_SPLIT, COMMON_OFFASM, COMMON_UNALIGN,
    SPLIT_OFFASM, SPLIT_UNALIGN,
    OFFASM_UNALIGN

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

        while(iterator.hasNext()) {
            val currentRecord = iterator.next()
            if(currentRecord.readName != currentReadId) {
                //process the reads
                processReads(sampleName, recordsForRead, hapIdToRefRangeMap, hapIdsToSampleGametes)
                //reset the records
                recordsForRead = mutableListOf()
                currentReadId = currentRecord.readName
            }
            recordsForRead.add(currentRecord)
        }
        //process the last read alignments
        processReads(sampleName, recordsForRead, hapIdToRefRangeMap, hapIdsToSampleGametes)
    }

    fun processReads(sampleName:String, recordsForRead: List<SAMRecord>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, hapIdToSampleGamete: Map<String,List<SampleGamete>>) {
        //Pair off the reads by their alignment to haplotype ids
        val recordsGroupedByContig = recordsForRead.groupBy { record -> hapIdToSampleGamete[record.contig]!! }.map { filterAlignmentToPair(it.value) }
        //For the pair only keep track of the best ones based on edit distance
        //We are looking for various edge cases
        classifyAlignments(sampleName, recordsGroupedByContig, hapIdToRefRangeMap)
    }

    fun filterAlignmentToPair(records: List<SAMRecord>): Pair<SAMRecord?,SAMRecord?> {
        //these records are all hitting the same contig.  Need to split them by first in pair and second in pair
        val bestAlignments = records.groupBy { it.firstOfPairFlag }.map { keepBestAlignment(it.value) }
        return Pair(bestAlignments[0], bestAlignments[1])
    }

    fun keepBestAlignment(records : List<SAMRecord>) : SAMRecord? {
        return records.minByOrNull { it.getIntegerAttribute("NM") }
    }

    fun classifyAlignments(sampleName: String, records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>): AlignmentClass {
        return when {
            records.size == 1 && records.first().first == null && records.first().second == null -> AlignmentClass.PAIRUNALIGN
            records.size == 1 && (records.first().first == null || records.first().second == null) -> AlignmentClass.UNIQUE_UNALIGN
            records.size == 1 -> AlignmentClass.PAIRUNIQUE
            records.size in 2 .. 5 -> classifyRareAlignments(sampleName, records, hapIdToRefRangeMap )

            else -> AlignmentClass.PAIRUNALIGN
        }

    }

    fun classifyRareAlignments(sampleName: String, records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>): AlignmentClass {
        TODO("Finish implementing this.")
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