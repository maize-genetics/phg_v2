package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import java.io.File

enum class AlignmentClass {
    PAIRUNIQUE, PAIRRARE, PAIRCOMMON, PAIRSPLIT, PAIROFFASM, PAIRUNALIGN,
    UNIQUE_RARE, UNIQUE_COMMON, UNIQUE_SPLIT, UNIQUE_OFFASM, UNIQUE_UNALIGN,
    RARE_COMMON, RARE_SPLIT, RARE_OFFASM, RARE_UNALIGN,
    COMMON_SPLIT, COMMON_OFFASM, COMMON_UNALIGN,
    SPLIT_OFFASM, SPLIT_UNALIGN,
    OFFASM_UNALIGN

}
class ExtractEdgeReads : CliktCommand( help = "Extract out Edge Case reads from SAMs/BAMs") {
    val bamDir by option(help = "Folder name where The BAM/SAM files are")
        .required()
    val hvcfDir by option(help = "Folder name where the HVCFS are")
        .required()

    override fun run() {

        val hvcfFiles = File(hvcfDir).listFiles { file -> file.name.endsWith(".h.vcf") || file.name.endsWith(".h.vcf.gz") }.map { it.path }

        //load the graph in
        val graph = HaplotypeGraph(hvcfFiles)

        val bamFiles = File(bamDir).listFiles { file -> file.name.endsWith(".bam") || file.name.endsWith(".sam") }.map { it.path }

        //extract the edge reads
        for(bamFile in bamFiles) {
            extractReads(bamFile, graph)
        }
        TODO("Not yet implemented")
    }

    fun extractReads(bamFile: String, graph: HaplotypeGraph) {
        //load in the reads
        val samReader = SamReaderFactory.makeDefault().open(File(bamFile))
        val iterator = samReader.iterator()
        var currentReadId = ""
        var recordsForRead = mutableListOf<SAMRecord>()

        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()

        while(iterator.hasNext()) {
            val currentRecord = iterator.next()
            if(currentRecord.readName != currentReadId) {
                //process the reads
                processReads(recordsForRead, hapIdToRefRangeMap)
                //reset the records
                recordsForRead = mutableListOf()
                currentReadId = currentRecord.readName
            }
            recordsForRead.add(currentRecord)
        }
        //process the last read alignments
        processReads(recordsForRead, hapIdToRefRangeMap)
    }

    fun processReads(recordsForRead: List<SAMRecord>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>) {
        //Pair off the reads by their alignment to haplotype ids
        val recordsGroupedByContig = recordsForRead.groupBy { record -> record.contig }.map { filterAlignmentToPair(it.value) }
        //For the pair only keep track of the best ones based on edit distance
        //We are looking for various edge cases
        classifyAlignments(recordsGroupedByContig)
    }

    fun filterAlignmentToPair(records: List<SAMRecord>): Pair<SAMRecord?,SAMRecord?> {
        //these records are all hitting the same contig.  Need to split them by first in pair and second in pair
        val bestAlignments = records.groupBy { it.firstOfPairFlag }.map { keepBestAlignment(it.value) }
        return Pair(bestAlignments[0], bestAlignments[1])
    }

    fun keepBestAlignment(records : List<SAMRecord>) : SAMRecord? {
        return records.minByOrNull { it.getIntegerAttribute("NM") }
    }

    fun classifyAlignments(records: List<Pair<SAMRecord?,SAMRecord?>>): AlignmentClass {
        return when {
            records.size == 1 && records.first().first == null && records.first().second == null -> AlignmentClass.PAIRUNALIGN
            records.size == 1 && (records.first().first == null || records.first().second == null) -> AlignmentClass.UNIQUE_UNALIGN
            records.size == 1 -> AlignmentClass.PAIRUNIQUE

            else -> AlignmentClass.PAIRUNALIGN
        }

    }
}