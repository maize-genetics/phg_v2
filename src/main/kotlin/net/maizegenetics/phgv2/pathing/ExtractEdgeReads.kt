package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import net.maizegenetics.phgv2.api.HaplotypeGraph
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
        while(iterator.hasNext()) {
            val currentRecord = iterator.next()
            if(currentRecord.readName != currentReadId) {
                //process the reads
                processReads(recordsForRead, graph)
                //reset the records
                recordsForRead = mutableListOf()
                currentReadId = currentRecord.readName
            }
            recordsForRead.add(currentRecord)
        }
        //process the last read alignments
        processReads(recordsForRead, graph)
    }

    fun processReads(recordsForRead: List<SAMRecord>, graph: HaplotypeGraph) {
        //Pair off the reads by their alignment to haplotype ids
        //For the pair only keep track of the best ones based on edit distance
        //We are looking for various edge cases
    }

    fun classifyAlignments(records: List<Pair<SAMRecord,SAMRecord>>): AlignmentClass {
        return when {
            records.size == 1 -> AlignmentClass.PAIRUNIQUE

            else -> AlignmentClass.PAIRUNALIGN
        }

    }
}