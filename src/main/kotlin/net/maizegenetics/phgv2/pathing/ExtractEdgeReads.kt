package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import org.apache.logging.log4j.LogManager
import java.io.File

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

    val outputFileDir by option(help = "Output file directory")
        .required()

    val maxClassNum by option(help = "Maximum number of classes to use")
        .int()
        .default(10)


    private val myLogger = LogManager.getLogger(ExtractEdgeReads::class.java)


    override fun run() {

        val hvcfFiles = File(hvcfDir).listFiles { file -> file.name.endsWith(".h.vcf") || file.name.endsWith(".h.vcf.gz") }.map { it.path }

        //load the graph in
        val graph = HaplotypeGraph(hvcfFiles)

        val bamFiles = File(bamDir).listFiles { file -> file.name.endsWith(".bam") || file.name.endsWith(".sam") }.map { it.path }

        //extract the edge reads
        for(bamFile in bamFiles) {
            myLogger.info("Extracting reads from $bamFile")
            val (countMap, recordsToExport) = extractReads(sampleName, bamFile, graph, maxClassNum)
            //Print out the countMap
            printCountMap(countMap)
            val outputFileName = "${outputFileDir}/${File(bamFile).nameWithoutExtension}"
            myLogger.info("Writing output to $outputFileName")
            outputReadsAndHapIdSets("${outputFileName}_table.txt", "${outputFileName}.fastq", recordsToExport)
        }
    }

    fun extractReads(sampleName: String, bamFile: String, graph: HaplotypeGraph, maxClassNum : Int = 10) : Pair<Map<AlignmentClass,Int>, Map<String, List<Pair<SAMRecord?,SAMRecord?>>>> {
        //load in the reads
        val samReader = SamReaderFactory.makeDefault().open(File(bamFile))
        val iterator = samReader.iterator()
        var currentReadId = ""
        var recordsForRead = mutableListOf<SAMRecord>()

        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()
        val hapIdsToSampleGametes = graph.hapIdsToSampleGametes()

        val numSampleGametes = graph.numSampleGametes()
        val refRangeToIndexMap = graph.refRangeStrToIndexMap()

        val countMap = AlignmentClass.entries.associateWith { 0 }.toMutableMap()
        val recordsToExport = mutableMapOf<String, List<Pair<SAMRecord?,SAMRecord?>>>()

        var readCounter = 0
        while(iterator.hasNext()) {
            val currentRecord = iterator.next()

            if(currentReadId == "") {
                currentReadId = currentRecord.readName
            }

            if(currentRecord.readName != currentReadId) {
                //process the reads
                val (classification, records) = processReads(sampleName, numSampleGametes, recordsForRead, hapIdToRefRangeMap, hapIdsToSampleGametes,refRangeToIndexMap)

                countMap[classification] = countMap[classification]!! + 1
                if(countMap[classification]!! < maxClassNum) {
                    recordsToExport[currentReadId] = records
                }

//                if(classification == AlignmentClass.PAIROFFASM && readCounter < 50000) {
//                    println(recordsForRead.first().readName)
//                }


                readCounter++
                if(readCounter % 10000 == 0) {
                    myLogger.info("Processed $readCounter reads")
                }

                //reset the records
                recordsForRead = mutableListOf()
                currentReadId = currentRecord.readName
            }
            recordsForRead.add(currentRecord)
        }
        //process the last read alignments
        val (classification, records) =  processReads(sampleName, numSampleGametes, recordsForRead, hapIdToRefRangeMap, hapIdsToSampleGametes, refRangeToIndexMap)
        countMap[classification] = countMap[classification]!! + 1
        if(countMap[classification]!! < maxClassNum) {
            recordsToExport[currentReadId] = records
        }

        return Pair(countMap, recordsToExport)
    }

    fun processReads(sampleName:String, numSampleGametes: Int, recordsForRead: List<SAMRecord>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, hapIdToSampleGamete: Map<String,List<SampleGamete>>, refRangeToIndexMap: Map<String, Int>) : Pair<AlignmentClass, List<Pair<SAMRecord?,SAMRecord?>>> {
        //Before we pair off the reads we should pull out all the truely unaligned classes
        //This should never happen but we should check for it
        if(recordsForRead.isEmpty()) {
            return Pair(AlignmentClass.UNALIGN, listOf(Pair(null,null)))
        }

        //check to see that both reads are either null or are unaligned
        if( recordsForRead.size == 2 && recordsForRead[0].readUnmappedFlag && recordsForRead[0].mateUnmappedFlag) {
            //Its unaligned
            return Pair(AlignmentClass.UNALIGN, listOf(Pair(recordsForRead[0],recordsForRead[1])))
        }

        //Then need to group by sampleGamete by hapId
        //Then need to pair off correctly
        //For the pair only keep track of the best ones based on edit distance
        val recordsGroupedByContig = filterAlignmentByPairFlag(recordsForRead)
            .groupBy { record -> hapIdToSampleGamete[record.contig]!! }
            .map { pairOffAlignments(it.value) }

        //We are looking for various edge cases
        val classification = classifyAlignments(sampleName, numSampleGametes,recordsGroupedByContig, hapIdToRefRangeMap, hapIdToSampleGamete, refRangeToIndexMap)
        return Pair(classification, recordsGroupedByContig)
    }

    fun filterAlignmentByPairFlag(records: List<SAMRecord>): List<SAMRecord> {
        return records.filter{!it.readUnmappedFlag}.groupBy { it.firstOfPairFlag }.filter{ it.value.isNotEmpty() }.map { keepBestAlignments(it.value) }.flatten()
    }

    fun pairOffAlignments(records: List<SAMRecord>): Pair<SAMRecord?,SAMRecord?> {
        //these records are all hitting the same contig.  Need to split them by first in pair and second in pair
        val bestAlignments = records.groupBy { it.firstOfPairFlag }.map { keepBestAlignment(it.value) } //Group and filter just to be safe
        if (bestAlignments.size == 1) { return Pair(bestAlignments[0], null)}
        return Pair(bestAlignments[0], bestAlignments[1])
    }

    fun keepBestAlignment(records : List<SAMRecord>) : SAMRecord? {
        return records.minByOrNull { it.getIntegerAttribute("NM") }
    }

    fun keepBestAlignments(records: List<SAMRecord>) : List<SAMRecord> {
        val bestAlignmentScore = records.minOf { it.getIntegerAttribute("NM") }
        return records.filter { it.getIntegerAttribute("NM") == bestAlignmentScore }
    }

    fun classifyAlignments(sampleName: String, numSampleGametes: Int ,records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, hapIdToSampleGamete: Map<String, List<SampleGamete>>, refRangeToIndexMap: Map<String, Int>): AlignmentClass {
        return when {
            isSingle(records) -> classifySingleAlignments(sampleName, numSampleGametes, records, hapIdToRefRangeMap, hapIdToSampleGamete, refRangeToIndexMap)
            isPaired(records) -> classifyPairedAlignments(sampleName, numSampleGametes, records, hapIdToRefRangeMap, hapIdToSampleGamete, refRangeToIndexMap)
            else -> AlignmentClass.UNALIGN
        }
    }

    fun isSingle(records: List<Pair<SAMRecord?, SAMRecord?>>): Boolean {
        //check to see if all the records are single ended
        return records.map { it.first == null || it.second == null }.all { it }
    }

    fun isPaired(records: List<Pair<SAMRecord?, SAMRecord?>>): Boolean {
        //check to see if  there is a  record are paired ended
        //If we have one we have enough, As long as we have 1 paired we ignore the rest
        return records.map { it.first != null && it.second != null }.any { it }
    }

    fun classifySingleAlignments(sampleName: String, numSampleGametes: Int, records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, hapIdToSampleGamete: Map<String, List<SampleGamete>>, refRangeToIndexMap: Map<String, Int>) :AlignmentClass {
        //Check unique first then the splits and offASM and then rare vs common
        //We cant have readSplit as these are all single ended reads
        return when{
            records.size == 1 -> AlignmentClass.SINGLEUNIQUE
            isSingleOffASM(sampleName, records, hapIdToSampleGamete) -> AlignmentClass.SINGLEOFFASM
            isSingleAlignConsec(records,hapIdToRefRangeMap, refRangeToIndexMap) -> AlignmentClass.SINGLEALIGNSPLITCONSEC //Check consec first as its a subclass of AlignSplit
            isSingleAlignSplit(records,hapIdToRefRangeMap) -> AlignmentClass.SINGLEALIGNSPLIT
            records.size in 2 .. numSampleGametes/2 -> AlignmentClass.SINGLERARE
            else -> AlignmentClass.SINGLECOMMON
        }
    }

    fun classifyPairedAlignments(sampleName: String, numSampleGametes: Int, records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, hapIdToSampleGamete: Map<String, List<SampleGamete>>, refRangeToIndexMap: Map<String, Int>) :AlignmentClass {
        //Filter out the single ended records as we have some well formed paired records
        val onlyPairedRecords = records.filter { it.first != null && it.second != null }

        assert(onlyPairedRecords.isNotEmpty()) { "Only paired records should be passed to this function" }

        //Check unique first then the splits and offASM and then rare vs common
        return when{
            isPairedAlignConsec(onlyPairedRecords,hapIdToRefRangeMap, refRangeToIndexMap) -> AlignmentClass.PAIRALIGNSPLITCONSEC //Check consec first as its a subclass of AlignSplit
            isPairedAlignSplit(onlyPairedRecords,hapIdToRefRangeMap) -> AlignmentClass.PAIRALIGNSPLIT
            isPairedOffASM(sampleName, onlyPairedRecords, hapIdToSampleGamete) -> AlignmentClass.PAIROFFASM
            isPairedReadSplitConsec(onlyPairedRecords, hapIdToRefRangeMap, refRangeToIndexMap) -> AlignmentClass.PAIRREADSPLITCONSEC
            isPairedReadSplit(onlyPairedRecords, hapIdToRefRangeMap) -> AlignmentClass.PAIRREADSPLIT
            records.size == 1 -> AlignmentClass.PAIRUNIQUE
            records.size in 2 .. numSampleGametes/2 -> AlignmentClass.PAIRRARE
            else -> AlignmentClass.PAIRCOMMON
        }
    }

    fun isSingleAlignSplit(records: List<Pair<SAMRecord?, SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>) : Boolean {
        //If one record has a different refRange than the others then it is an align split
        val refRangeSet = getRefRangesHit(records, hapIdToRefRangeMap)
        return refRangeSet.size > 1
    }

    private fun getRefRangesHit(
        records: List<Pair<SAMRecord?, SAMRecord?>>,
        hapIdToRefRangeMap: Map<String, List<ReferenceRange>>
    ): MutableSet<ReferenceRange> {
        val refRangeSet = mutableSetOf<ReferenceRange>()
        for (record in records) {
            val readOneHapId = record.first?.contig
            val readTwoHapId = record.second?.contig

            val currentRangeSet = mutableSetOf<ReferenceRange>()
            //Need to check nulls here
            if (readOneHapId != null ) currentRangeSet.addAll(hapIdToRefRangeMap[readOneHapId]!!)
            if (readTwoHapId != null ) currentRangeSet.addAll(hapIdToRefRangeMap[readTwoHapId]!!)

            refRangeSet.addAll(currentRangeSet)

        }
        return refRangeSet
    }

    private fun getRefRangesHitPaired(
        records: List<Pair<SAMRecord?, SAMRecord?>>,
        hapIdToRefRangeMap: Map<String, List<ReferenceRange>>,
        checkReadSame : Boolean
    ): MutableSet<ReferenceRange> {
        val refRangeSet = mutableSetOf<ReferenceRange>()
        for (record in records) {
            val readOneHapId = record.first?.contig
            val readTwoHapId = record.second?.contig

            //This is useful for checking if its readSplit or not
            if(checkReadSame && readOneHapId != readTwoHapId) {
                continue
            }

            val currentRangeSet = mutableSetOf<ReferenceRange>()
            //Need to check nulls here
            currentRangeSet.addAll(hapIdToRefRangeMap[readOneHapId]!!)
            currentRangeSet.addAll(hapIdToRefRangeMap[readTwoHapId]!!)

            refRangeSet.addAll(currentRangeSet)

        }
        return refRangeSet
    }


    fun isSingleAlignConsec(records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, refRangeToIndexMap: Map<String,Int>) : Boolean {
        //get all the referenceRanges
        val refRangesHit = getRefRangesHit(records, hapIdToRefRangeMap)
        //if we have more than 2 refRanges then we are not consecutive
        return when {
            (refRangesHit.size > 2) -> false
            (refRangesHit.size == 1) -> false
            (refRangesHit.size == 2) -> {
                val refRangeIndexes = refRangesHit.map { refRangeToIndexMap[it.toString()]!! }
                (refRangeIndexes[0] + 1 == refRangeIndexes[1] || refRangeIndexes[0] - 1 == refRangeIndexes[1])
            }
            else -> false
        }
    }


    fun isSingleOffASM(sampleName: String, records: List<Pair<SAMRecord?, SAMRecord?>>, hapIdToSampleGamete: Map<String, List<SampleGamete>> ) :Boolean {
        check(records.isNotEmpty()) { "Records must not be empty" }
        //loop through the records and get out the hapIds
        val samplesHit = mutableSetOf<String>()
        for (record in records) {
            val readOneHapId = record.first?.contig
            val readTwoHapId = record.second?.contig

            //Need to check nulls here
            if (readOneHapId != null) samplesHit.addAll(hapIdToSampleGamete[readOneHapId]!!.map { it.name })
            if (readTwoHapId != null) samplesHit.addAll(hapIdToSampleGamete[readTwoHapId]!!.map { it.name })
        }
        //If that set does not contain the Sample Name return true as no alignment has our sample
        return !samplesHit.contains(sampleName)
    }

    fun isPairedAlignConsec(records: List<Pair<SAMRecord?, SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, refRangeToIndexMap: Map<String, Int>) : Boolean {
        //get all the referenceRanges
        //We want to make sure that the reads are hitting the same haplotype
        val refRangesHit = getRefRangesHitPaired(records, hapIdToRefRangeMap,true)
        //if we have more than 2 refRanges then we are not consecutive
        return when {
            (refRangesHit.size > 2) -> false
            (refRangesHit.size == 1) -> false
            (refRangesHit.size == 2) -> {
                val refRangeIndexes = refRangesHit.map { refRangeToIndexMap[it.toString()]!! }
                (refRangeIndexes[0] + 1 == refRangeIndexes[1] || refRangeIndexes[0] - 1 == refRangeIndexes[1])
            }
            else -> false
        }
    }

    fun isPairedAlignSplit(records: List<Pair<SAMRecord?, SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>) : Boolean {
        //If one record has a different refRange than the others then it is an align split
        val refRangeSet = getRefRangesHitPaired(records, hapIdToRefRangeMap,true)
        return refRangeSet.size > 1
    }

    fun isPairedOffASM(sampleName: String, records: List<Pair<SAMRecord?, SAMRecord?>>, hapIdToSampleGamete: Map<String, List<SampleGamete>> ) :Boolean {
        check(records.isNotEmpty()) { "Records must not be empty" }
        //loop through the records and get out the hapIds
        val samplesHit = mutableSetOf<String>()
        for (record in records) {
            val readOneHapId = record.first?.contig
            val readTwoHapId = record.second?.contig

            //Need to check nulls here
            samplesHit.addAll(hapIdToSampleGamete[readOneHapId]!!.map { it.name })
            samplesHit.addAll(hapIdToSampleGamete[readTwoHapId]!!.map { it.name })
        }
        //If that set does not contain the Sample Name return true as no alignment has our sample
        return !samplesHit.contains(sampleName)
    }

    fun isPairedReadSplitConsec(records: List<Pair<SAMRecord?, SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, refRangeToIndexMap: Map<String, Int>) : Boolean {
        //Check to see that the reads are not all hitting same haplotype as their pair
        //If the reads are not split they cannot be split consecutively
        if (!isPairedReadSplit(records, hapIdToRefRangeMap)) {
            return false
        }
        //get all the referenceRanges
        val refRangesHit = getRefRangesHitPaired(records, hapIdToRefRangeMap,false)
        //if we have more than 2 refRanges then we are not consecutive
        return when {
            (refRangesHit.size > 2) -> false
            (refRangesHit.size == 1) -> false
            (refRangesHit.size == 2) -> {
                val refRangeIndexes = refRangesHit.map { refRangeToIndexMap[it.toString()]!! }
                (refRangeIndexes[0] + 1 == refRangeIndexes[1] || refRangeIndexes[0] - 1 == refRangeIndexes[1])
            }
            else -> false
        }
    }

    fun isPairedReadSplit(records: List<Pair<SAMRecord?, SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>) : Boolean {
        //If a pair of records hit different haplotypes then it is a read split
        for(recordPair in records) {
            if(recordPair.first?.contig != recordPair.second?.contig) {
                return true
            }
        }
        return false
    }




    private fun writeFastq(fastqFileName : String, fastq : List<SAMRecord>, order : String) {
        bufferedWriter("$order$fastqFileName").use { writer ->
            for (record in fastq) {
                writer.write("@${record.readName}\n")
                writer.write("${record.readString}\n")
                writer.write("+\n")
                writer.write("${record.baseQualityString}\n")
            }
        }
    }

    /**
     * Function that exports the filtered down fastq files and a table of [readID] -> [HapID hits] to outputFileName
     */
    fun outputReadsAndHapIdSets(
        tableFileName: String,
        fastqFileName: String,
        // key = readId, values = alignments paired off by the correct strand that we want to export
        alignments: Map<String, List<Pair<SAMRecord?, SAMRecord?>>>
    ) {

        val fastq1 : MutableList<SAMRecord> = mutableListOf()
        val fastq2 : MutableList<SAMRecord> = mutableListOf()

        // Write the table to the file as tab-delimited text
        bufferedWriter(tableFileName).use { writer ->
            writer.write("readID\tHapIDHits\n")
            for ((readID, samlist) in alignments) {
                // add all HapIDs to hits
                val hits : MutableSet<String> = mutableSetOf<String>()
                for (pair in samlist) {
                    // if pair has base quality string, add to the fastq file list
                    if(pair.first != null && pair.first!!.baseQualityString.isNotEmpty()) {
                        fastq1.add(pair.first!!)
                        // SAMRecord.contig gets HapID
                        hits.add(pair.first!!.contig)
                    }
                    if(pair.second != null && pair.second!!.baseQualityString.isNotEmpty()) {
                        fastq2.add(pair.second!!)
                        hits.add(pair.second!!.contig)
                    }
                }
                writer.write("$readID\t${hits.joinToString(separator = ", ")}\n")
            }
        }
        // write fastq files
        writeFastq(fastqFileName, fastq1, "1")
        writeFastq(fastqFileName, fastq2, "2")
    }

    fun printCountMap(countMap: Map<AlignmentClass, Int>) {
        myLogger.info("Number of reads per class:")
        for ((key, value) in countMap) {
            myLogger.info("$key: $value")
        }
    }
}