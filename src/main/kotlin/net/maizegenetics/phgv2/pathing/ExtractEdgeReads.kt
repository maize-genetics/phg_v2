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
import kotlin.math.abs

enum class AlignmentClass {
    PAIRUNIQUE, PAIRRARE, PAIRCOMMON, PAIRREADSPLIT, PAIRREADSPLITCONSEC, PAIRALIGNSPLIT, PAIRALIGNSPLITCONSEC, PAIROFFASM, UNALIGN,
    SINGLEPARTIALUNIQUE, SINGLEPARTIALRARE, SINGLEPARTIALCOMMON, SINGLEPARTIALALIGNSPLIT, SINGLEPARTIALALIGNSPLITCONSEC, SINGLEPARTIALOFFASM,
    SINGLEUNIQUE, SINGLERARE, SINGLECOMMON, SINGLEALIGNSPLIT, SINGLEALIGNSPLITCONSEC,SINGLEOFFASM

}

/**
 * This class is to find edge case reads for use as integration tests for the read mapping pipeline
 *
 * IT IS CURRENTLY A WIP AND IS UNDER ACTIVE DEVELOPMENT
 */
class ExtractEdgeReads : CliktCommand( help = "Extract out Edge Case reads from SAMs/BAMs") {
    val bamFile by option(help = "BAM/SAM file to use")
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

        myLogger.info("Extracting reads from $bamFile")
        val (countMap, recordsToExport, readToClassificaiton) = extractReads(sampleName, bamFile, graph, maxClassNum)
        //Print out the countMap
        printCountMap(countMap)
        val outputFileName = "${outputFileDir}/${File(bamFile).nameWithoutExtension}"
        myLogger.info("Writing output to $outputFileName")
        outputReadsAndHapIdSets("${outputFileName}_table.txt", "${outputFileName}", recordsToExport, readToClassificaiton)
    }

    fun extractReads(sampleName: String, bamFile: String, graph: HaplotypeGraph, maxClassNum : Int = 10) : Triple<Map<AlignmentClass,Int>, Map<String, List<Pair<SAMRecord?,SAMRecord?>>>, Map<String,AlignmentClass>> {
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
        val readToClassification = mutableMapOf<String, AlignmentClass>()
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
                if(countMap[classification]!! <= maxClassNum) {
                    recordsToExport[currentReadId] = records
                    readToClassification[currentReadId] = classification
                }

                readCounter++
                if(readCounter % 100000 == 0) {
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
        if(countMap[classification]!! <= maxClassNum) {
            recordsToExport[currentReadId] = records
            readToClassification[currentReadId] = classification
        }

        return Triple(countMap, recordsToExport, readToClassification)
    }

    fun processReads(sampleName:String, numSampleGametes: Int, recordsForRead: List<SAMRecord>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, hapIdToSampleGamete: Map<String,List<SampleGamete>>, refRangeToIndexMap: Map<String, Int>) : Pair<AlignmentClass, List<Pair<SAMRecord?,SAMRecord?>>> {

        val primaryAlignments = listOf(getPrimaryAlignments(recordsForRead))


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
        val bestRecordsForAlignment = filterAlignmentByPairFlag(recordsForRead)
        if(bestRecordsForAlignment.isEmpty()) {
            return Pair(AlignmentClass.UNALIGN, primaryAlignments)
        }

        //Then need to group by sampleGamete by hapId
        //Then need to pair off correctly
        //For the pair only keep track of the best ones based on edit distance
        val recordsGroupedByContig = bestRecordsForAlignment
            .groupBy { record -> hapIdToSampleGamete[record.contig]!! }
            .map { pairOffAlignmentsByHapId(it.value) }
            .map { pairOffAlignmentsBySample(it,hapIdToRefRangeMap,hapIdToSampleGamete, refRangeToIndexMap) }
            .flatten()

        //We are looking for various edge cases
        val classification = classifyAlignments(sampleName, numSampleGametes,recordsGroupedByContig, hapIdToRefRangeMap, hapIdToSampleGamete, refRangeToIndexMap)


        //we need to add in the primary alignments here as we need to make sure that the sequence is included
        //The primary alignment might be missing if it was not paired off correctly.
        return Pair(classification, recordsGroupedByContig + primaryAlignments)
    }

    fun filterAlignmentByPairFlag(records: List<SAMRecord>): List<SAMRecord> {
        return records.filter{!it.readUnmappedFlag}.filter { !it.cigar.isClipped }.groupBy { it.firstOfPairFlag }.filter{ it.value.isNotEmpty() }.map { keepBestAlignments(it.value) }.flatten()
    }

    fun pairOffAlignmentsByHapId(records: List<SAMRecord>): List<Pair<SAMRecord?,SAMRecord?>> {
        //these records are all hitting the same contig.  Need to split them by first in pair and second in pair
        val bestAlignments = records.groupBy { it.firstOfPairFlag }.map { Pair(it.key,keepBestAlignments(it.value)) }.toMap()//Group and filter just to be safe

        //Go through best alignments and try to match them up by hapId
        val trueHaps = bestAlignments[true]?.groupBy { it.contig } ?: mapOf()
        val falseHaps = bestAlignments[false]?.groupBy { it.contig } ?: mapOf()

        //loop through trueHaps and try to match them up with falseHaps
        val pairedBestAlignments = mutableListOf<Pair<SAMRecord?,SAMRecord?>>()
        buildFirstPairs(trueHaps, falseHaps, pairedBestAlignments)
        buildSecondPairs(trueHaps, falseHaps, pairedBestAlignments)

        return pairedBestAlignments
    }

    fun pairOffAlignmentsBySample(records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>,
                                  hapIdToSampleGamete: Map<String, List<SampleGamete>>, refRangeToIndexMap: Map<String, Int> ): List<Pair<SAMRecord?,SAMRecord?>> {

        val allRecords = records.map { listOf(it.first,it.second) }.flatten().filterNotNull()

        val sampleGametesToRecords = mutableMapOf<SampleGamete,MutableList<SAMRecord>>()

        for(record in allRecords) {
            if(record.contig !in hapIdToSampleGamete.keys) {
                continue
            }

            val sampleGametes = hapIdToSampleGamete[record.contig]
            for(sampleGamete in sampleGametes!!) {
                if (sampleGametesToRecords.containsKey(sampleGamete)) {
                    sampleGametesToRecords[sampleGamete]!!.add(record)
                } else {
                    sampleGametesToRecords[sampleGamete] = mutableListOf(record)
                }
            }
        }

        val pairedRecords = mutableListOf<Pair<SAMRecord?,SAMRecord?>>()

        for(sampleGamete in sampleGametesToRecords.keys) {
            val recordsFromGamete = sampleGametesToRecords[sampleGamete]!!

            //filter out paired records and add to list
            val pairedByHapId = pairOffAlignmentsByHapId(recordsFromGamete)
            val onlyPaired = pairedByHapId.filter { isPaired(listOf(it)) }
            if(onlyPaired.isNotEmpty()) {
                pairedRecords.addAll(onlyPaired)
                continue
            }

            val remainingRecords = pairedByHapId.filter { !isPaired(listOf(it)) }
            if(remainingRecords.size == 1) {
                pairedRecords.add(remainingRecords[0])
            }
            else {
                val firstHits = remainingRecords.filter { it.first != null }
                val secondHits = remainingRecords.filter { it.second != null }

                //if one of the two is empty add the other
                if(firstHits.isEmpty()) {
                    pairedRecords.add(secondHits[0])
                }
                else if(secondHits.isEmpty()) {
                    pairedRecords.add(firstHits[0])
                }
                else {
                    //If we have both we need to decide how to pair them off
                    var pairedOff = false
                    for(firstHit in firstHits) {
                        val firstRefRange = hapIdToRefRangeMap[firstHit.first!!.contig]!!
                        val firstRefRangeIdx = refRangeToIndexMap[firstRefRange[0].toString()]!!
                        for(secondHit in secondHits) {
                            val secondRefRange = hapIdToRefRangeMap[secondHit.second!!.contig]!!
                            val secondRefRangeIdx = refRangeToIndexMap[secondRefRange[0].toString()]!!
                            if(abs(firstRefRangeIdx - secondRefRangeIdx) == 1) {
                                pairedRecords.add(Pair(firstHit.first, secondHit.second))
                                pairedOff = true
                                break
                            }
                        }
                    }
                    if(!pairedOff) {
                        //If we have not paired off any of the records we just add the first one from each
                        pairedRecords.add(Pair(firstHits[0].first, secondHits[0].second))
                    }
                }
            }
        }
        return pairedRecords
    }

    fun buildFirstPairs(trueHaps: Map<String, List<SAMRecord>>, falseHaps: Map<String, List<SAMRecord>>, pairedBestAlignments: MutableList<Pair<SAMRecord?,SAMRecord?>>) {
        for(key in trueHaps.keys) {
            val trueRecords = trueHaps[key]!!
            val falseRecords = falseHaps[key] ?: listOf()
            //Find the primary alignments in both the true and false records
            val truePrimaryOrNull = getPrimaryAlignments(trueRecords)
            val falsePrimaryOrNull = getPrimaryAlignments(falseRecords)

            val trueAlignment = if(truePrimaryOrNull == Pair(null,null)) {
                trueRecords.first()
            }
            else {
                truePrimaryOrNull.first
            }

            val falseAlignment = if(falseRecords.isEmpty()) {
                null
            }
            else if(falsePrimaryOrNull == Pair(null,null)) {
                falseRecords.first()
            }
            else {
                falsePrimaryOrNull.first
            }

            if(pairedBestAlignments.contains(Pair(trueAlignment, falseAlignment))) {
                continue
            }
            //If we have a list of true records it does not matter we just need one of them
            pairedBestAlignments.add(Pair(trueAlignment, falseAlignment))
        }
    }

    fun buildSecondPairs(trueHaps: Map<String, List<SAMRecord>>, falseHaps: Map<String, List<SAMRecord>>, pairedBestAlignments: MutableList<Pair<SAMRecord?,SAMRecord?>>) {
        for(key in falseHaps.keys) {
            val trueRecords = trueHaps[key] ?: listOf()
            val falseRecords = falseHaps[key]!!
            //Find the primary alignments in both the true and false records
            val truePrimaryOrNull = getPrimaryAlignments(trueRecords)
            val falsePrimaryOrNull = getPrimaryAlignments(falseRecords)

            val trueAlignment = if(trueRecords.isEmpty()) {
                null
            }
            else if(truePrimaryOrNull == Pair(null,null)) {
                trueRecords.first()
            }
            else {
                truePrimaryOrNull.first
            }

            val falseAlignment = if(falsePrimaryOrNull == Pair(null,null)) {
                falseRecords.first()
            }
            else {
                falsePrimaryOrNull.first
            }

            if(pairedBestAlignments.contains(Pair(trueAlignment, falseAlignment))) {
                continue
            }

            //If we have a list of true records it does not matter we just need one of them
            pairedBestAlignments.add(Pair(trueAlignment, falseAlignment))
        }
    }

    /**
     * Function to find the best alignment from the list of SAMRecords
     */
    fun getBestAlignment(records : List<SAMRecord>) : SAMRecord? {
        return records.minByOrNull { it.getIntegerAttribute("NM") }
    }

    /**
     * Function to filter out the bad alignments based on edit distance
     */
    fun keepBestAlignments(records: List<SAMRecord>) : List<SAMRecord> {
        val bestAlignmentScore = records.minOf { it.getIntegerAttribute("NM") }
        return records.filter { it.getIntegerAttribute("NM") == bestAlignmentScore }
    }

    /**
     * Function to get the primary alignments from the list of SAMRecords
     * TODO Unit test this fully
     */
    fun getPrimaryAlignments(records: List<SAMRecord>) : Pair<SAMRecord?,SAMRecord?> {
        val primaryAlignments = records.filter { !it.isSecondaryOrSupplementary }

        return when {
            primaryAlignments.isEmpty() -> Pair(null,null)
            primaryAlignments.size == 1 -> Pair(primaryAlignments[0], null)
            else -> {
                val firstInPair = primaryAlignments.filter { it.firstOfPairFlag }
                val secondInPair = primaryAlignments.filter { !it.firstOfPairFlag }
                Pair(firstInPair.first(), secondInPair.first())
            }
        }
    }

    /**
     * Function to classify the alignment into one of the various classes
     */
    fun classifyAlignments(sampleName: String, numSampleGametes: Int ,records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, hapIdToSampleGamete: Map<String, List<SampleGamete>>, refRangeToIndexMap: Map<String, Int>): AlignmentClass {
        return when {
            isSingle(records) -> classifySingleAlignments(sampleName, numSampleGametes, records, hapIdToRefRangeMap, hapIdToSampleGamete, refRangeToIndexMap)
            isPaired(records) -> classifyPairedAlignments(sampleName, numSampleGametes, records, hapIdToRefRangeMap, hapIdToSampleGamete, refRangeToIndexMap)
            else -> AlignmentClass.UNALIGN
        }
    }


    /**
     * Function to see if the [records] are unable to be paired off.
     * If a single set of alignments can be paired this will be returned false
     */
    fun isSingle(records: List<Pair<SAMRecord?, SAMRecord?>>): Boolean {
        //check to see if all the records are single ended
        return records.map { it.first == null || it.second == null }.all { it }
    }

    /**
     * Function to see if  any of the [records] are paired off.
     * We will filter out the single ended reads later on if this returns true.
     */
    fun isPaired(records: List<Pair<SAMRecord?, SAMRecord?>>): Boolean {
        //check to see if  there is a  record are paired ended
        //If we have one we have enough, As long as we have 1 paired we ignore the rest
        return records.map { it.first != null && it.second != null }.any { it }
    }

    /**
     * Function to classify the reads into the various alignment classes.  This will only work for non-paired reads
     */
    fun classifySingleAlignments(sampleName: String, numSampleGametes: Int, records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, hapIdToSampleGamete: Map<String, List<SampleGamete>>, refRangeToIndexMap: Map<String, Int>) :AlignmentClass {
        //Check offASM then the splits  and then unique rare and common
        //We cant have readSplit as these are all single ended reads
        return when{
            isSingleOffASM(sampleName, records, hapIdToSampleGamete) -> AlignmentClass.SINGLEOFFASM
            isSingleAlignConsec(records,hapIdToRefRangeMap, refRangeToIndexMap) -> AlignmentClass.SINGLEALIGNSPLITCONSEC //Check consec first as its a subclass of AlignSplit
            isSingleAlignSplit(records,hapIdToRefRangeMap) -> AlignmentClass.SINGLEALIGNSPLIT
            records.size == 1 -> AlignmentClass.SINGLEUNIQUE
            records.size in 2 .. numSampleGametes/2 -> AlignmentClass.SINGLERARE
            else -> AlignmentClass.SINGLECOMMON
        }
    }

    /**
     * Function to classify the paired reads into the various alignment classes.  This will first filter out the single ended reads.
     * If that list is empty it will throw an error.
     */
    fun classifyPairedAlignments(sampleName: String, numSampleGametes: Int, records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, hapIdToSampleGamete: Map<String, List<SampleGamete>>, refRangeToIndexMap: Map<String, Int>) :AlignmentClass {
        //Filter out the single ended records as we have some well formed paired records
        val onlyPairedRecords = records.filter { it.first != null && it.second != null }

        assert(onlyPairedRecords.isNotEmpty()) { "Only paired records should be passed to this function" }

        //Check the Split and off asm first as they could still be that way but also common/unique
        return when{
            isPairedOffASM(sampleName, onlyPairedRecords, hapIdToSampleGamete) -> AlignmentClass.PAIROFFASM
            isPairedAlignConsec(onlyPairedRecords,hapIdToRefRangeMap, refRangeToIndexMap) -> AlignmentClass.PAIRALIGNSPLITCONSEC //Check consec first as its a subclass of AlignSplit
            isPairedAlignSplit(onlyPairedRecords,hapIdToRefRangeMap) -> AlignmentClass.PAIRALIGNSPLIT
            isPairedReadSplitConsec(onlyPairedRecords, hapIdToRefRangeMap, refRangeToIndexMap) -> AlignmentClass.PAIRREADSPLITCONSEC
            isPairedReadSplit(onlyPairedRecords) -> AlignmentClass.PAIRREADSPLIT
            records.size == 1 -> AlignmentClass.PAIRUNIQUE
            records.size in 2 .. numSampleGametes/2 -> AlignmentClass.PAIRRARE
            else -> AlignmentClass.PAIRCOMMON
        }
    }

    /**
     * Function to check the maps to see which reference ranges got hit by the alignments.
     * This is designed to work for single ended reads
     */
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

    /**
     * Function to get the set of reference ranges which are hit by the paired alignments
     */
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

    /**
     * Function to check to see if we have an alignment split - Some alignments hit one ref range and others hit a different one
     */
    fun isSingleAlignSplit(records: List<Pair<SAMRecord?, SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>) : Boolean {
        //If one record has a different refRange than the others then it is an align split
        val refRangeSet = getRefRangesHit(records, hapIdToRefRangeMap)
        return refRangeSet.size > 1
    }


    /**
     * Function to check to see if we have an alignment split but the reference ranges are consecutively next to each other
     */
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

    /**
     * Function to determine if we have only unpaired non-target alignment hits.
     * We are looking to see if none of the alignments hit the target [sampleName]
     */
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

    /**
     * Function to see if we have paired alignments but they are splitting multiple reference ranges but those reference
     * ranges are consecutive
     */
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

    /**
     * Function to check to see if our paired alignments are split over multiple reference ranges but the two reads in the pair hit the same haplotype
     */
    fun isPairedAlignSplit(records: List<Pair<SAMRecord?, SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>) : Boolean {
        //If one record has a different refRange than the others then it is an align split
        val refRangeSet = getRefRangesHitPaired(records, hapIdToRefRangeMap,true)
        return refRangeSet.size > 1
    }

    /**
     * Function to check to see if we have paired reads but none of them hit the target sample
     */
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

    /**
     * Function to check to see if the paired alignments hit consecutive reference ranges, and the reads themselves are split over the consecutive haplotypes
     */
    fun isPairedReadSplitConsec(records: List<Pair<SAMRecord?, SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, refRangeToIndexMap: Map<String, Int>) : Boolean {
        //Check to see that the reads are not all hitting same haplotype as their pair
        //If the reads are not split they cannot be split consecutively
        if (!isPairedReadSplit(records)) {
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

    /**
     * Function to check to see if the paired alignments hit different reference ranges and the reads themselves are split over the different haplotypes
     */
    fun isPairedReadSplit(records: List<Pair<SAMRecord?, SAMRecord?>>) : Boolean {
        //If a pair of records hit different haplotypes then it is a read split
        for(recordPair in records) {
            if(recordPair.first?.contig != recordPair.second?.contig) {
                return true
            }
        }
        return false
    }


    /**
     * Function to write out fastq files based on the [fastq] list
     */
    private fun writeFastq(fastqFileName : String, fastq : List<SAMRecord>, readNumber : String) {
        bufferedWriter("${fastqFileName}_$readNumber.fastq").use { writer ->
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
        alignments: Map<String, List<Pair<SAMRecord?, SAMRecord?>>>,
        readToClassification: Map<String, AlignmentClass>
    ) {

        val fastq1 = mutableListOf<SAMRecord>()
        val fastq2 = mutableListOf<SAMRecord>()
        val readProcessedFastq1 = mutableSetOf<String>()
        val readProcessedFastq2 = mutableSetOf<String>()

        // Write the table to the file as tab-delimited text
        bufferedWriter(tableFileName).use { writer ->
            writer.write("readID\tClass\tHapIDHits\n")
            for ((readID, samlist) in alignments) {
                // add all HapIDs to hits
                val hits : MutableSet<String> = mutableSetOf<String>()
                for (pair in samlist) {
                    // if pair has base quality string, add to the fastq file list
                    if(pair.first != null) {
                        // SAMRecord.contig gets HapID
                        if(!pair.first!!.readUnmappedFlag) { //The first read might be unmapped which means we would not have a contig hit
                            hits.add(pair.first!!.contig)
                        }
                        val baseQualString = pair.first!!.baseQualityString
                        if(baseQualString.isNotEmpty() && baseQualString != "*" && !readProcessedFastq1.contains(readID)) {
                            fastq1.add(pair.first!!)
                            readProcessedFastq1.add(readID)
                        }
                    }
                    if(pair.second != null) {
                        if(!pair.second!!.readUnmappedFlag) { //The second read might be unmapped which means we would not have a contig hit
                            hits.add(pair.second!!.contig)
                        }

                        val baseQualString = pair.second!!.baseQualityString
                        if(baseQualString.isNotEmpty() && baseQualString != "*" && !readProcessedFastq2.contains(readID)) {
                            fastq2.add(pair.second!!)
                            readProcessedFastq2.add(readID)
                        }
                    }
                }
                writer.write("$readID\t${readToClassification[readID]}\t${hits.joinToString(separator = ",")}\n")
            }
        }

        // write fastq files
        writeFastq(fastqFileName, fastq1, "0")
        writeFastq(fastqFileName, fastq2, "1")

        for(readID in readProcessedFastq1) {
            if(!readProcessedFastq2.contains(readID)) {
                println("Read $readID was only processed in the first read")
            }
        }
        for(readID in readProcessedFastq2) {
            if(!readProcessedFastq1.contains(readID)) {
                println("Read $readID was only processed in the second read")
            }
        }
    }


    /**
     * Function to print out the total count map to the logger.
     */
    fun printCountMap(countMap: Map<AlignmentClass, Int>) {
        myLogger.info("Number of reads per class:")
        for ((key, value) in countMap) {
            myLogger.info("$key: $value")
        }
    }
}