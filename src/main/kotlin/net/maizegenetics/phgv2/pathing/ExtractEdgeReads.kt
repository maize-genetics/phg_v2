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
                if(countMap[classification]!! < maxClassNum) {
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
        if(countMap[classification]!! < maxClassNum) {
            recordsToExport[currentReadId] = records
            readToClassification[currentReadId] = classification
        }

        return Triple(countMap, recordsToExport, readToClassification)
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
            .flatten()

        //We are looking for various edge cases
        val classification = classifyAlignments(sampleName, numSampleGametes,recordsGroupedByContig, hapIdToRefRangeMap, hapIdToSampleGamete, refRangeToIndexMap)

        val primaryAlignments = listOf(getPrimaryAlignments(recordsForRead))

//        return Pair(classification, recordsGroupedByContig)
        return Pair(classification, primaryAlignments)
    }

    fun filterAlignmentByPairFlag(records: List<SAMRecord>): List<SAMRecord> {
        return records.filter{!it.readUnmappedFlag}.groupBy { it.firstOfPairFlag }.filter{ it.value.isNotEmpty() }.map { keepBestAlignments(it.value) }.flatten()
    }

    fun pairOffAlignments(records: List<SAMRecord>): List<Pair<SAMRecord?,SAMRecord?>> {
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
//        if (bestAlignments.size == 1) { return Pair(bestAlignments[0], null)}
//        return Pair(bestAlignments[0], bestAlignments[1])
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

            //If we have a list of true records it does not matter we just need one of them
            pairedBestAlignments.add(Pair(trueAlignment, falseAlignment))
        }
    }

    fun keepBestAlignment(records : List<SAMRecord>) : SAMRecord? {
        return records.minByOrNull { it.getIntegerAttribute("NM") }
    }

    fun keepBestAlignments(records: List<SAMRecord>) : List<SAMRecord> {
        val bestAlignmentScore = records.minOf { it.getIntegerAttribute("NM") }
        return records.filter { it.getIntegerAttribute("NM") == bestAlignmentScore }
    }

    fun getPrimaryAlignments(records: List<SAMRecord>) : Pair<SAMRecord?,SAMRecord?> {
        val primaryAlignments = records.filter { !it.isSecondaryOrSupplementary }

        return when {
            primaryAlignments.isEmpty() -> Pair(null,null)
            primaryAlignments.size == 1 -> Pair(primaryAlignments[0], null)
            else -> {
                val firstInPair = primaryAlignments.filter { it.firstOfPairFlag }
                val secondInPair = primaryAlignments.filter { !it.firstOfPairFlag }
//                Pair(primaryAlignments[0], primaryAlignments[1])
                Pair(firstInPair.first(), secondInPair.first())
            }
        }
    }

    fun classifyAlignments(sampleName: String, numSampleGametes: Int ,records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, hapIdToSampleGamete: Map<String, List<SampleGamete>>, refRangeToIndexMap: Map<String, Int>): AlignmentClass {
        if(records.first().first?.readName == "ST-E00317:129:HVMFTCCXX:7:1101:11160:1309" || records.first().second?.readName == "ST-E00317:129:HVMFTCCXX:7:1101:11160:1309") {
            println("Here")
        }

        return when {
            isPartialSingle(records) -> classifySingleAlignments(true, sampleName, numSampleGametes, records, hapIdToRefRangeMap, hapIdToSampleGamete, refRangeToIndexMap)
            isSingle(records) -> classifySingleAlignments(false, sampleName, numSampleGametes, records, hapIdToRefRangeMap, hapIdToSampleGamete, refRangeToIndexMap)
            isPaired(records) -> classifyPairedAlignments(sampleName, numSampleGametes, records, hapIdToRefRangeMap, hapIdToSampleGamete, refRangeToIndexMap)
            else -> AlignmentClass.UNALIGN
        }
    }

    /**
     * Function to see if we have one of the reads hitting a set of haplotypes and the other hitting a different set
     * of haplotypes where neither of the two sets share sample names.  We are calling this a partial alignment as it is
     * different from a normal single alignment where only one read is aligning.
     */
    fun isPartialSingle(records: List<Pair<SAMRecord?, SAMRecord?>>): Boolean {
        if (isSingle(records)) {
            //Split the records first or second in pair
            val firstInPair = records.filter { it.first != null }
            val secondInPair = records.filter { it.second != null }
            return !(firstInPair.isEmpty() || secondInPair.isEmpty()) //If its single and each of the 2 reads have single alignments but not just one of the reads it is a partial
        }
        else {
            return false
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

    fun classifySingleAlignments(isPartial: Boolean, sampleName: String, numSampleGametes: Int, records: List<Pair<SAMRecord?,SAMRecord?>>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, hapIdToSampleGamete: Map<String, List<SampleGamete>>, refRangeToIndexMap: Map<String, Int>) :AlignmentClass {
        //Check unique first then the splits and offASM and then rare vs common
        //We cant have readSplit as these are all single ended reads
        return when{
            !isPartial && records.size == 1 -> AlignmentClass.SINGLEUNIQUE
            !isPartial && isSingleOffASM(sampleName, records, hapIdToSampleGamete) -> AlignmentClass.SINGLEOFFASM
            !isPartial && isSingleAlignConsec(records,hapIdToRefRangeMap, refRangeToIndexMap) -> AlignmentClass.SINGLEALIGNSPLITCONSEC //Check consec first as its a subclass of AlignSplit
            !isPartial && isSingleAlignSplit(records,hapIdToRefRangeMap) -> AlignmentClass.SINGLEALIGNSPLIT
            !isPartial && records.size in 2 .. numSampleGametes/2 -> AlignmentClass.SINGLERARE
            !isPartial -> AlignmentClass.SINGLECOMMON
            isPartial && records.size == 1 -> AlignmentClass.SINGLEPARTIALUNIQUE
            isPartial && isSingleOffASM(sampleName, records, hapIdToSampleGamete) -> AlignmentClass.SINGLEPARTIALOFFASM
            isPartial && isSingleAlignConsec(records,hapIdToRefRangeMap, refRangeToIndexMap) -> AlignmentClass.SINGLEPARTIALALIGNSPLITCONSEC //Check consec first as its a subclass of AlignSplit
            isPartial && isSingleAlignSplit(records,hapIdToRefRangeMap) -> AlignmentClass.SINGLEPARTIALALIGNSPLIT
            isPartial && records.size in 2 .. numSampleGametes/2 -> AlignmentClass.SINGLEPARTIALRARE
            else -> AlignmentClass.SINGLEPARTIALCOMMON

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

    fun printCountMap(countMap: Map<AlignmentClass, Int>) {
        myLogger.info("Number of reads per class:")
        for ((key, value) in countMap) {
            myLogger.info("$key: $value")
        }
    }
}