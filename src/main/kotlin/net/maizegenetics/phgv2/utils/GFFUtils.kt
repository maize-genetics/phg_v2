package net.maizegenetics.phgv2.utils

import biokotlin.util.bufferedReader
import com.google.common.collect.Range
import com.google.common.collect.RangeSet
import com.google.common.collect.TreeRangeSet
import htsjdk.tribble.AbstractFeatureReader
import htsjdk.tribble.TribbleException
import htsjdk.tribble.annotation.Strand
import htsjdk.tribble.gff.*
import htsjdk.variant.vcf.VCFFileReader
import java.io.IOException
import java.nio.file.Paths
import net.maizegenetics.phgv2.api.ReferenceRange
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.*
import kotlin.Comparator
import kotlin.collections.ArrayList
import kotlin.collections.HashMap
import kotlin.collections.HashSet

private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.Utils.GFFUtils")
/**
 * GFFUtils class:  Contains functions used to process GFF3 files.
 * These files are processed using  htsjdk GFF3 functions
 *
 * Restrictions:  The gff files read by the htsjdk reader must end
 * with .gff3, .gff3.gz, .gff or .gff.gz
 *
 * @author lcj34
 */


/**
 * Function to use htsjdk to read gff into memory.
 * Returns a Set<Gff3Feature> .
 *
 * NOTE:  gff files must end with .gff3, .gff, .gff3.gz or .gff.gz
 *        Any other extension causes htsjdk feature reader to throw an exception
 *
 * The returned feature set may be printed using htsdjk Gff3Writer
 */
fun readGFFtoGff3Feature(gffFile: String):List<Gff3Feature> {

    val gffFeatures: MutableSet<Gff3Feature> = HashSet<Gff3Feature>()
    AbstractFeatureReader.getFeatureReader(gffFile, null, Gff3Codec(), false).use{ reader ->
        for (feature in reader.iterator()) {
            // This grabs a line. It skips over commented lines
            gffFeatures.add(feature)
        }
    }

    // While gff files are generally sorted by chr, they are not always sorted by start/end.
    // The start/end postions start over when they are in a second transcript.
    // For our processing, in order to get all overlap efficiently, these need have all entries
    // sorted by chr/start/end.  This means T0002 will be interspered with T0001.

    val sortedGffEntries = gffFeatures.toList().sortedWith(Comparator { first: Gff3Feature, second: Gff3Feature ->
        var contigCompare = first.contig.compareTo(second.contig)
        if (contigCompare == 0) {
            var startCompare = first.start.compareTo(second.start)
            if (startCompare == 0) {
                first.end.compareTo(second.end)
            } else {
                first.start.compareTo(second.start)
            }
        } else {
            first.contig.compareTo(second.contig)
        }
    })
    return sortedGffEntries
}

/**
 * Function takes a key file with columns taxon, gffFile.
 * For each taxon, read the gff file using htsjdk gff feature reader.
 * Return a map of taxon to associated Gff3Feature entries
 */
fun loadGFFsToGff3Feature(keyFile: String): Map<String, TreeMap<Position,ArrayList<Gff3Feature>>> {

    val resultsTree:MutableMap<String, TreeMap<Position,ArrayList<Gff3Feature>>> = mutableMapOf()
    //val taxonToGFFfileMap:MutableMap<String, String> = HashMap<String, String>()
    val taxonToGFFfileMap:MutableMap<String, String> = mutableMapOf()
    bufferedReader(keyFile)
            .forEachLine{
                val parsedLine = it.split("\t");
                taxonToGFFfileMap.put(parsedLine[0], parsedLine[1])
            }

    // process each gff file
    for ((taxon, gffFile) in taxonToGFFfileMap) {
        val time = System.nanoTime()
        val features = readGFFtoGff3Feature(gffFile)
        val featureTreeMap = createTreeMapFromFeaturesCenter(features)
        resultsTree.put(taxon, featureTreeMap)
        val endTime = (System.nanoTime() - time)/1e9
        myLogger.info("loadGffsToGff3Feature: time to load ${taxon} gff file: ${endTime}")

    }
    return resultsTree
}

/**
 * Reads a key file with columns for "taxon" and "Path/name of gff file"
 * Returns a map of taxon->fileName
 */
fun getTaxonToGffFileMap(keyFile:String): Map<String,String>  {
    val taxonToGFFfileMap:MutableMap<String, String> = HashMap<String, String>()
    bufferedReader(keyFile)
            .forEachLine{
                val parsedLine = it.split("\t")
                taxonToGFFfileMap.put(parsedLine[0], parsedLine[1])
            }
    return taxonToGFFfileMap
}
/**
 * This method creates a mapping of the feature center position (posSTart + posEnd)/2 to
 * list of Gff3Features
 */
fun createTreeMapFromFeaturesCenter(features: List<Gff3Feature>): TreeMap<Position,ArrayList<Gff3Feature>> {
    val treeMap = TreeMap<Position,ArrayList<Gff3Feature>>()

    for (feature in features) {

        // Don't include the "chromosome" line.  This will overlap with everything for that chromosome
        if (feature.type.equals("chromosome")) continue
        // for gffs, the start is always <= end per the specs
        val center = (feature.end + feature.start)/2
        val centerPos = Position(feature.contig,center)
        var posList = treeMap.get(centerPos)
        if (posList == null) {
                posList = ArrayList<Gff3Feature>()
        }
        posList.add(feature)
        treeMap.put(centerPos,posList)
    }

    return treeMap
}

/**
 * This method takes a path ( a list of integer haplotype ids),
 *   a phg graph that is based on the haplotypeIds in the path,
 *   and an optional output file name.
 *
 *  From the graph it  pulls the asm contig and coordinates for each
 *  item on the path, finds regions in the Gff3 entries that overlap
 *  with the graph haplotype asm coordinate entries, and creates a set
 *  of Gff3Feature entries.
 *
 *  If an outputFile is specified, the gff file is written to the specified path.
 *
 * Return:  A set of Gff3Features
 *
 * LCJ update for phgv2:  send int ALT headers instead of graph?  These will have sample
 * name and regions for the haplotype nodes.
 *
 */
fun makeGffFromHvcf(hvcfFile: String, centerGffs: Map<String, TreeMap<Position,ArrayList<Gff3Feature>>>, outputFile:String): Set<Gff3Feature> {

    var time = System.nanoTime()
    val pseudoGenomeGff3Features :MutableSet<Gff3Feature> = mutableSetOf()

    val vcfFileReader = VCFFileReader(File(hvcfFile), false)
    val altHeaders = parseALTHeader(vcfFileReader.header)

    val hvcfRecords = vcfFileReader.iterator().asSequence().toList()

    vcfFileReader.close()

    var curChrom = "NIL"
    var prevChrom = curChrom
    var offset = 0;

    // There will be multiple reference ranges per chromosome.
    var count = 0
    var total = 0

    // traverse the hvcf records, creating the entries
    for (hvcfRecord in hvcfRecords) {
        val hapId = hvcfRecord.alternateAlleles.first().displayString.removeSurrounding("<", ">")
        val altData = altHeaders[hapId] ?: continue
        val sampleName = altData.sampleName()
        val refRange = altData.refRange
        curChrom = hvcfRecord.contig
        if (!curChrom.equals(prevChrom) && !prevChrom.equals("NIL")) {
            // When score is < 0, it will print Gff3Constants.UNDEFINED_FIELD_VALUEm, which is "."
            // Create "chromosome" Entry for gff
            val gff3FeatureEntry = createGffChromosomeEntry(prevChrom, offset)

            // This unfortunately comes at the end, not the beginning. of this chromosome list. But we don't
            // know the size until we've finished with haplotypes for this chrom.
            pseudoGenomeGff3Features.add(gff3FeatureEntry)
            offset = 0 // starting a new chromosome, offset becomes 0
        };
        prevChrom = curChrom

        // Get GFF entries for this sample
        val asmCenterGffEntries = centerGffs.get(sampleName)
        if (asmCenterGffEntries == null || asmCenterGffEntries.size == 0) {
            // myLogger.info("makeGffFromPath: WARNING: no gff entries for taxon ${taxon.name}")
            continue // go on to next one
        }
        // get regions that make up this hapid
        val regions = altData.regions

        // Create a region list, flipping any negative strand entries
        // This is needed for the getPseudoGFFCoordsMultipleRegions() function
        // where all the regions need to be on the forward strand.
        val regionList = mutableListOf<IntRange>()
        for (region in regions){
            if (region.first > region.second) {
                regionList.add(region.second.position .. region.first.position)
            } else {
                regionList.add(region.first.position .. region.second.position)
            }
        }


        // Could do this in 1`shot - just get all GFF entries that fall somewhere between
        // the start of the first region and the end of the second region.
        // However, if we do this, we could pick up an entry that falls completely in a gap between regions
        // and that GFF annotation is actually not represented here.
        // val regionGffEntries = getOverlappingEntriesFromGff(hvcfRecord.contig, regionList[0].first..regionList[regions.size-1].last, asmCenterGffEntries)

        // Add overlaps from all regions to a single Set,
        // this will prevent duplicates
        val regionOverlapLists = mutableSetOf<Gff3Feature>()
        for (region in regions) {
            val refRange = ReferenceRange(hvcfRecord.contig, region.first.position, region.second.position)

            // Grab the gff entries that overlap with the asm region used to create this haplotype
            // there may be multiple entries that overlap, and there may be multiple regions from
            // this Aseembly which make up the haplotype.
            val regionGffEntries = getOverlappingEntriesFromGff(hvcfRecord.contig, refRange.start..refRange.end, asmCenterGffEntries)
            regionOverlapLists.addAll(regionGffEntries)
        }

        // Sort the set of Gff3Features by contig, start, end
        val sortedGffEntries = regionOverlapLists.sortedWith(Comparator { first: Gff3Feature, second: Gff3Feature ->
            var contigCompare = first.contig.compareTo(second.contig)
            if (contigCompare == 0) {
                var startCompare = first.start.compareTo(second.start)
                if (startCompare == 0) {
                    first.end.compareTo(second.end)
                } else {
                    first.start.compareTo(second.start)
                }
            } else {
                first.contig.compareTo(second.contig)
            }
        })

        // Create the Gff3Feature entry, add to list
        for (entry in sortedGffEntries) {
            val newRange = getPseudoGFFCoordsMultipleRegions( entry.start..entry.end, regionList, offset)

            // Update the entry to be specific to gff coordinates.
            // Add an annotations for:
            //   reference_range_id
            //   haplotypes_id
            //   haplotype taxon
            //   asm coordinates for haplotype
            // The haplotype asm coordinates - can be used to distinguish partial
            // vs full mapping based on size of asm sequence and where the coordinates overlapped
            val annotations:MutableMap<String,List<String>> = mutableMapOf()
            annotations.putAll(entry.attributes)
            val refRangeList = listOf("${refRange}")
            annotations.put("referenceRangeID",refRangeList)
            val haplotypeIDList = listOf("${hapId}")
            annotations.put("haplotypeID",haplotypeIDList)
            val taxonList = listOf("${sampleName}") // the taxon from which these entries came.
            annotations.put("taxon",taxonList)
            // regions is a List<Pair<Position,Position>> - convert to List<String>, add to annotations
            val annoList = regions.map { "${it.first.contig}:${it.first.position}-${it.second.position}" }
            annotations.put("halotypeAsmCoordinates", annoList)

            // Create the Gff3Feature entry, add to list
            val gff3FeatureEntry = Gff3FeatureImpl(entry.contig, entry.source, entry.type, newRange.first, newRange.last, entry.score, entry.strand, entry.phase, annotations)
            pseudoGenomeGff3Features.add(gff3FeatureEntry)
        }
        // set offset to be offset plus length of all regions in the regionList
        offset += regionList.sumOf { it.last - it.first + 1 } // add 1 because are inclusive/inclusive, we need a correct count

    }

    // Create the last "chromosome" entry for this file
    val gff3FeatureEntry = createGffChromosomeEntry(curChrom, offset)
    pseudoGenomeGff3Features.add(gff3FeatureEntry)

    if (count > 0) total += count
    myLogger.info("makeGffFromPath: all entries finished, total: ${total}")
    if (outputFile != null && outputFile !="") {
        myLogger.info("makeGffFromPath: writing pseudo gff file to ${outputFile}")
        writeGffFile(outputFile, pseudoGenomeGff3Features,null,null)
    }
    val endingTime = (System.nanoTime() - time)/1e9
    myLogger.info("makeGffFromPath: time to process path: ${endingTime}")
    return pseudoGenomeGff3Features
}



/**
 * Using htsjdk classes,  writeGffFile will separately write comments, sequenceRegions, and features.
 * The list of comments (List<String>) and regions (Set<SequenceRegion>) are optional
 * The list of features (Set<Gff3Feature>) is required.
 *
 * The data is written to a GFF3 formatted file.
 */
fun writeGffFile(outputFile: String, features: Set<Gff3Feature>, comments: List<String>?, regions: Set<SequenceRegion>?) {

    try {

        Gff3Writer(Paths.get(outputFile)).use { writer ->
            // comments are optional
            comments?.forEach {
                writer.addComment(it)
            }

            // regions are optional
            regions?.forEach {
                writer.addDirective(Gff3Codec.Gff3Directive.SEQUENCE_REGION_DIRECTIVE, it)
            }
            // There should be features!
            features.forEach {
                writer.addFeature(it)
            }

        }
    } catch (ioe: IOException) {
        throw TribbleException("Error writing to file $outputFile", ioe)
    }
}

/**
 * Search the gff entry map for overlaps with the haplotypeNode asm coordinates.
 * Requires as input a map keyed by a Position object as well as the contig
 * and coordinates from the haplotypeNode.
 */
fun getOverlappingEntriesFromGff(contig: String, haplotypeRange:IntRange, asmCenterGffs: TreeMap<Position,ArrayList<Gff3Feature>>): Set<Gff3Feature> {

    val gff3Features = HashSet<Gff3Feature>() // no ordering here!
    var hapStartPos = Position(contig,haplotypeRange.first)
    var hapEndPos = Position(contig,haplotypeRange.endInclusive)

    // during assembly processing, we did not change the direction of the
    // assembly coordinates.  If the start > end, it means it was aligned
    // on the reverse strand.  Flip them here so we don't have problems
    // pulling the submap below.
    if (hapStartPos > hapEndPos) {
        var temp = hapStartPos
        hapStartPos = hapEndPos
        hapEndPos = temp
    }
    // The asmCenterGffs has a key that is a Position object whose "position" is created from the
    // middle point of the coordinate range.  We grab a submap of the gff that contains this endpoint.
    // the submap will be based on the haplotype node's asm coordinates.  "floorKey" and "ceilingKey"
    // are used as the specific haplotype node start/end positions may not be keys in the map.

    // Start with floorKey 1 less than the haplotype node start position to get more
    // gff3 entries around the center position.
    var floorKeyTemp = asmCenterGffs.floorKey(hapStartPos)
    var floorKey = floorKeyTemp
    if (floorKeyTemp != null) {
        floorKey = Position(floorKeyTemp.contig, floorKeyTemp.position-1)
    }

    // Because submap below grabs the values in an inclusive/exclusive manner,
    // the ceilingKey value must be 1 greater than what we want.
    var ceilingKeyTemp = asmCenterGffs.ceilingKey(hapEndPos)
    var ceilingKey = ceilingKeyTemp
    if (ceilingKeyTemp != null) {
        ceilingKey = Position(ceilingKeyTemp.contig,ceilingKeyTemp.position+1)
    }

    if (floorKey == null) {
        if (ceilingKey != null) {
            floorKey = Position(ceilingKey.contig,ceilingKey.position-1)
        }
    }
    if (ceilingKey == null) {
        if (floorKey != null) {
            ceilingKey = Position(floorKey.contig,floorKey.position+1)
        }
    }

    val subMap1 = if (floorKey != null ) asmCenterGffs.subMap(floorKey,ceilingKey).flatMap { it -> it.value } else null

    // There are many instances where the coordinates from floorKey or ceilingKey do not overlap the
    // haplotypeNode's asm coordinates.  So the entries returned in the submap must be checked for overlap.
    if (subMap1 != null) {
        for (gff3entry in subMap1) {
            if (gff3entry.contig.equals(contig)) { // do I need this conditional ??
                if ((haplotypeRange.first <= gff3entry.end) && (haplotypeRange.endInclusive >= gff3entry.start)) {
                    gff3Features.add(gff3entry)
                }
            }
        }
    }
    return gff3Features
}

/**
 * This function determines which regions in a set overlap
 * a specific GFF entry.  The regions are a list of IntRanges
 * that comprise a haplotype.  Only some of them may overlap
 * the GFF entry.
 *
 * It returns a list of regions that overlap, the gap size between the regions,
 * and an offset from the start of the haplotype regions if the initial regions
 * in the list do not over lap the GFF entry.
 *
 * If there is only 1 region that overlaps, the gapSize will be 0.
 *
 */
fun findRegionsOverlappingGFF(asmGffRange: IntRange, regions:List<IntRange>):Triple<List<IntRange>,Int,Int> {
    // this function is used to determine which of the regions
    // in the haplotype node overlap with the GFF entries.
    var gapSize=0
    val overlappingRegions = mutableListOf<IntRange>()
    var offSet = 0
    // check each region to see if it overlaps the GFF entry
    // Until you reach the first region that overlaps, count both the
    // number of bases in the region and the number in the gap before
    // the next region.  Add these values to the "offSet" value.
    for (idx in 0 until  regions.size ) {
        if ((asmGffRange.start <= regions[idx].endInclusive) && (asmGffRange.endInclusive >= regions[idx].start)) {
            overlappingRegions.add(regions[idx])
        } else {
            if (overlappingRegions.size == 0) {
                // add the size of the non-overlapping region to the initial offset
                offSet += regions[idx].endInclusive - regions[idx].start +1 // plus 1 as is inclusive/inclusive
                // We don't need to calculate gap size as they are not part of the sequence.
            }
        }
    }
    // if there is more than 1 overlapping region, calculate the gap size
    if (overlappingRegions.size > 1) {
        for (idx in 0 until overlappingRegions.size-1) {
            gapSize += overlappingRegions[idx+1].start - overlappingRegions[idx].endInclusive -1
        }
    }
    return Triple(overlappingRegions,gapSize,offSet)
}

/**
 * This function takes a GFF3 entry,  list of Assembly regions and a pseudo assembly
 * fasta offset to determine coordinates for a pseudo assembly GFF entry.
 *
 * The colling procedure must ensure all region entries are converted to the forward strand
 * coordinates.
 *
 * Not all regions in the list may overlap the GFF entry.  The function handles
 * overlapping with the initial regions, just regions in the middle, or just the end region.
 * If more than 1 region in the list overlaps the GFF entry coordinates, those regiones
 * should be consecutive in the list.
 */
fun getPseudoGFFCoordsMultipleRegions(asmGffRange: IntRange, regions:List<IntRange>, baseOffset: Int):IntRange {
    // Because the range is inclusive/inclusive, the size is actually + 1
    // But the value we need to add is just the difference between the start and end.
    val featureSize = asmGffRange.endInclusive - asmGffRange.start

    val (overlappingRegions,gapSize,regionStartOffset) = findRegionsOverlappingGFF(asmGffRange,regions)
    // The regionStartOffset returned above indicates number of basepairs in the regions
    // on the list that do not overlap the GFF entry.  This is only regions prior to the
    // first region that overlaps the GFF entry.  This value is added to the offset indicating
    // haplotype sequence start position in the pseudo-genome.
    val offset = baseOffset + regionStartOffset

    // startDiff is negative if the haplotype start coordinate begins after
    // the assembly gff entry start value.
    var hapStart = overlappingRegions[0].start
    var hapEnd = overlappingRegions[overlappingRegions.size-1].endInclusive

    val startDiff = asmGffRange.start - hapStart
    val endDiff = hapEnd - asmGffRange.endInclusive
    // endDiff and endAdjust are negative if the haplotype node coordinate end value
    // is less than the assembly gff end value, meaning it doesn't cover the full GFF entry
    val endAdjust = if (endDiff > 0) 0 else endDiff

    // pgStart is > 1 if the haplotypeNode starts  prior to the gff entry
    val pgStart = if (startDiff > 0) startDiff else 1

    // If startDiff is 0, or pgStart > 1, the haplotype node includes the beginning
    // of the assembly gff entry.  The end is then the size of this feature, minus
    // any part of the end that is not included in the haplotype node.
    // Otherwise, if pgStart == 1, it means the beginning of the feature was not included.
    // The end is then adjusted by subtracting a count of the missing bps from both
    // ends from the feature size.

    // GapSize is subtracted from the end coordinates to account for the missing bps
    // between regions that comprise the haplotype node.
    val pgEnd = if (startDiff == 0 || pgStart > 1) pgStart + featureSize + endAdjust - gapSize
        else pgStart + featureSize + endAdjust + startDiff - gapSize
    val
            returnStart = pgStart+offset
    val returnEnd = pgEnd+offset

    return pgStart+offset..pgEnd+offset
}

/**
 * Create the "chromosome" type line in the gff file
 */
fun createGffChromosomeEntry(prevChrom:String, offset:Int):Gff3Feature {
    var annotations:MutableMap<String,List<String>> = mutableMapOf()
    annotations.put("ID",listOf(prevChrom))
    annotations.put("Name",listOf("chromosome:PHG-pseudo_genome"))
    var chromGffSource = "PHG"
    // When score or phase values are < 0, Gff3Feature writer will print Gff3Constants.UNDEFINED_FIELD_VALUE,
    // which is "." and what we want here
    val minus1 = -1
    val minus1Double = minus1.toDouble()
    val gff3FeatureEntry = Gff3FeatureImpl(prevChrom,chromGffSource,"chromosome",1,
            offset,minus1Double, Strand.NONE, -1, annotations)

    return gff3FeatureEntry
}

/**
 * Counts the number of times a single GFF ID appears.
 * User may specify the full ID or just part of it.  The
 * code checks for GFF entries where the ID contains the
 * user specified string.
 *
 * FOr example, if the user wanted to check for and ID of "Zm00001e000002":
 * this could show up in the GFF3 file attributes column as:
 *  ID=gene:Zm00001e000002
 *  ID=transcript:Zm00001e000002_T001
 *  ID=Zm00001e000002
 *
 *  Or, there may be no ID, but a Parent attribute, e.g.:
 *  Parent=transcript:Zm00001e000002_T001
 *
 * Base on the above, the code checks each feature first for an ID
 * field in the attributes column, and if not found, for a Parent fiels
 * in the attributes column.  Then verifies is the user "idToMatch" string
 * is contained in the ID or Parent value.
 *
 * THe original code assumed there would always be an ID, but the ID is
 * only required for attributes that have children, it is optional otherwise.
 *
 */
fun gffSingleIDcount(idToMatch:String, gffSet:Set<Gff3Feature>):Int {

    var idCount = 0
    for (gff3Feature in gffSet) {
        val featureId = gff3Feature.getBaseData().id ?: ""
        if (featureId.contains(idToMatch)) {
            idCount++
        } else {
            val parentId = gff3Feature.attributes["Parent"] ?: null
            if (parentId != null && parentId[0].contains(idToMatch)) {
                idCount++
            }
        }
    }
    return idCount
}

/**
 * Takes a set of Gff3Features, creates a mapping of
 * GFF id (ID from the attributes column) to number of times it appears
 * The "count" parameter tells by what amount to shorten the ID.
 *
 * For example:  If the user only wants the first 8 characters of the
 * ID field considered, count would be 8.   This would allow for counting
 * something like all instances of RST00015*.
 *
 * If "count" is > than the length of the ID, the full value is used as the key.
 *
 * returns:  a mapping of ID (subsetted to first "count" letters) to number of
 * times it occurs.
 */
fun gffAllIDcount(gffSet:Set<Gff3Feature>, count:Int):Map<String,Int> {
    // filter for  null ID before mapping
    val idToCountMap = gffSet.asSequence().filter {it.id != null}.groupingBy { it.id.take(count)}.eachCount()
    return idToCountMap
}

/**
 * Counts the gff entries per chromosome
 *
 * returns:  Map of chromsome -> number of entries
 */
fun getGFFEntriesPerChrom(gffSet:Set<Gff3Feature>):Map<String,Int> {
    val chromCounts = gffSet.groupingBy{it.contig}.eachCount()
    return chromCounts
}

/**
 * Counts the number of distinct ID values either across the full
 * GFF3Feature set (if contig==null) or for a specific contig
 */
fun countDistinctID(gffSet:Set<Gff3Feature>,contig:String?): Int {
    var distinctCount = 0
    if (contig == null) {
        // count all distinct IDs
       distinctCount= gffSet.asSequence()
           .filter{it.id != null}
           .distinctBy {gff3Feature -> gff3Feature.id }
           .count()
    } else {
        distinctCount = gffSet.asSequence()
                .filter{it.contig.equals(contig)}
                .filter{it.id != null}
                .distinctBy {it.id}
                .count()
    }
    return distinctCount
}

/**
 * Get list of distinct GFF IDs in the full Gff3Feature Set
 *  NOTE: what this returns is based on how the GFF3 has the IDs defined.
 *  If user wants a count  of something like "Zm000a2" but the ID is:
 *    "ID=gene:Zm000a2", it won't be what they want.  Because "ID="gene:Zm000a2"
 *    and "ID=transcript:Zm000a2" are different IDs.
 *
 * returns: List<String> where String is the ID value
 */
fun getDistinctGffIds(gffSet:Set<Gff3Feature>): List<String> {
    // Checking for null as ID is not required
    val distinctIDs = gffSet.asSequence()
            .filter{it.id !=null}
            .map{it -> it.id}
            .distinct().toList()
    return distinctIDs
}

/**
 * Creates map of contig(chrom) to set of distinct Gff3Feature ids
 * associated with that contig
 */
fun getDistinctGffIdsByChrom(gffSet:Set<Gff3Feature>): Map<String,Set<String>> {
    val chromToIds = TreeMap<String,MutableSet<String>>()

    for (feature in gffSet) {
        val contig = feature.contig
        val idSet =  chromToIds?.get(contig)?:mutableSetOf<String>()
        if (feature.id != null) {
            idSet.add(feature.id)
        }
        chromToIds.put(contig,idSet)
    }
    return chromToIds
}

/**
 * Calculates the number of BPs for each contig/chromosome in the GFF file
 * Overlapping positions are only counted once.
 *
 * return:  A map of chromosome to number of included positions
 */
fun sumPerChromGFFBasePairs(gffSet:Set<Gff3Feature>): Map<String,Int> {

    val chromSizeMap = mutableMapOf<String,Int>() // holds total length of each chromsome

    // Using RangeSets as GFFs may have overlapping positions among different types (Exon, CDS, etc)
    val chromRangeSet: RangeSet<Position> = TreeRangeSet.create()
    gffSet.asSequence().forEach {
        //val chrom = Chromosome.instance(it.contig) // because db used Chromosome object, must get gff to match
        if (it.type.equals("chromosome")) {
            chromSizeMap.put(it.contig,it.end)
        } else {
            val range= Range.closed(Position(it.contig,it.start),Position(it.contig,it.end))
            chromRangeSet.add(range)
        }

    }

    // process the overlaps for each chrom
    // The number of base pairs represented in the GFF is the total number of
    // non-overlapping basepairs that occur for each chromosome.  The lower bound
    // is 1 and the upper bound is the size of the chromosome.
    var chromBPsumMap = chromSizeMap.map {(chr,size) ->
        val testRange = Range.closed(Position(chr,1),Position(chr,size))
        val includedBP = chromRangeSet.subRangeSet(testRange).asRanges()
                .map { range ->
                    range.upperEndpoint().position - range.lowerEndpoint().position +1
                }.sum()
        Pair(chr, includedBP)
    }.toMap()


    return chromBPsumMap
}

fun percentPerChromGFFBasePairs(gffSet:Set<Gff3Feature>): Map<String,Double> {
    val chromSizeMap = mutableMapOf<String,Int>() // holds total length of each chromsome

    // Using RangeSets as GFFs may have overlapping positions among different types (Exon, CDS, etc)
    val chromRangeSet: RangeSet<Position> = TreeRangeSet.create()
    gffSet.asSequence().forEach {
        //val chrom = Chromosome.instance(it.contig) // because db used Chromosome object, must get gff to match
        if (it.type.equals("chromosome")) {
            chromSizeMap.put(it.contig,it.end)
        } else {
            val range= Range.closed(Position(it.contig,it.start),Position(it.contig,it.end))
            chromRangeSet.add(range)
        }

    }

    // process the overlaps for each chrom
    // The number of base pairs represented in the GFF is the total number of
    // non-overlapping basepairs that occur for each chromosome.  The lower bound
    // is 1 and the upper bound is the size of the chromosome.
    var chromBPsumMap = chromSizeMap.map {(chr,size) ->
        //val chrom = Chromosome.instance(chr)
        val testRange = Range.closed(Position(chr,1),Position(chr,size))
        val includedBP = chromRangeSet.subRangeSet(testRange).asRanges()
            .map { range ->
                range.upperEndpoint().position - range.lowerEndpoint().position +1
            }.sum()
        // When multiplying by 100 (because of double times int?), 29 of 100 becomes percent 28.99999 instead of 29
        val percent = (includedBP.toDouble()/size)
        Pair(chr, percent)
    }.toMap()


    return chromBPsumMap
}
/**
 * Calculates the number of BPs for each contig/chromosome that were not
 * present in the GFF file.  For example:  If the chromosome is of length
 * 2000, and there are only 850 total base pairs for that chromosome
 * in the GFF, then the non-represented number is 1150
 *
 * return: a map of chromosome to number of bps NOT represented in the GFF file
 */
fun sumPerChromNonGFFBasePairs(gffSet:Set<Gff3Feature>): Map<String,Int> {

    val chromSizeMap = mutableMapOf<String,Int>() // holds total length of each chromsome

    // Using RangeSets as GFFs may have overlapping positions among different types (Exon, CDS, etc)
    val chromRangeSet: RangeSet<Position> = TreeRangeSet.create()
    gffSet.asSequence().forEach {
       // val chrom = Chromosome.instance(it.contig) // because db used Chromosome object, must get gff to match
        if (it.type.equals("chromosome")) {
            chromSizeMap.put(it.contig,it.end)
        } else {
            val range= Range.closed(Position(it.contig,it.start),Position(it.contig,it.end))
            chromRangeSet.add(range)
        }
    }

    // process the overlaps for each chrom
    // The number of base pairs represented in the GFF is the total number of
    // non-overlapping basepairs that occur for each chromosome.  The upper and
    // lower bound are 1 and the size of the chromosome.
    // Subtract the overlapping bps from the total bps in the chromosome to get
    // the number of chromsome base-pairs that are not represented in this GFF set.
    var chromBPsumMap = chromSizeMap.map { (chr,size) ->
        //val chrom = Chromosome.instance(chr)
        val testRange = Range.closed(Position(chr,1),Position(chr,size))
        val includedBP = chromRangeSet.subRangeSet(testRange).asRanges()
                .map { range ->
                    range.upperEndpoint().position - range.lowerEndpoint().position +1
                }.sum()

        // This line differentiates this method from sumPerChromGFFBasePairs()
        Pair(chr,size-includedBP)
    }.toMap()

    return chromBPsumMap
}

fun percentPerChromNonGFFBasePairs(gffSet:Set<Gff3Feature>): Map<String,Double> {

    val chromSizeMap = mutableMapOf<String,Int>() // holds total length of each chromsome

    // Using RangeSets as GFFs may have overlapping positions among different types (Exon, CDS, etc)
    val chromRangeSet: RangeSet<Position> = TreeRangeSet.create()
    gffSet.asSequence().forEach {
        //val chrom = Chromosome.instance(it.contig) // because db used Chromosome object, must get gff to match
        if (it.type.equals("chromosome")) {
            chromSizeMap.put(it.contig,it.end)
        } else {
            val range= Range.closed(Position(it.contig,it.start),Position(it.contig,it.end))
            chromRangeSet.add(range)
        }
    }

    // process the overlaps for each chrom
    // The number of base pairs represented in the GFF is the total number of
    // non-overlapping basepairs that occur for each chromosome.  The upper and
    // lower bound are 1 and the size of the chromosome.
    // Subtract the overlapping bps from the total bps in the chromosome to get
    // the number of chromsome base-pairs that are not represented in this GFF set.
    var chromBPsumMap = chromSizeMap.map { (chr,size) ->
        //val chrom = Chromosome.instance(chr)
        val testRange = Range.closed(Position(chr,1),Position(chr,size))
        val includedBP = chromRangeSet.subRangeSet(testRange).asRanges()
            .map { range ->
                range.upperEndpoint().position - range.lowerEndpoint().position +1
            }.sum()

        // When multiplying by 100, the values get off slightly.  29 out of 100 becomes 28.99999999 instead of 29
        // and 87 of 200 becomes 56.499999999 instead of 56.4
        val percent = ((size-includedBP).toDouble()/size)
        // This line differentiates this method from sumPerChromGFFBasePairs()
        Pair(chr,percent)
    }.toMap()

    return chromBPsumMap
}

// This function is used when we call the PHG plugin.  That test is
// temporarily removed while we determine how to do integration testing
// of PHGdbAccess
//fun createTaxaListFromFileOrString(taxa: String): TaxaList {
//    if (taxa.endsWith(".txt")) {
//        try {
//            val builder = TaxaListBuilder()
//            Utils.getBufferedReader(taxa).use { br ->
//                var line = br.readLine()
//                val sep = Pattern.compile("\\s+")
//                while (line != null) {
//                    line = line.trim()
//                    val parsedline = sep.split(line)
//                    for (idx in parsedline.indices) {
//                        if (parsedline[idx] != null || parsedline[idx]!!.length != 0) {
//                            builder.add(parsedline[idx])
//                        }
//                    }
//                    line = br.readLine()
//                }
//            }
//            return builder.build()
//        } catch (exc: Exception) {
//            throw IllegalArgumentException("Error reading taxa file ${taxa}")
//        }
//    } else {
//        val taxa = taxa.trim().split(",")
//        return TaxaListBuilder().addAll(taxa).build()
//    }
//}