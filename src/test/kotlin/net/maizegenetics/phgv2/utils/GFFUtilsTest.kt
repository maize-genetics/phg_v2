package net.maizegenetics.phgv2.utils

import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.testing.test
import com.google.common.collect.Range
import com.google.common.collect.RangeSet
import com.google.common.collect.TreeRangeSet
import htsjdk.tribble.annotation.Strand
import htsjdk.tribble.gff.Gff3Feature
import htsjdk.tribble.gff.SequenceRegion
import junit.framework.TestCase
import net.maizegenetics.phgv2.cli.PathsToGff
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.Assertions.assertFalse
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import java.io.IOException
import java.nio.file.Files
import java.nio.file.Paths
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertFails

private const val i = 134374722

@ExtendWith(TestExtension::class)
class GFFUtilsTest {

    companion object {
        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"
        val outputDir = "${TestExtension.tempDir}/GFFTests/"
        val dataDir = outputDir + "data/"
        val gffFileWithHeadersB73 = "${dataDir}/smallFullGff_withHeaders_withTE_B73.txt"
        val gffFileWithHeadersCML103 = "${dataDir}/smallFullGff_withHeaders_CML103.txt"
        val gffFileShortBPTesting = "${dataDir}/shortGFF_bpCovTesting.txt"
        val keyFile = "${dataDir}/keyFileJunit.txt"
        val testOutputFile = "${dataDir}/gffOutputTestFile.gff3"

        @JvmStatic
        @BeforeAll
        fun setup() {
            // Delete recursively the tempDir
            File(TestExtension.tempDir).deleteRecursively()
            File(dataDir).mkdirs()
            try {

                Files.createDirectories(Paths.get(dataDir))
            } catch (ioe: IOException) {
                throw IllegalStateException("GFFTests:setup: error deleting/creating folders: " + ioe.message)
            }
            createGFFFileWithHeadersB73(gffFileWithHeadersB73)
            createGFFFileWithHeadersCML103(gffFileWithHeadersCML103)
            createGFFFileForBPTests(gffFileShortBPTesting)
            createGFFKeyfile(keyFile,gffFileWithHeadersB73,gffFileWithHeadersCML103)
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(TestExtension.tempDir).deleteRecursively()
        }
    }


    @Test
    fun testPathsToGff() {
        // This is testing the actual clikt class.  Everything below here are tests of
        // the individual functions that are called by the clikt class.
        // THis is the keyfile we had fortesting with Michelle's data in phgv1:
        //    /Users/lcj34/notes_files/phg_2018/new_features/phg493_GFF_plugin_fromAsmCoords/testing/keyFile.txt
        // TODO: LCJ - remove hard coding - store files somewhere
        // TODO : this test still fails, and it doesn't pick up the early genes - debug that !
        // It isn't finindht them in getOverlappingEntriesFromGff()
        val keyFile = "/Users/lcj34/notes_files/phg_v2/newFeatures/pathsToGFF/testKeyFile.txt"
        val hvcfFile = "/Users/lcj34/notes_files/phg_v2/newFeatures/pathsToGFF/Imputation.h.vcf"
        val outputFile = "/Users/lcj34/notes_files/phg_v2/newFeatures/pathsToGff/junit_tests/testPathsToGffOutput.gff3"

        val pathsToGff = PathsToGff()
        //val goodParamsTest = pathsToGff.test("--key-file ${keyFile} --hvcf-file ${hvcfFile} --output-file ${testOutputFile}")
        val goodParamsTest = pathsToGff.test("--key-file ${keyFile} --hvcf-file ${hvcfFile} --output-file ${outputFile}")

        assertEquals(0, goodParamsTest.statusCode)
        // Verify specific entries
        // read the outputFile using bufferedReader
        val newGff = bufferedReader(outputFile).readLines()

        // The imputed haplotype had positions less that 27000 belonging to LineA, and greater than that
        // belonging to LineA.  So the gff should be picking up entries from LineA that have positions
        // less than 2700.  Entries with starting positions > 2700 should be entries from LineB GFF
        // THese entries are distinguished by LineA or LineB in their attribute annotations
        // Because there may be offset differences in the composite genome, we cannot check based on the position alone
        // so we're checking based on the attribute string, which will contain either a LineA or LineB identifier

        // verify there is a line in the file that contains "ID=TE:LineB_TE22;Name=Copia_LineB" in the line
        val lineWithTE = newGff.filter { it.contains("ID=TE:LineB_TE22;Name=Copia_LineB") }
        assertTrue(lineWithTE.size > 0)
        //verify there is no line in the file that contains "ID=TE:LineA_TE22;Name=Copia_LineA"
        val lineWithTEA = newGff.filter { it.contains("ID=TE:LineA_TE22;Name=Copia_LineA") }
        assertTrue(lineWithTEA.size == 0)

        // Check beginning entries are lineA specific not lineB
        // this gff's mikado gene is in the 300 range, so should show LIneA_mikado_gene,
        // and not LineB_mikado_gene
        val lineAEntries = newGff.filter { it.contains("LineA_mikado_gene") }
        assertTrue(lineAEntries.size > 0)
        val lineBEntries = newGff.filter { it.contains("LineB_mikado_gene") }
        assertTrue(lineBEntries.size == 0)

        // TODO - finish this

        // Need to add some assertions.  How many entries do I expect? etc
        // Verify specific entries are in the output file.
        println("Finished testPathsToGff")


    }

    @Test
    fun loadGffsToGFFFeatureTest() {
        // this returns a Map<String, TreeMap<Position,ArrayList<Gff3Feature>>> where
        // the key is the taxon, and the value is a TreeMap of mid-center position to featues that overlap that position

        // NOte more detailed verificaion of the treeMap itself are performed in below
        // in test testCreateTreeMapFromFeatureCenter
        val mapTreeMap = loadGFFsToGff3Feature(keyFile)
        assertTrue(mapTreeMap.containsKey("B73"))
        assertTrue(mapTreeMap.containsKey("CML103"))
        assertFalse(mapTreeMap.containsKey("KY59"))

        // B73 GFF file has 56 non-comment entries, CML103 has 35
        assertEquals(56, mapTreeMap["B73"]!!.size)
        assertEquals(35, mapTreeMap["CML103"]!!.size)
    }

    @Test
    fun testWriteGFFFile() {
        val gffFeatures = readGFFtoGff3Feature(gffFileShortBPTesting)

        val comments = listOf("##gff-version   3","# genebuild 2019-cshl")

        // regions: Set<SequenceRegion>
        val regions = HashSet<SequenceRegion>()
        // SequenceRegion will not accept the constructor call, though the class is Public
        //regions.add(SequenceRegion("chr1",1,789789))
        writeGffFile(testOutputFile, gffFeatures.toSet(), comments, regions)

        val gffFeaturesReturned = readGFFtoGff3Feature(testOutputFile)
        assertEquals(gffFeatures.size, gffFeaturesReturned.size)
        val aFeature = gffFeatures.first()
        TestCase.assertTrue(gffFeaturesReturned.contains(aFeature))
        val lastFeature = gffFeatures.last()
        TestCase.assertTrue(gffFeaturesReturned.contains(lastFeature))
    }

    @Test
    fun testCreateGFFChromosomeEntry() {
        val chrom = "myChrom"
        val offset = 303303
        val chromEntry = createGffChromosomeEntry(chrom,offset)

        assertEquals(chromEntry.contig,chrom)
        assertEquals(chrom,chromEntry.id)
        assertEquals(1,chromEntry.start)
        assertEquals(303303,chromEntry.end)
        assertEquals(Strand.NONE,chromEntry.strand)
    }

    @Test
    fun testGffAllIDcount() {
        val gffFeatures = readGFFtoGff3Feature(gffFileShortBPTesting).toSet()
        val idsWithCount = gffAllIDcount(gffFeatures,4)
        assertEquals(4,idsWithCount.keys.size)
        assertTrue(idsWithCount.containsKey("gene"))
        assertTrue(idsWithCount.containsKey("tran"))
        assertEquals(4, idsWithCount["gene"])
        assertEquals(2, idsWithCount["tran"])
        assertTrue(idsWithCount.containsKey("1"))
        assertTrue(idsWithCount.containsKey("2"))
        assertEquals(1, idsWithCount["1"])
        assertEquals(1, idsWithCount["2"])

    }

    @Test
    fun testCountDistinctIds() {
        val gffFeatures = readGFFtoGff3Feature(gffFileShortBPTesting).toSet()
        println("Size of gffFeatures: ${gffFeatures.size}")
        val distinctIds = countDistinctID(gffFeatures,null)
        assertEquals(8, distinctIds)

        // chr1 and chr2 in this short set have the same number of distinct ids
        val chr1DistinctIds = countDistinctID(gffFeatures,"chr1")
        assertEquals(4, chr1DistinctIds)
        val chr2DistinctIds = countDistinctID(gffFeatures,"chr2")
        assertEquals(4,chr2DistinctIds)
    }

    @Test
    fun testGetTaxonToGffFileMap() {
        val taxaFileMap = getTaxonToGffFileMap(keyFile)
        assertEquals(2,taxaFileMap.size)
        TestCase.assertTrue(taxaFileMap.containsKey("B73"))
        TestCase.assertTrue(taxaFileMap.containsKey("CML103"))
        TestCase.assertTrue(taxaFileMap["B73"]!!.contains("smallFullGff_withHeaders_withTE_B73.txt"))
        TestCase.assertTrue(taxaFileMap["CML103"]!!.contains("smallFullGff_withHeaders_CML103.txt"))
    }

    @Test
    fun readGFFhtsjdk() {

        // THis skips over the headers.
        // This test exercises the code with a small GFF that is created programmatically.
        // Additional non-junit tests have been run with a full sized GFF file (B73 v5 version)
        val gffFeatureSet = readGFFtoGff3Feature(gffFileWithHeadersB73)

        //test what the attributes look like
        assertEquals(96,gffFeatureSet.size)

        // verify the first attribute of the first feature has ID 1 and the second attribute is Name and the mapped value contains "chromosome"
        val firstFeature = gffFeatureSet.first()
        val firstAttributes = firstFeature.attributes
        val firstAttribute = firstAttributes.entries.first()
        assertEquals("ID", firstAttribute.key)
        assertEquals("1", firstAttribute.value[0])
        val secondAttribute = firstAttributes.entries.elementAt(1)
        assertEquals("Name", secondAttribute.key)
        assertEquals("chromosome:Zm-B73-REFERENCE-NAM-5.0:1:1:308452471:1", secondAttribute.value[0])

        // verify the last feature has chromsome chr2 with start 1005303 and end 1006303,
        // and the first attribute is ID and the value is TE:Lynn_TE1
        val lastFeature = gffFeatureSet.last()
        assertEquals("chr2", lastFeature.contig)
        assertEquals(1005303, lastFeature.start)
        assertEquals(1006303, lastFeature.end)
        val lastAttributes = lastFeature.attributes
        val lastAttribute = lastAttributes.entries.first()
        assertEquals("ID", lastAttribute.key)
        assertEquals("TE:Lynn_TE22", lastAttribute.value[0])

    }

    @Test
    fun testCreateTreeMapFromFeatureCenter() {
        // NOTE - the readGFFtoGff3Feature() was tested in the case above (readGFFhtsjdk() )
        // we are separating that test from this so the treemap can be tested separately
        // The test gff created here has 96 data entries
        val features = readGFFtoGff3Feature(gffFileWithHeadersB73)
        assertEquals(96,features.size)
        // This returns TreeMap<Position,ArrayList<Gff3Feature>>
        val featureTreeMap = createTreeMapFromFeaturesCenter(features)

        // featureTreeMap will have fewer entries than features, as some of the features are
        // combined into a list and stored against a mid-point position
        TestCase.assertTrue(features.size > featureTreeMap.entries.size)

        val pos43988 = Position("chr1",43988)
        TestCase.assertTrue(featureTreeMap.keys.contains(pos43988))

        // There are 2 GFF file entries with start/stop of 41214/46762 - mid point = 43988 (on chromosome 1)
        // They are for ID=gene:Zm00001e000002 and ID=transcript:Zm00001e000002_T002
        TestCase.assertTrue(featureTreeMap.get(pos43988)?.size!! == 2)
        val entries43988 = featureTreeMap.get(pos43988)

        // The same 2 entries will be on the list, but their order isn't guaranteed
        // so verify the ID is one of those from the list.
        val id43988List = listOf("transcript:Zm00001e000002_T002","gene:Zm00001e000002")
        TestCase.assertTrue(id43988List.contains(entries43988!!.get(0).baseData.id))
        TestCase.assertTrue(id43988List.contains(entries43988!!.get(1).baseData.id))

        // There is 1 GFF File entry with start/stop of 1005003/1006003 - mid point = 1005803 (on chromosome 2)
        // Verify this mid-point position appears in the featureTreeMap with a single value against it
        val pos1005803 = Position("chr2", 1005803)
        TestCase.assertTrue(featureTreeMap.keys.contains(pos1005803))

        val entries1005803 = featureTreeMap.get(pos1005803)
        TestCase.assertTrue(entries1005803?.size!! == 1)
        TestCase.assertTrue(entries1005803!!.get(0).baseData.id == "TE:Lynn_TE22")


        // THere are 3 entries in the gff file for chrom 2 with start/end of 29923/31419, midpoint=30671
        val pos30671 = Position("chr2", 30671)
        val entries30671 = featureTreeMap.get(pos30671)

        // Not all entries have an ID (gff only requires ID for features that have children,
        // e.g. gene and mRNAs) and for those that span multiple lines.  The "exon" entry here has
        // a Parent attribute but no ID, so we'll verify the type of the 3 entries instead of ID
        val type30671List = listOf("exon","gene","mRNA")
        TestCase.assertTrue(entries30671!!.size == 3)
        TestCase.assertTrue(type30671List.contains(entries30671!!.get(0).baseData.type))
        TestCase.assertTrue(type30671List.contains(entries30671!!.get(1).baseData.type))
        TestCase.assertTrue(type30671List.contains(entries30671!!.get(2).baseData.type))

    }

    @Test
    fun testGff3FeatureOverlaps() {
        val time = System.nanoTime()
        val taxonToGFFfileMap = getTaxonToGffFileMap(keyFile)

        val centerGffs:MutableMap<String, TreeMap<Position, ArrayList<Gff3Feature>>> = mutableMapOf()

        for ((taxon, gffFile) in taxonToGFFfileMap) {
            val time = System.nanoTime()
            val features = readGFFtoGff3Feature(gffFile)
            val featureTreeMapCenter = createTreeMapFromFeaturesCenter(features)
            centerGffs.put(taxon, featureTreeMapCenter)
            val endTime = (System.nanoTime() - time)/1e9
            println("createTreeMapFromFeaturesCenter: time to load ${taxon} gff file: ${endTime}")
        }

        val name = "B73"
        // THe center of this range will not be a key in the map.
        // This verifies the code picks the next highest key and
        // that appropriate entries are identified.
        val testRange = 43067..43100
        val asmCenterGffEntries = centerGffs.get(name)

        if (asmCenterGffEntries == null || asmCenterGffEntries.size == 0) {
            assertFails{"No values in the asmCenterGffEntries map !!!"}
        }

        val pseudoGffEntries = getOverlappingEntriesFromGff("chr1", testRange, asmCenterGffEntries!!)
        for (entry in pseudoGffEntries) {
            println("${entry.contig} ${entry.start} ${entry.end} ${entry.type}")
        }

        // Result should be 4 based on the keyfile and gffFile created in setup()
        assertEquals(4, pseudoGffEntries.size)

    }

    @Test
    fun testGetPseudoGFFCoordsSingleRegion() {
        // ALl of the test cases here have only a single region for the haplotype.
        // Verify getPseudoGFFCoordsMultipleRegions() works with a single region

        // haplotype asm starts before the gff entry, asm ends before gff entry, 0 offset
        var gffCoords = 134374229..134374722
        var hapAsmCoords = 134370278..134374620
        var regions = mutableListOf(hapAsmCoords)
        var offset = 0

        var expectedResult = 3951..4342
        var pseudoGenomeCoords = getPseudoGFFCoordsMultipleRegions(gffCoords, regions, offset)
        assertEquals(expectedResult, pseudoGenomeCoords)

        // same as above, but offset from start of chrom is 1000
        offset = 1000
        expectedResult = 4951..5342
        pseudoGenomeCoords = getPseudoGFFCoordsMultipleRegions(gffCoords, regions, offset)
        assertEquals(expectedResult, pseudoGenomeCoords)

        // haplotype asm start is after the gff entry, ends after the gff entry (beginning
        // is chopped off, but includes all the end).  no offset
        offset = 0
        gffCoords = 155822773..155823423
        hapAsmCoords = 155823019..155827364
        regions = mutableListOf(hapAsmCoords)
        expectedResult = 1..405
        pseudoGenomeCoords = getPseudoGFFCoordsMultipleRegions(gffCoords, regions, offset)
        assertEquals(expectedResult, pseudoGenomeCoords)


        // full gff entry is embedded in the haplotype node sequence (haplotype start is
        // less than gff start, and haplotype end is greated than gff end)
        gffCoords = 142564613..142564838
        hapAsmCoords = 142564609..142570362
        regions = mutableListOf(hapAsmCoords)
        expectedResult = 4..229
        pseudoGenomeCoords = getPseudoGFFCoordsMultipleRegions(gffCoords, regions, offset)
        assertEquals(expectedResult, pseudoGenomeCoords)

        // asm gff and haplotype asm coordinates are the same
        gffCoords = 183334395..183338622
        hapAsmCoords = 183334395..183338622
        regions = mutableListOf(hapAsmCoords)
        expectedResult = 1..4228
        pseudoGenomeCoords = getPseudoGFFCoordsMultipleRegions(gffCoords, regions, offset)
        assertEquals(expectedResult, pseudoGenomeCoords)

        // as above, offset is 250
        offset = 250
        expectedResult = 251..4478
        pseudoGenomeCoords = getPseudoGFFCoordsMultipleRegions(gffCoords, regions, offset)
        assertEquals(expectedResult, pseudoGenomeCoords)

    }

    @Test
    fun testGetPseudoGFFCoordsMultipleRegions() {

        // simple example where first 2 in region list overlap asm,
        // but 3rd one is outside the gene boundaries.
        var gffCoords = 871..1000
        var regions = mutableListOf(800..900, 950..1100, 1200..1300)
        var offset = 0
        // Total number of bps in the gff entry is 1000 - 871 + 1  = 130 (because is inclusive/inclusive)
        // The start will be 71 (offset=0 means this is the start of the haplotype sequence) and
        // the gene starts 71 bps into this sequence.
        // The first 2 regions overlap the gff entry, so the end range of the pseudo genome entry
        // should be 100 (900-800+1) + 51 (1000-950+1) = 151
        var expectedResult = 71..151  //
        var pseudoGenomeCoords = getPseudoGFFCoordsMultipleRegions(gffCoords, regions, offset)
        assertEquals(expectedResult, pseudoGenomeCoords)

        // same as above, but offset is 50
        offset = 50
        expectedResult = 121..201
        pseudoGenomeCoords = getPseudoGFFCoordsMultipleRegions(gffCoords, regions, offset)
        assertEquals(expectedResult, pseudoGenomeCoords)

        // only the last region is within the gffCoords
        regions = mutableListOf(50..100, 150..300, 800..1300)
        offset = 0
        // The last region is the only one that overlaps the gffCoords and it overlaps
        // all of it.  So the size will be 1000-871+1 = 130
        // There is no offset (this is the beginning of the haplotype equence) but there
        // are 2 regions of sequence before the gene starts.
        // SO the start is 100-50+1 = 51, plus 300-150 + 1 = 151 = 202, plus 871-800 = 71 (it starts on 871 so don't add 1)
        // Start is hence 51+151+71 = 273
        // end is 273+130(size of gene) -1 = 402
        expectedResult = 273..402
        pseudoGenomeCoords = getPseudoGFFCoordsMultipleRegions(gffCoords, regions, offset)
        assertEquals(expectedResult, pseudoGenomeCoords)

    }

        @Test
    fun testCountGffEntriesPerChrom() {

        // Michelle would like to create the gff's in memory, then count
        // how many gffs entries occur for each chromosome, and do some other
        // metrics.  THis test case verifies the call to utilities that create
        // these metrics

        val taxonGffMap = createTaxonToGFFMap(gffFileWithHeadersB73,gffFileWithHeadersCML103)

        // Get the number of ids that start with Zm00021ab000040 from the taxonGffMap for CML103

        // Checking CMl103 for Zm00021ab000040 - should be 5
        var cmlGffSet = taxonGffMap.get("CML103")!!.get(0)
        var idToMatch = "Zm00021ab000040"
        var idCount = gffSingleIDcount(idToMatch, cmlGffSet)
        // CML gff file has 1 gene, 1 mRNA, 1 exon, 1 CDS and 1 three-prime-UTR entry for Zm00021ab000040
        assertEquals(5,idCount)

        idToMatch ="Zm00021ab000090"
        idCount = gffSingleIDcount(idToMatch, cmlGffSet)
        // the CML gff file has 1 gene, 1 mRN. 4 exon and 4 CDS entries for  Zm00021ab000090
        assertEquals(10,idCount)

        // CML gff file has 36 chr1 entries and 20 chr2 entries
        var chromCount = getGFFEntriesPerChrom(cmlGffSet)
        assertEquals(36,chromCount.get("chr1"))
        assertEquals(20, chromCount.get("chr2"))

        var B73GffSet = taxonGffMap.get("B73")!!.get(0)
        idToMatch = "Zm00001e000002"
        idCount = gffSingleIDcount(idToMatch,B73GffSet)
        assertEquals(34, idCount)

        idToMatch = "Zm00001e006484"
        idCount = gffSingleIDcount(idToMatch,B73GffSet)
        assertEquals(6,idCount)

        // Check for an ID that doesn't exist in the gff file
        idToMatch = "RS23481"
        idCount = gffSingleIDcount(idToMatch,B73GffSet)
        assertEquals(0,idCount)

    }

    @Test   fun testSumPerChromBPs() {
        val gffFeatures = readGFFtoGff3Feature(gffFileShortBPTesting).toSet()
        var chromCount = sumPerChromGFFBasePairs(gffFeatures)
        TestCase.assertTrue(chromCount.containsKey("chr1"))
        TestCase.assertTrue(chromCount.containsKey("chr2"))

        assertEquals(87, chromCount["chr1"])
        assertEquals(29,chromCount["chr2"])

        chromCount = sumPerChromNonGFFBasePairs(gffFeatures)
        TestCase.assertTrue(chromCount.containsKey("chr1"))
        TestCase.assertTrue(chromCount.containsKey("chr2"))

        assertEquals(113, chromCount["chr1"])
        assertEquals(71,chromCount["chr2"])
    }

    @Test
    fun testCountingRanges() {

        val chromRangeSet: RangeSet<Position> = TreeRangeSet.create()

        // Create the ranges for testing
        var pos1 = Position("chr1", 6)
        var pos2 = Position("chr1", 22)

        chromRangeSet.add(Range.closed(pos1, pos2))
        pos1 = Position("chr1", 18)
        pos2 = Position("chr1", 28)
        chromRangeSet.add(Range.closed(pos1, pos2))

        pos1 = Position("chr2", 20)
        pos2 = Position("chr2", 28)
        chromRangeSet.add(Range.closed(pos1, pos2))
        chromRangeSet.add(Range.closed(Position("chr2", 10), Position("chr2", 14)))

        var testRange1 = Range.closed(Position("chr1", 1), Position("chr1", 100))
        var testRange2 = Range.closed(Position("chr2", 1), Position("chr2", 100))
        val subRange1 = chromRangeSet.subRangeSet(testRange1)
        val subRange2 = chromRangeSet.subRangeSet(testRange2)
        println("subRange1: ${subRange1.toString()}")
        println("\nsubrange2: ${subRange2}")

        val includedBP = chromRangeSet.subRangeSet(testRange2).asRanges()
            .map { range ->
                range.upperEndpoint().position - range.lowerEndpoint().position +1
            }.sum()

        assertEquals(14, includedBP)
    }

    @Test
    fun testBPpercentages() {

        val gffFeatures = readGFFtoGff3Feature(gffFileShortBPTesting).toSet()
        val chromCount = sumPerChromGFFBasePairs(gffFeatures)
        TestCase.assertTrue(chromCount.containsKey("chr1"))
        TestCase.assertTrue(chromCount.containsKey("chr2"))

        assertEquals(87, chromCount["chr1"])
        assertEquals(29,chromCount["chr2"])

        var chromosomes = chromCount.keys.toList()
        println("chromosomes keys as list:\n${chromosomes}")
        val counts = chromCount.values.toList()
        println("counts as list:\n${counts}")
        val mapData = mapOf("chromosome" to chromosomes,
            "count" to counts)


        val pctBPIncludedInGFF = percentPerChromGFFBasePairs(gffFeatures)
        val percentIn = pctBPIncludedInGFF.values.toList()

        assertEquals(0.435,pctBPIncludedInGFF["chr1"]!!,0.00001)
        assertEquals(0.29,pctBPIncludedInGFF["chr2"]!!, 0.00001)

        val pctBPNOTIncludedInGFF = percentPerChromNonGFFBasePairs(gffFeatures)
        val percentOut = pctBPNOTIncludedInGFF.values.toList()

        assertEquals(0.565,pctBPNOTIncludedInGFF["chr1"]!!,0.00001)
        assertEquals(0.71,pctBPNOTIncludedInGFF["chr2"]!!,0.00001)

    }

}

fun createGFFKeyfile(gffKeyFile:String, gffFileWithHeadersB73:String, gffFileWithHeadersCML103:String) {
    // write a keyfile to go with the gff file for these tests
    try {
        bufferedWriter(gffKeyFile).use { bw ->
            bw.write("B73\t${gffFileWithHeadersB73}\n")
            bw.write("CML103\t${gffFileWithHeadersCML103}")
        }

    } catch (exc:Exception) {
        throw IllegalStateException(" createGFFKeyfile: error writing test gff key file: " + exc.message)
    }
}

fun createTaxonToGFFMap(gffB73:String, gffCML103:String): Map<String,List<Set<Gff3Feature>>> {
    val gffFeatureSetB73 = readGFFtoGff3Feature(gffB73).toSet()
    val gffFeatureSetCML103 = readGFFtoGff3Feature(gffCML103).toSet()
    val taxonGffMap = mutableMapOf<String,List<Set<Gff3Feature>>>()
    taxonGffMap.put("B73",listOf(gffFeatureSetB73))
    taxonGffMap.put("CML103",listOf(gffFeatureSetCML103))
    return taxonGffMap
}

// This file is very short, will be used for summing and finding percentage
// of BPs covered for each chromosome
// For chrom 1: the chrom range is 1-200
//      The covered BPs are 4-37, 4-37, 26-50, 120-159
// For Chrom 2: the chrom range is 1-100
//      The covered BPs are 60-79, 60-79, 62-80, 85-92
fun createGFFFileForBPTests(gffFile:String) {
    println("createGFFFileForBPTests - writing to file ${gffFile}")
    try {
        bufferedWriter(gffFile).use { bw ->
            bw.write("##gff-version\t3\n")
            bw.write("chr1\tassembly\tchromosome\t1\t200\t.\t.\t.\tID=1;Name=chromosome:Zm-B73-REFERENCE-NAM-5.0:1:1:308452471:1\n")
            bw.write("chr1\tNAM\tgene\t4\t37\t.\t-\t.\tID=gene:Zm00001e000002;biotype=protein_coding;logic_name=mikado_gene\n")
            bw.write("chr1\tNAM\tmRNA\t4\t37\t.\t-\t.\tID=transcript:Zm00001e000002_T001;Parent=gene:Zm00001e000002;biotype=protein_coding;transcript_id=Zm00001e000002_T001;_AED=0.22;_QI=2|1|0.75|1|0|0|8|313|357\n")
            bw.write("chr1\tNAM\tthree_prime_UTR\t26\t50\t.\t-\t.\tParent=transcript:Zm00001e000002_T001\n")
            bw.write("chr1\tNAM\tgene\t120\t159\t.\t-\t.\tID=gene:Zm00001e006481;biotype=protein_coding;logic_name=mikado_gene\n")
            bw.write("chr2\tassembly\tchromosome\t1\t100\t.\t.\t.\tID=2;Name=chromosome:Zm-B73-REFERENCE-NAM-5.0:2:1:243675191:1\n")
            bw.write("chr2\tNAM\tgene\t60\t79\t.\t-\t.\tID=gene:Zm00001e006478;biotype=protein_coding;logic_name=mikado_gene\n")
            bw.write("chr2\tNAM\tmRNA\t60\t79\t.\t-\t.\tID=transcript:Zm00001e006478_T001;Parent=gene:Zm00001e006478;biotype=protein_coding;transcript_id=Zm00001e006478_T001;canonical_transcript=1;_AED=0.01;_QI=188|-1|0|1|-1|0|1|439|289;Dbxref=InterPro:IPR021148,Pfam:PF04669\n")
            bw.write("chr2\tNAM\tthree_prime_UTR\t62\t80\t.\t-\t.\tParent=transcript:Zm00001e006478_T001\n")
            bw.write("chr2\tNAM\tgene\t85\t92\t.\t-\t.\tID=gene:Zm00001e0030303;biotype=protein_coding;logic_name=denver_gene\n")
        }

    } catch (exc:Exception) {
        throw IllegalStateException("Error creating file ${gffFile}: ${exc.message}")
    }
}

fun createGFFFileWithHeadersB73(gffFileWithHeaders:String) {

    // Write a small test gff file.   These values were taken from a
    // B73 gff version with maybe a few changes.
    println("createGFFFileWithHeadersB73: writing to file ${gffFileWithHeaders}")
    try {
        bufferedWriter(gffFileWithHeaders).use { bw ->
            // write header data

            bw.write("##gff-version\t3\n")
            bw.write("##sequence-region\tchr1\t1\t308452471\n")
            bw.write("##sequence-region\tchr2\t1\t243675191\n")
            bw.write("chr1\tassembly\tchromosome\t1\t308452471\t.\t.\t.\tID=1;Name=chromosome:Zm-B73-REFERENCE-NAM-5.0:1:1:308452471:1\n")
            bw.write("###\n")
            bw.write("chr1\tNAM\tgene\t41214\t46762\t.\t-\t.\tID=gene:Zm00001e000002;biotype=protein_coding;logic_name=mikado_gene\n")
            bw.write("chr1\tNAM\tmRNA\t41214\t43902\t.\t-\t.\tID=transcript:Zm00001e000002_T001;Parent=gene:Zm00001e000002;biotype=protein_coding;transcript_id=Zm00001e000002_T001;_AED=0.22;_QI=2|1|0.75|1|0|0|8|313|357\n")
            bw.write("chr1\tNAM\tthree_prime_UTR\t41214\t41526\t.\t-\t.\tParent=transcript:Zm00001e000002_T001\n")
            bw.write("chr1\tNAM\texon\t41214\t41588\t.\t-\t.\tParent=transcript:Zm00001e000002_T001;Name=Zm00001e000002_T001.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00001e000002_T001.exon.1;rank=8\n")
            bw.write("chr1\tNAM\tCDS\t41527\t41588\t.\t-\t2\tID=CDS:Zm00001e000002_P001;Parent=transcript:Zm00001e000002_T001;protein_id=Zm00001e000002_P001\n")
            bw.write("chr1\tNAM\texon\t41742\t41807\t.\t-\t.\tParent=transcript:Zm00001e000002_T001;Name=Zm00001e000002_T001.exon.2;ensembl_end_phase=1;ensembl_phase=1;exon_id=Zm00001e000002_T001.exon.2;rank=7\n")
            bw.write("chr1\tNAM\tCDS\t41742\t41807\t.\t-\t2\tID=CDS:Zm00001e000002_P001;Parent=transcript:Zm00001e000002_T001;protein_id=Zm00001e000002_P001\n")
            bw.write("chr1\tNAM\texon\t42257\t42375\t.\t-\t.\tParent=transcript:Zm00001e000002_T001;Name=Zm00001e000002_T001.exon.3;ensembl_end_phase=1;ensembl_phase=1;exon_id=Zm00001e000002_T001.exon.3;rank=6\n")
            bw.write("chr1\tNAM\tCDS\t42257\t42375\t.\t-\t1\tID=CDS:Zm00001e000002_P001;Parent=transcript:Zm00001e000002_T001;protein_id=Zm00001e000002_P001\n")
            bw.write("chr1\tNAM\texon\t42508\t42665\t.\t-\t.\tParent=transcript:Zm00001e000002_T001;Name=Zm00001e000002_T001.exon.4;ensembl_end_phase=2;ensembl_phase=2;exon_id=Zm00001e000002_T001.exon.4;rank=5\n")
            bw.write("chr1\tNAM\tCDS\t42508\t42665\t.\t-\t0\tID=CDS:Zm00001e000002_P001;Parent=transcript:Zm00001e000002_T001;protein_id=Zm00001e000002_P001\n")
            bw.write("chr1\tNAM\texon\t42762\t42917\t.\t-\t.\tParent=transcript:Zm00001e000002_T001;Name=Zm00001e000002_T001.exon.5;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e000002_T001.exon.5;rank=4\n")
            bw.write("chr1\tNAM\tCDS\t42762\t42917\t.\t-\t0\tID=CDS:Zm00001e000002_P001;Parent=transcript:Zm00001e000002_T001;protein_id=Zm00001e000002_P001\n")
            bw.write("chr1\tNAM\texon\t43039\t43197\t.\t-\t.\tParent=transcript:Zm00001e000002_T001;Name=Zm00001e000002_T001.exon.6;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e000002_T001.exon.6;rank=3\n")
            bw.write("chr1\tNAM\tCDS\t43039\t43197\t.\t-\t0\tID=CDS:Zm00001e000002_P001;Parent=transcript:Zm00001e000002_T001;protein_id=Zm00001e000002_P001\n")
            bw.write("chr1\tNAM\texon\t43318\t43520\t.\t-\t.\tParent=transcript:Zm00001e000002_T001;Name=Zm00001e000002_T001.exon.7;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e000002_T001.exon.7;rank=2\n")
            bw.write("chr1\tNAM\tCDS\t43318\t43520\t.\t-\t2\tID=CDS:Zm00001e000002_P001;Parent=transcript:Zm00001e000002_T001;protein_id=Zm00001e000002_P001\n")
            bw.write("chr1\tNAM\tCDS\t43750\t43900\t.\t-\t0\tID=CDS:Zm00001e000002_P001;Parent=transcript:Zm00001e000002_T001;protein_id=Zm00001e000002_P001\n")
            bw.write("chr1\tNAM\texon\t43750\t43902\t.\t-\t.\tParent=transcript:Zm00001e000002_T001;Name=Zm00001e000002_T001.exon.8;ensembl_end_phase=1;ensembl_phase=1;exon_id=Zm00001e000002_T001.exon.8;rank=1\n")
            bw.write("chr1\tNAM\tfive_prime_UTR\t43901\t43902\t.\t-\t.\tParent=transcript:Zm00001e000002_T001\n")
            bw.write("chr1\tNAM\tmRNA\t41214\t46762\t.\t-\t.\tID=transcript:Zm00001e000002_T002;Parent=gene:Zm00001e000002;biotype=protein_coding;transcript_id=Zm00001e000002_T002;canonical_transcript=1;_AED=0.20;_QI=843|1|0.77|1|0|0|9|313|519;Dbxref=InterPro:IPR016072,InterPro:IPR016073,Pfam:PF01466,Pfam:PF03931;Ontology_term=GO:0006511\n")
            bw.write("chr1\tNAM\tthree_prime_UTR\t41214\t41526\t.\t-\t.\tParent=transcript:Zm00001e000002_T002\n")
            bw.write("chr1\tNAM\texon\t41214\t41588\t.\t-\t.\tParent=transcript:Zm00001e000002_T002;Name=Zm00001e000002_T002.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00001e000002_T002.exon.1;rank=9\n")
            bw.write("chr1\tNAM\tCDS\t41527\t41588\t.\t-\t2\tID=CDS:Zm00001e000002_P002;Parent=transcript:Zm00001e000002_T002;protein_id=Zm00001e000002_P002\n")
            bw.write("chr1\tNAM\texon\t41742\t41807\t.\t-\t.\tParent=transcript:Zm00001e000002_T002;Name=Zm00001e000002_T002.exon.2;ensembl_end_phase=1;ensembl_phase=1;exon_id=Zm00001e000002_T002.exon.2;rank=8\n")
            bw.write("chr1\tNAM\tCDS\t41742\t41807\t.\t-\t2\tID=CDS:Zm00001e000002_P002;Parent=transcript:Zm00001e000002_T002;protein_id=Zm00001e000002_P002\n")
            bw.write("chr1\tNAM\texon\t42257\t42375\t.\t-\t.\tParent=transcript:Zm00001e000002_T002;Name=Zm00001e000002_T002.exon.3;ensembl_end_phase=1;ensembl_phase=1;exon_id=Zm00001e000002_T002.exon.3;rank=7\n")
            bw.write("chr1\tNAM\tCDS\t42257\t42375\t.\t-\t1\tID=CDS:Zm00001e000002_P002;Parent=transcript:Zm00001e000002_T002;protein_id=Zm00001e000002_P002\n")
            bw.write("chr1\tNAM\texon\t42508\t42665\t.\t-\t.\tParent=transcript:Zm00001e000002_T002;Name=Zm00001e000002_T002.exon.4;ensembl_end_phase=2;ensembl_phase=2;exon_id=Zm00001e000002_T002.exon.4;rank=6\n")
            bw.write("chr1\tNAM\tCDS\t42508\t42665\t.\t-\t0\tID=CDS:Zm00001e000002_P002;Parent=transcript:Zm00001e000002_T002;protein_id=Zm00001e000002_P002\n")
            bw.write("chr1\tNAM\texon\t42762\t42917\t.\t-\t.\tParent=transcript:Zm00001e000002_T002;Name=Zm00001e000002_T002.exon.5;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e000002_T002.exon.5;rank=5\n")
            bw.write("chr1\tNAM\tCDS\t42762\t42917\t.\t-\t0\tID=CDS:Zm00001e000002_P002;Parent=transcript:Zm00001e000002_T002;protein_id=Zm00001e000002_P002\n")
            bw.write("chr1\tNAM\texon\t43039\t43197\t.\t-\t.\tParent=transcript:Zm00001e000002_T002;Name=Zm00001e000002_T002.exon.6;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e000002_T002.exon.6;rank=4\n")
            bw.write("chr1\tNAM\tCDS\t43039\t43197\t.\t-\t0\tID=CDS:Zm00001e000002_P002;Parent=transcript:Zm00001e000002_T002;protein_id=Zm00001e000002_P002;product=lynnProduct1,MichelleProduct2\n")
            bw.write("chr1\tNAM\tTE\t206700\t206800\t.\t-\t0\tID=TE:Lynn_TE1;Name=Copia\n")
            bw.write("chr2\tassembly\tchromosome\t1\t243675191\t.\t.\t.\tID=2;Name=chromosome:Zm-B73-REFERENCE-NAM-5.0:2:1:243675191:1\n")
            bw.write("###\n")
            bw.write("chr2\tNAM\tgene\t29923\t31419\t.\t-\t.\tID=gene:Zm00001e006478;biotype=protein_coding;logic_name=mikado_gene\n")
            bw.write("chr2\tNAM\tmRNA\t29923\t31419\t.\t-\t.\tID=transcript:Zm00001e006478_T001;Parent=gene:Zm00001e006478;biotype=protein_coding;transcript_id=Zm00001e006478_T001;canonical_transcript=1;_AED=0.01;_QI=188|-1|0|1|-1|0|1|439|289;Dbxref=InterPro:IPR021148,Pfam:PF04669\n")
            bw.write("chr2\tNAM\tthree_prime_UTR\t29923\t30361\t.\t-\t.\tParent=transcript:Zm00001e006478_T001\n")
            bw.write("chr2\tNAM\texon\t29923\t31419\t.\t-\t.\tParent=transcript:Zm00001e006478_T001;Name=Zm00001e006478_T001.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00001e006478_T001.exon.1;rank=1\n")
            bw.write("chr2\tNAM\tCDS\t30362\t31231\t.\t-\t0\tID=CDS:Zm00001e006478_P001;Parent=transcript:Zm00001e006478_T001;protein_id=Zm00001e006478_P001\n")
            bw.write("chr2\tNAM\tfive_prime_UTR\t31232\t31419\t.\t-\t.\tParent=transcript:Zm00001e006478_T001\n")
            bw.write("###\n")
            bw.write("chr2\tNAM\tgene\t37920\t41238\t.\t+\t.\tID=gene:Zm00001e006479;biotype=protein_coding;logic_name=mikado_gene\n")
            bw.write("chr2\tNAM\tmRNA\t37920\t41238\t.\t+\t.\tID=transcript:Zm00001e006479_T001;Parent=gene:Zm00001e006479;biotype=protein_coding;transcript_id=Zm00001e006479_T001;canonical_transcript=1;_AED=1.00;_QI=0|0|0|0|0|0|2|2|1082\n")
            bw.write("chr2\tNAM\texon\t37920\t41090\t.\t+\t.\tParent=transcript:Zm00001e006479_T001;Name=Zm00001e006479_T001.exon.1;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e006479_T001.exon.1;rank=1\n")
            bw.write("chr2\tNAM\tCDS\t37920\t41090\t.\t+\t0\tID=CDS:Zm00001e006479_P001;Parent=transcript:Zm00001e006479_T001;protein_id=Zm00001e006479_P001\n")
            bw.write("chr2\tNAM\tCDS\t41162\t41236\t.\t+\t0\tID=CDS:Zm00001e006479_P001;Parent=transcript:Zm00001e006479_T001;protein_id=Zm00001e006479_P001\n")
            bw.write("chr2\tNAM\texon\t41162\t41238\t.\t+\t.\tParent=transcript:Zm00001e006479_T001;Name=Zm00001e006479_T001.exon.2;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00001e006479_T001.exon.2;rank=2\n")
            bw.write("chr2\tNAM\tthree_prime_UTR\t41237\t41238\t.\t+\t.\tParent=transcript:Zm00001e006479_T001\n")
            bw.write("###\n")
            bw.write("chr2\tNAM\tgene\t166443\t167966\t.\t-\t.\tID=gene:Zm00001e006481;biotype=protein_coding;logic_name=mikado_gene\n")
            bw.write("chr2\tNAM\tmRNA\t166443\t167966\t.\t-\t.\tID=transcript:Zm00001e006481_T001;Parent=gene:Zm00001e006481;biotype=protein_coding;transcript_id=Zm00001e006481_T001;canonical_transcript=1;_AED=0.25;_QI=1|0|0|1|0|0|2|663|252\n")
            bw.write("chr2\tNAM\tthree_prime_UTR\t166443\t167105\t.\t-\t.\tParent=transcript:Zm00001e006481_T001\n")
            bw.write("chr2\tNAM\texon\t166443\t167109\t.\t-\t.\tParent=transcript:Zm00001e006481_T001;Name=Zm00001e006481_T001.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00001e006481_T001.exon.1;rank=2\n")
            bw.write("chr2\tNAM\tCDS\t167106\t167109\t.\t-\t1\tID=CDS:Zm00001e006481_P001;Parent=transcript:Zm00001e006481_T001;protein_id=Zm00001e006481_P001\n")
            bw.write("chr2\tNAM\tCDS\t167211\t167965\t.\t-\t0\tID=CDS:Zm00001e006481_P001;Parent=transcript:Zm00001e006481_T001;protein_id=Zm00001e006481_P001\n")
            bw.write("chr2\tNAM\texon\t167211\t167966\t.\t-\t.\tParent=transcript:Zm00001e006481_T001;Name=Zm00001e006481_T001.exon.2;ensembl_end_phase=2;ensembl_phase=2;exon_id=Zm00001e006481_T001.exon.2;rank=1\n")
            bw.write("chr2\tNAM\tfive_prime_UTR\t167966\t167966\t.\t-\t.\tParent=transcript:Zm00001e006481_T001\n")
            bw.write("###\n")
            bw.write("chr2\tNAM\tgene\t167108\t168713\t.\t+\t.\tID=gene:Zm00001e006482;biotype=protein_coding;logic_name=mikado_gene\n")
            bw.write("chr2\tNAM\tmRNA\t167108\t168713\t.\t+\t.\tID=transcript:Zm00001e006482_T001;Parent=gene:Zm00001e006482;biotype=protein_coding;transcript_id=Zm00001e006482_T001;canonical_transcript=1;_AED=0.28;_QI=0|1|0.5|1|0|0|2|517|336\n")
            bw.write("chr2\tNAM\texon\t167108\t167548\t.\t+\t.\tParent=transcript:Zm00001e006482_T001;Name=Zm00001e006482_T001.exon.1;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e006482_T001.exon.1;rank=1\n")
            bw.write("chr2\tNAM\tCDS\t167108\t167548\t.\t+\t0\tID=CDS:Zm00001e006482_P001;Parent=transcript:Zm00001e006482_T001;protein_id=Zm00001e006482_P001\n")
            bw.write("chr2\tNAM\tCDS\t167627\t168196\t.\t+\t0\tID=CDS:Zm00001e006482_P001;Parent=transcript:Zm00001e006482_T001;protein_id=Zm00001e006482_P001\n")
            bw.write("chr2\tNAM\texon\t167627\t168713\t.\t+\t.\tParent=transcript:Zm00001e006482_T001;Name=Zm00001e006482_T001.exon.2;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00001e006482_T001.exon.2;rank=2\n")
            bw.write("chr2\tNAM\tthree_prime_UTR\t168197\t168713\t.\t+\t.\tParent=transcript:Zm00001e006482_T001\n")
            bw.write("###\n")
            bw.write("chr2\tNAM\tgene\t249250\t250710\t.\t+\t.\tID=gene:Zm00001e006483;biotype=protein_coding;logic_name=mikado_gene\n")
            bw.write("chr2\tNAM\tmRNA\t249250\t250710\t.\t+\t.\tID=transcript:Zm00001e006483_T001;Parent=gene:Zm00001e006483;biotype=protein_coding;transcript_id=Zm00001e006483_T001;canonical_transcript=1;_AED=0.04;_QI=2|-1|0|1|-1|0|1|286|390;Dbxref=InterPro:IPR000210,InterPro:IPR002083,Pfam:PF00651,Pfam:PF00917;Ontology_term=GO:0005515\n")
            bw.write("chr2\tNAM\tfive_prime_UTR\t249250\t249251\t.\t+\t.\tParent=transcript:Zm00001e006483_T001\n")
            bw.write("chr2\tNAM\texon\t249250\t250710\t.\t+\t.\tParent=transcript:Zm00001e006483_T001;Name=Zm00001e006483_T001.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00001e006483_T001.exon.1;rank=1\n")
            bw.write("chr2\tNAM\tCDS\t249252\t250424\t.\t+\t0\tID=CDS:Zm00001e006483_P001;Parent=transcript:Zm00001e006483_T001;protein_id=Zm00001e006483_P001\n")
            bw.write("chr2\tNAM\tthree_prime_UTR\t250425\t250710\t.\t+\t.\tParent=transcript:Zm00001e006483_T001\n")
            bw.write("###\n")
            bw.write("chr2\tNAM\tgene\t342932\t344528\t.\t-\t.\tID=gene:Zm00001e006484;biotype=protein_coding;logic_name=mikado_gene\n")
            bw.write("chr2\tNAM\tmRNA\t342932\t344528\t.\t-\t.\tID=transcript:Zm00001e006484_T001;Parent=gene:Zm00001e006484;biotype=protein_coding;transcript_id=Zm00001e006484_T001;canonical_transcript=1;_AED=0.07;_QI=113|-1|0|1|-1|0|1|2|494\n")
            bw.write("chr2\tNAM\tthree_prime_UTR\t342932\t342933\t.\t-\t.\tParent=transcript:Zm00001e006484_T001\n")
            bw.write("chr2\tNAM\texon\t342932\t344528\t.\t-\t.\tParent=transcript:Zm00001e006484_T001;Name=Zm00001e006484_T001.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00001e006484_T001.exon.1;rank=1\n")
            bw.write("chr2\tNAM\tCDS\t342934\t344415\t.\t-\t0\tID=CDS:Zm00001e006484_P001;Parent=transcript:Zm00001e006484_T001;protein_id=Zm00001e006484_P001\n")
            bw.write("chr2\tNAM\tfive_prime_UTR\t344416\t344528\t.\t-\t.\tParent=transcript:Zm00001e006484_T001\n")
            bw.write("###\n")
            bw.write("chr2\tNAM\tgene\t999575\t1007836\t.\t-\t.\tID=gene:Zm00001e006492;biotype=protein_coding;logic_name=mikado_gene\n")
            bw.write("chr2\tNAM\tmRNA\t999575\t1007822\t.\t-\t.\tID=transcript:Zm00001e006492_T001;Parent=gene:Zm00001e006492;biotype=protein_coding;transcript_id=Zm00001e006492_T001;_AED=0.20;_QI=356|0.90|0.75|1|0|0|12|978|463;Dbxref=InterPro:IPR001245,Pfam:PF07714;Ontology_term=GO:0004672,GO:0006468\n")
            bw.write("chr2\tNAM\texon\t999575\t1000510\t.\t-\t.\tParent=transcript:Zm00001e006492_T001;Name=Zm00001e006492_T001.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00001e006492_T001.exon.1;rank=12\n")
            bw.write("chr2\tNAM\tthree_prime_UTR\t999575\t1000510\t.\t-\t.\tParent=transcript:Zm00001e006492_T001\n")
            bw.write("chr2\tNAM\tthree_prime_UTR\t1001043\t1001084\t.\t-\t.\tParent=transcript:Zm00001e006492_T001\n")
            bw.write("chr2\tNAM\texon\t1001043\t1001249\t.\t-\t.\tParent=transcript:Zm00001e006492_T001;Name=Zm00001e006492_T001.exon.2;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00001e006492_T001.exon.2;rank=11\n")
            bw.write("chr2\tNAM\tCDS\t1001085\t1001249\t.\t-\t0\tID=CDS:Zm00001e006492_P001;Parent=transcript:Zm00001e006492_T001;protein_id=Zm00001e006492_P001\n")
            bw.write("chr2\tNAM\texon\t1001573\t1001686\t.\t-\t.\tParent=transcript:Zm00001e006492_T001;Name=Zm00001e006492_T001.exon.3;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e006492_T001.exon.3;rank=10\n")
            bw.write("chr2\tNAM\tCDS\t1001573\t1001686\t.\t-\t0\tID=CDS:Zm00001e006492_P001;Parent=transcript:Zm00001e006492_T001;protein_id=Zm00001e006492_P001\n")
            bw.write("chr2\tNAM\texon\t1002097\t1002261\t.\t-\t.\tParent=transcript:Zm00001e006492_T001;Name=Zm00001e006492_T001.exon.4;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e006492_T001.exon.4;rank=9\n")
            bw.write("chr2\tNAM\tCDS\t1002097\t1002261\t.\t-\t0\tID=CDS:Zm00001e006492_P001;Parent=transcript:Zm00001e006492_T001;protein_id=Zm00001e006492_P001\n")
            bw.write("chr2\tNAM\texon\t1003055\t1003243\t.\t-\t.\tParent=transcript:Zm00001e006492_T001;Name=Zm00001e006492_T001.exon.5;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e006492_T001.exon.5;rank=8\n")
            bw.write("chr2\tNAM\tCDS\t1003055\t1003243\t.\t-\t0\tID=CDS:Zm00001e006492_P001;Parent=transcript:Zm00001e006492_T001;protein_id=Zm00001e006492_P001\n")
            bw.write("chr2\tNAM\texon\t1004129\t1004214\t.\t-\t.\tParent=transcript:Zm00001e006492_T001;Name=Zm00001e006492_T001.exon.6;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e006492_T001.exon.6;rank=7\n")
            bw.write("chr2\tNAM\tCDS\t1004129\t1004214\t.\t-\t2\tID=CDS:Zm00001e006492_P001;Parent=transcript:Zm00001e006492_T001;protein_id=Zm00001e006492_P001\n")
            bw.write("chr2\tNAM\texon\t1004818\t1004923\t.\t-\t.\tParent=transcript:Zm00001e006492_T001;Name=Zm00001e006492_T001.exon.7;ensembl_end_phase=1;ensembl_phase=1;exon_id=Zm00001e006492_T001.exon.7;rank=6\n")
            bw.write("chr2\tNAM\tCDS\t1004818\t1004923\t.\t-\t0\tID=CDS:Zm00001e006492_P001;Parent=transcript:Zm00001e006492_T001;protein_id=Zm00001e006492_P001\n")
            bw.write("chr2\tNAM\texon\t1005003\t1005136\t.\t-\t.\tParent=transcript:Zm00001e006492_T001;Name=Zm00001e006492_T001.exon.8;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001e006492_T001.exon.8;rank=5\n")
            bw.write("chr2\tNAM\tTE\t1005303\t1006303\t.\t-\t0\tID=TE:Lynn_TE22;Name=Copia\n")


        }

    } catch (exc: Exception) {
        throw IllegalStateException(" createGFFFileWithHeaders: error writing test gff file: " + exc.message)
    }
}
fun createGFFFileWithHeadersCML103(gffFileWithHeaders:String) {

    // Write a small test gff file.   These values were taken from a
    // CML103 GFF from Michelle used when testing PathsToGFFPlugin
    println("createGFFFileWithHeadersCML103: writing to file ${gffFileWithHeaders}")
    try {
        bufferedWriter(gffFileWithHeaders).use { bw ->
            bw.write("##gff-version\t3\n")
            bw.write("#\tgenerated\ton\tWed\tAug\t26\t07:59:04\t2020\tby\t/sonas-hs/ware/hpc/home/olson/src/warelab/gramene-ensembl/scripts/dump-scripts/sharon_gff_dump_by_logicname_nam.pl\n")
            bw.write("#\tfor\tspecies\tMaize\tcml103\n")
            bw.write("#\tgenebuild\t2019-cshl\n")
            bw.write("chr1\tassembly\tchromosome\t1\t307163243\t.\t.\t.\tID=1;Name=chromosome:Zm-CML103-REFERENCE-NAM-1.0:1:1:307163243:1\n")
            bw.write("chr1\tNAM\tgene\t97153\t97537\t.\t+\t.\tID=Zm00021ab000040;biotype=protein_coding;logic_name=cshl_gene\n")
            bw.write("chr1\tNAM\tmRNA\t97153\t97537\t.\t+\t.\tID=Zm00021ab000040_T001;Parent=Zm00021ab000040;biotype=protein_coding;transcript_id=Zm00021ab000040_T001;canonical_transcript=1\n")
            bw.write("chr1\tNAM\texon\t97153\t97537\t.\t+\t.\tParent=Zm00021ab000040_T001;Name=Zm00021ab000040_T001.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00021ab000040_T001.exon.1;rank=1\n")
            bw.write("chr1\tNAM\tCDS\t97153\t97458\t.\t+\t0\tID=Zm00021ab000040_P001;Parent=Zm00021ab000040_T001;protein_id=Zm00021ab000040_P001\n")
            bw.write("chr1\tNAM\tthree_prime_UTR\t97459\t97537\t.\t+\t.\tParent=Zm00021ab000040_T001\n")
            bw.write("chr1\tNAM\tgene\t100827\t101813\t.\t-\t.\tID=Zm00021ab000050;biotype=protein_coding;logic_name=cshl_gene\n")
            bw.write("chr1\tNAM\tmRNA\t100827\t101813\t.\t-\t.\tID=Zm00021ab000050_T001;Parent=Zm00021ab000050;biotype=protein_coding;transcript_id=Zm00021ab000050_T001;canonical_transcript=1\n")
            bw.write("chr1\tNAM\tfive_prime_UTR\t101701\t101813\t.\t-\t.\tParent=Zm00021ab000050_T001\n")
            bw.write("chr1\tNAM\texon\t100827\t101813\t.\t-\t.\tParent=Zm00021ab000050_T001;Name=Zm00021ab000050_T001.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00021ab000050_T001.exon.1;rank=1\n")
            bw.write("chr1\tNAM\tCDS\t101095\t101700\t.\t-\t0\tID=Zm00021ab000050_P001;Parent=Zm00021ab000050_T001;protein_id=Zm00021ab000050_P001\n")
            bw.write("chr1\tNAM\tthree_prime_UTR\t100827\t101094\t.\t-\t.\tParent=Zm00021ab000050_T001\n")
            bw.write("chr1\tNAM\tgene\t109169\t109528\t.\t+\t.\tID=Zm00021ab000070;biotype=protein_coding;logic_name=cshl_gene\n")
            bw.write("chr1\tNAM\tmRNA\t109169\t109528\t.\t+\t.\tID=Zm00021ab000070_T001;Parent=Zm00021ab000070;biotype=protein_coding;transcript_id=Zm00021ab000070_T001;canonical_transcript=1\n")
            bw.write("chr1\tNAM\texon\t109169\t109528\t.\t+\t.\tParent=Zm00021ab000070_T001;Name=Zm00021ab000070_T001.exon.1;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00021ab000070_T001.exon.1;rank=1\n")
            bw.write("chr1\tNAM\tCDS\t109169\t109528\t.\t+\t0\tID=Zm00021ab000070_P001;Parent=Zm00021ab000070_T001;protein_id=Zm00021ab000070_P001\n")
            bw.write("chr1\tNAM\tgene\t174601\t176014\t.\t-\t.\tID=Zm00021ab000090;biotype=protein_coding;logic_name=cshl_gene\n")
            bw.write("chr1\tNAM\tmRNA\t174601\t176014\t.\t-\t.\tID=Zm00021ab000090_T001;Parent=Zm00021ab000090;biotype=protein_coding;transcript_id=Zm00021ab000090_T001;canonical_transcript=1\n")
            bw.write("chr1\tNAM\texon\t174601\t174608\t.\t-\t.\tParent=Zm00021ab000090_T001;Name=Zm00021ab000090_T001.exon.1;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00021ab000090_T001.exon.1;rank=4\n")
            bw.write("chr1\tNAM\texon\t174710\t174970\t.\t-\t.\tParent=Zm00021ab000090_T001;Name=Zm00021ab000090_T001.exon.2;ensembl_end_phase=1;ensembl_phase=1;exon_id=Zm00021ab000090_T001.exon.2;rank=3\n")
            bw.write("chr1\tNAM\texon\t175857\t175880\t.\t-\t.\tParent=Zm00021ab000090_T001;Name=Zm00021ab000090_T001.exon.3;ensembl_end_phase=1;ensembl_phase=1;exon_id=Zm00021ab000090_T001.exon.3;rank=2\n")
            bw.write("chr1\tNAM\texon\t176008\t176014\t.\t-\t.\tParent=Zm00021ab000090_T001;Name=Zm00021ab000090_T001.exon.4;ensembl_end_phase=1;ensembl_phase=1;exon_id=Zm00021ab000090_T001.exon.4;rank=1\n")
            bw.write("chr1\tNAM\tCDS\t174601\t174608\t.\t-\t2\tID=Zm00021ab000090_P001;Parent=Zm00021ab000090_T001;protein_id=Zm00021ab000090_P001\n")
            bw.write("chr1\tNAM\tCDS\t174710\t174970\t.\t-\t2\tID=Zm00021ab000090_P001;Parent=Zm00021ab000090_T001;protein_id=Zm00021ab000090_P001\n")
            bw.write("chr1\tNAM\tCDS\t175857\t175880\t.\t-\t2\tID=Zm00021ab000090_P001;Parent=Zm00021ab000090_T001;protein_id=Zm00021ab000090_P001\n")
            bw.write("chr1\tNAM\tCDS\t176008\t176014\t.\t-\t0\tID=Zm00021ab000090_P001;Parent=Zm00021ab000090_T001;protein_id=Zm00021ab000090_P001\n")
            bw.write("chr1\tNAM\tgene\t189498\t190330\t.\t-\t.\tID=Zm00021ab000110;biotype=protein_coding;logic_name=cshl_gene\n")
            bw.write("chr1\tNAM\tmRNA\t189498\t190330\t.\t-\t.\tID=Zm00021ab000110_T001;Parent=Zm00021ab000110;biotype=protein_coding;transcript_id=Zm00021ab000110_T001;canonical_transcript=1\n")
            bw.write("chr1\tNAM\tfive_prime_UTR\t189957\t190059\t.\t-\t.\tParent=Zm00021ab000110_T001\n")
            bw.write("chr1\tNAM\tfive_prime_UTR\t190140\t190330\t.\t-\t.\tParent=Zm00021ab000110_T001\n")
            bw.write("chr1\tNAM\texon\t189498\t189636\t.\t-\t.\tParent=Zm00021ab000110_T001;Name=Zm00021ab000110_T001.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00021ab000110_T001.exon.1;rank=3\n")
            bw.write("chr1\tNAM\texon\t189879\t190059\t.\t-\t.\tParent=Zm00021ab000110_T001;Name=Zm00021ab000110_T001.exon.2;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00021ab000110_T001.exon.2;rank=2\n")
            bw.write("chr1\tNAM\texon\t190140\t190330\t.\t-\t.\tParent=Zm00021ab000110_T001;Name=Zm00021ab000110_T001.exon.3;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00021ab000110_T001.exon.3;rank=1\n")
            bw.write("chr1\tNAM\tCDS\t189508\t189636\t.\t-\t0\tID=Zm00021ab000110_P001;Parent=Zm00021ab000110_T001;protein_id=Zm00021ab000110_P001\n")
            bw.write("chr1\tNAM\tCDS\t189879\t189956\t.\t-\t0\tID=Zm00021ab000110_P001;Parent=Zm00021ab000110_T001;protein_id=Zm00021ab000110_P001\n")
            bw.write("chr1\tNAM\tthree_prime_UTR\t189498\t189507\t.\t-\t.\tParent=Zm00021ab000110_T001\n")
            bw.write("chr2\tNAM\tgene\t17451\t19030\t.\t-\t.\tID=Zm00021ab065850;biotype=protein_coding;logic_name=cshl_gene\n")
            bw.write("chr2\tNAM\tmRNA\t17451\t19030\t.\t-\t.\tID=Zm00021ab065850_T001;Parent=Zm00021ab065850;biotype=protein_coding;transcript_id=Zm00021ab065850_T001;canonical_transcript=1\n")
            bw.write("chr2\tNAM\tfive_prime_UTR\t18800\t19030\t.\t-\t.\tParent=Zm00021ab065850_T001\n")
            bw.write("chr2\tNAM\texon\t17451\t19030\t.\t-\t.\tParent=Zm00021ab065850_T001;Name=Zm00021ab065850_T001.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00021ab065850_T001.exon.1;rank=1\n")
            bw.write("chr2\tNAM\tCDS\t17930\t18799\t.\t-\t0\tID=Zm00021ab065850_P001;Parent=Zm00021ab065850_T001;protein_id=Zm00021ab065850_P001\n")
            bw.write("chr2\tNAM\tthree_prime_UTR\t17451\t17929\t.\t-\t.\tParent=Zm00021ab065850_T001\n")
            bw.write("chr2\tNAM\tgene\t139960\t142935\t.\t-\t.\tID=Zm00021ab065870;biotype=protein_coding;logic_name=cshl_gene\n")
            bw.write("chr2\tNAM\tmRNA\t139960\t142935\t.\t-\t.\tID=Zm00021ab065870_T001;Parent=Zm00021ab065870;biotype=protein_coding;transcript_id=Zm00021ab065870_T001;canonical_transcript=1\n")
            bw.write("chr2\tNAM\tfive_prime_UTR\t142700\t142935\t.\t-\t.\tParent=Zm00021ab065870_T001\n")
            bw.write("chr2\tNAM\texon\t139960\t140467\t.\t-\t.\tParent=Zm00021ab065870_T001;Name=Zm00021ab065870_T001.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00021ab065870_T001.exon.1;rank=2\n")
            bw.write("chr2\tNAM\texon\t140612\t142935\t.\t-\t.\tParent=Zm00021ab065870_T001;Name=Zm00021ab065870_T001.exon.2;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00021ab065870_T001.exon.2;rank=1\n")
            bw.write("chr2\tNAM\tCDS\t141764\t142699\t.\t-\t0\tID=Zm00021ab065870_P001;Parent=Zm00021ab065870_T001;protein_id=Zm00021ab065870_P001\n")
            bw.write("chr2\tNAM\tthree_prime_UTR\t139960\t140467\t.\t-\t.\tParent=Zm00021ab065870_T001\n")
            bw.write("chr2\tNAM\tthree_prime_UTR\t140612\t141763\t.\t-\t.\tParent=Zm00021ab065870_T001\n")
            bw.write("chr2\tNAM\tgene\t204555\t205834\t.\t+\t.\tID=Zm00021ab065880;biotype=protein_coding;logic_name=cshl_gene\n")
            bw.write("chr2\tNAM\tmRNA\t204555\t205834\t.\t+\t.\tID=Zm00021ab065880_T001;Parent=Zm00021ab065880;biotype=protein_coding;transcript_id=Zm00021ab065880_T001;canonical_transcript=1\n")
            bw.write("chr2\tNAM\tfive_prime_UTR\t204555\t204587\t.\t+\t.\tParent=Zm00021ab065880_T001\n")
            bw.write("chr2\tNAM\texon\t204555\t205834\t.\t+\t.\tParent=Zm00021ab065880_T001;Name=Zm00021ab065880_T001.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00021ab065880_T001.exon.1;rank=1\n")
            bw.write("chr2\tNAM\tCDS\t204588\t205667\t.\t+\t0\tID=Zm00021ab065880_P001;Parent=Zm00021ab065880_T001;protein_id=Zm00021ab065880_P001\n")
            bw.write("chr2\tNAM\tthree_prime_UTR\t205668\t205834\t.\t+\t.\tParent=Zm00021ab065880_T001\n")
        }

    } catch (exc: Exception) {
        throw IllegalStateException(" createGFFFileWithHeaders: error writing test gff file: " + exc.message)
    }
}
