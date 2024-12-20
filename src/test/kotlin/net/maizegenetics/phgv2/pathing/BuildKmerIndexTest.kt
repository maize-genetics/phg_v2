package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.testing.test
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.cli.AgcCompress
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.getBufferedReader
import net.maizegenetics.phgv2.utils.getBufferedWriter
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Disabled
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import org.junit.jupiter.api.extension.ExtendWith
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
class BuildKmerIndexTest {
    companion object {
        //Setup/download  files
        //Resetting on both setup and teardown just to be safe.
        @JvmStatic
        @BeforeAll
        fun setup() {
            resetDirs()
            setupAgc()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            resetDirs()
        }

        fun resetDirs() {
            val tempTestDir = "${TestExtension.tempDir}kmerTest/"
            val tempHvcfDir = "${tempTestDir}hvcfDir/"
            val tempDBPathDir = "${TestExtension.testOutputFastaDir}dbPath/"

            File(TestExtension.tempDir).deleteRecursively()
            File(TestExtension.testOutputFastaDir).deleteRecursively()
            File(TestExtension.testOutputDir).deleteRecursively()
            File(tempTestDir).deleteRecursively()
            File(tempDBPathDir).deleteRecursively()
            File(tempHvcfDir).deleteRecursively()

            File(TestExtension.tempDir).mkdirs()
            File(TestExtension.testOutputFastaDir).mkdirs()
            File(TestExtension.testOutputDir).mkdirs()
            File(tempTestDir).mkdirs()
            File(tempDBPathDir).mkdirs()
            File(tempHvcfDir).mkdirs()

        }

        private fun setupAgc() {
            //create an AGC record with the Ref in it
            val altFileListFile = TestExtension.testOutputFastaDir+"/agc_altList.txt"
            BufferedWriter(FileWriter(altFileListFile)).use { writer ->
                writer.write("data/test/smallseq/LineA.fa\n")
                writer.write("data/test/smallseq/LineB.fa\n")
                writer.write("data/test/smallseq/Ref.fa\n")
            }

            val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"
            File(dbPath).mkdirs()

            //Call AGCCompress to create the AGC file
            val agcCompress = AgcCompress()
            agcCompress.processAGCFiles(dbPath,altFileListFile,"data/test/smallseq/Ref.fa","")
        }
    }

    @Test
    fun testKmerUpdating() {
        val testString = "ACACGTGTAACCGGTTGTGACTGACGGTAACGTCGAATGACGTAACCGTCGA"

        var hashValues = Pair(0L, 0L)
        repeat(32) {
            hashValues = BuildKmerIndex.updateKmerHashAndReverseCompliment(hashValues, testString[it])
        }

        // The first 32-mer of testString
        // A  C  A  C  G  T  G  T  A  A  C  C  G  G  T  T  G  T  G  A  C  T  G  A  C  G  G  T  A  A  C  G
        // 00 01 00 01 10 11 10 11 00 00 01 01 10 10 11 11 10 11 10 00 01 11 10 00 01 10 10 11 00 00 01 10

        // The reverse compliment
        // G  C  A  A  T  G  G  C  A  G  T  C  A  G  T  G  T  T  G  G  C  C  A  A  T  G  T  G  C  A  C  A
        // 01 10 11 11 00 01 01 10 11 01 00 10 11 01 00 01 00 00 01 01 10 10 11 11 00 01 00 01 10 11 10 11

        val expected = 0b0001000110111011000001011010111110111000011110000110101100000110.toLong()
        assertEquals(expected, hashValues.first, "Discrepancy in Kmer hash")

        val expectedRC = 0b0110111100010110110100101101000100000101101011110001000110111011.toLong()
        assertEquals(expectedRC, hashValues.second, "Discrepancy in reverse compliment Kmer hash")

        //test trying to convert an invalid character
        val badTestString = "ACACGTGTAACCGGTTGTMACTGACGGTAACGTCGAATGACGTAACCGTCGA"
        var badHashValues = Pair(0L, 0L)

        val exception = assertThrows<IllegalArgumentException> {
            repeat(32) {
                badHashValues = BuildKmerIndex.updateKmerHashAndReverseCompliment(badHashValues, badTestString[it])
            }
        }

        assertEquals("Attempted to update kmer hash with an invalid nucleotide character(M). Must be one of A,G,C,T", exception.message)

    }

    @Test
    fun testCliktParams() {
        //parameters:
        //db-path - required
        //maxHaplotypeProportion - defaults to 0.75
        //hashMask - default to 3L
        //hashFilterValue - defauls to 1L
        //hvcfDir - required for now

        val testBuild = BuildKmerIndex()
        val noargResult = testBuild.test("--db-path ${TestExtension.testOutputFastaDir}/dbPath")
        assertEquals(1, noargResult.statusCode)
        assertEquals("Usage: build-kmer-index [<options>]\n\n" +
                "Error: missing option --hvcf-dir\n", noargResult.stderr)


        val noDbPathResult = testBuild.test("--hvcf-dir ${TestExtension.smallseqLineAHvcfFile}")
        assertEquals(1, noDbPathResult.statusCode)
        assertEquals("Usage: build-kmer-index [<options>]\n\n" +
                "Error: missing option --db-path\n", noDbPathResult.stderr)
    }


    @Test
    fun testProcessGraphKmers() {

        //set up temporary file names
        val tempTestDir = "${TestExtension.tempDir}kmerTest/"
        val tempHvcfDir = "${tempTestDir}hvcfDir/"
        val tempDBPathDir = "${TestExtension.testOutputFastaDir}/dbPath"

        //copy hvcf files to temp directory,
        // include the ref hvcf to test what happens when samples have no haplotype in some ref range
        File(tempHvcfDir).listFiles().forEach { file -> file.delete()}
        listOf(TestExtension.smallseqLineAHvcfFile,TestExtension.smallseqRefHvcfFile)
            .forEach { hvcfFile ->
            val dst = File("$tempHvcfDir${File(hvcfFile).name}")
            if (!dst.exists()) {
                File(hvcfFile).copyTo(dst)
            }
        }

        //reverse one of the alt headers in LineB hvcf, in order to cover lines testing for inverted sequence in the alt header
        getBufferedReader(TestExtension.smallseqLineBHvcfFile).use { hvcfReader ->
            getBufferedWriter("${tempHvcfDir}LineB.h.vcf").use { hvcfWriter ->
                var inputLine = hvcfReader.readLine()
                while (inputLine != null) {
                    val outputLine = inputLine.replace("Regions=1:28501-33000","Regions=1:33000-28501")
                    hvcfWriter.write(outputLine)
                    hvcfWriter.write("\n")
                    inputLine = hvcfReader.readLine()
                }
            }
        }

        //create a HaplotypeGraph from the hvcf files
        val buildIndexResult = BuildKmerIndex().test("--db-path $tempDBPathDir --hvcf-dir $tempHvcfDir --initial-keep-size 10000")

        //Was the index created?
        assertEquals(0, buildIndexResult.statusCode)
        assert(File("${tempHvcfDir}/kmerIndex.txt").exists())
    }

    @Test
    fun testSourceFromOtherChr() {

        //set up temporary file names
        val tempTestDir = "${TestExtension.tempDir}kmerTest/"
        val tempHvcfDir = "${tempTestDir}hvcfDir/"
        val tempDBPathDir = "${TestExtension.testOutputFastaDir}dbPath/"

        //copy hvcf files to temp directory,
        // include the ref hvcf to test what happens when samples have no haplotype in some ref range
        File("${tempHvcfDir}LineB.h.vcf").delete()
        listOf(TestExtension.smallseqLineAHvcfFile, "${TestExtension.smallSeqInputDir}LineB_agc_command_test.h.vcf", TestExtension.smallseqRefHvcfFile)
            .forEach { hvcfFile ->
                val dst = File("$tempHvcfDir${File(hvcfFile).name}")
                if (!dst.exists()) {
                    File(hvcfFile).copyTo(dst)
                }
            }

        //create a HaplotypeGraph from the hvcf files
        val buildIndexResult = BuildKmerIndex().test("--db-path $tempDBPathDir --hvcf-dir $tempHvcfDir " +
                "--max-arg-length 150 --index-file ${tempHvcfDir}kmerIndexOther.txt --initial-keep-size 10000")

        //Was the index created?
        assertEquals(0, buildIndexResult.statusCode)
        assert(File("${tempHvcfDir}/kmerIndexOther.txt").exists())

        //test using lines D and E to test condition: sampleContigList is empty
        //  caused by half or more of the haplotypes in a range having their source on a different contig
        File(tempHvcfDir).listFiles().forEach { it.delete() }
        listOf("${TestExtension.smallSeqInputDir}LineD.h.vcf", "${TestExtension.smallSeqInputDir}LineE.h.vcf")
            .forEach { hvcfFile ->
                val dst = File("$tempHvcfDir${File(hvcfFile).name}")
                if (!dst.exists()) {
                    File(hvcfFile).copyTo(dst)
                }
            }

        //create a HaplotypeGraph from the hvcf files
        val buildIndexResult2 = BuildKmerIndex().test("--db-path $tempDBPathDir --hvcf-dir $tempHvcfDir " +
                "--max-arg-length 150 --index-file ${tempHvcfDir}kmerIndexOther.txt --initial-keep-size 10000")

        //Was the index created?
        assertEquals(0, buildIndexResult2.statusCode)
        assert(File("${tempHvcfDir}/kmerIndexOther.txt").exists())

        //clean up the tempHvcfDir
        File(tempHvcfDir).listFiles().forEach { it.delete() }
    }

    @Test
    fun testKmersDiagnostics() {
        //set up temporary file names
        val tempTestDir = "${TestExtension.tempDir}kmerTest/"
        val tempHvcfDir = "${tempTestDir}hvcfDir/"
        val tempDBPathDir = "${TestExtension.testOutputFastaDir}/dbPath"

        //copy hvcf files to temp directory,
        // include the ref hvcf to test what happens when samples have no haplotype in some ref range
        File("${tempHvcfDir}LineB.h.vcf").delete()
        listOf(TestExtension.smallseqLineAHvcfFile,"${TestExtension.smallSeqInputDir}LineB_shiftedToAdjacentRange.h.vcf", TestExtension.smallseqRefHvcfFile)
            .forEach { hvcfFile ->
                val dst = File("$tempHvcfDir${File(hvcfFile).name}")
                if (!dst.exists()) {
                    File(hvcfFile).copyTo(dst)
                }
            }

        //create a HaplotypeGraph from the hvcf files
        val buildIndexResult = BuildKmerIndex().test("--db-path $tempDBPathDir --hvcf-dir $tempHvcfDir --initial-keep-size 10000")

        //Was the index created?
        assertEquals(0, buildIndexResult.statusCode)
        assert(File("${tempHvcfDir}/kmerIndex.txt").exists())

        getBufferedReader("${tempHvcfDir}kmerIndexStatistics.txt").use { myReader ->
            var inputLine = myReader.readLine()
            var lineCount = 1
            while (inputLine != null) {
                //diagnostic report header is "contig	start	end	length	kmerCount	adjacentCount"
                //only lines 6 and 28 should have adjacentCount > 0
                if (lineCount == 2) assertEquals("1\t1\t1000\t1000\t634\t0", inputLine)
                if (lineCount == 8) assertEquals("1\t16501\t17500\t1000\t775\t0", inputLine)
                if (lineCount == 28) assertEquals("2\t16501\t17500\t1000\t771\t0", inputLine)
                if (lineCount == 35) assertEquals("2\t34001\t38500\t4500\t2773\t0", inputLine)
                inputLine = myReader.readLine()
                lineCount++
            }
        }

        //delete the diagnostic file then rerun with noDiagnostic flag
        File("${tempHvcfDir}kmerIndexStatistics.txt").delete()
        val noDiagnoticResult = BuildKmerIndex().test("--db-path $tempDBPathDir --hvcf-dir $tempHvcfDir -n --initial-keep-size 10000")
        assertEquals(0, noDiagnoticResult.statusCode)
        assert( !File("${tempHvcfDir}kmerIndexStatistics.txt").exists()) {"The diagnostic file was written but option is -n"}

        //delete the shifted B hvcf to avoid problems with other tests
        File("${tempHvcfDir}LineB_shiftedToAdjacentRange.h.vcf").delete()
    }

    //Ignore for now as we are requiring both db-path and hvcf-dir currently.
    @Disabled
    @Test
    fun testTiledb() {
        //Setting the tiledb path but not the hvcf should generate a not implemented error
        val tempAGCDir = "${TestExtension.testOutputFastaDir}/dbPath"

        //try to build a graph from a (non-existent) tiledb database
        val exception = assertThrows<NotImplementedError> {
            BuildKmerIndex().test("--db-path $tempAGCDir --db-path ${TestExtension.testTileDBURI}")
        }
        assertEquals("An operation is not implemented: TileDB VCF Reader Not implemented yet.  Please run with --hvcf-dir", exception.message)

    }

    @Test
    fun testGetRefRangeToKmerSetMap() {
        //fun getRefRangeToKmerSetMap(
        val buildKmerIndex = BuildKmerIndex()

        //create a kmerMapToHapIds
        val kmerMapToHapIds = Long2ObjectOpenHashMap<Set<String>>()

        //create a hapIdToRefRangeMap
        val hapIdToRefRangeMap = mutableMapOf<String, List<ReferenceRange>>()

        //add in some KmerToHapIdSets
        kmerMapToHapIds[1L] = setOf("hap1","hap2")
        kmerMapToHapIds[2L] = setOf("hap3","hap4")
        kmerMapToHapIds[3L] = setOf("hap1","hap3")
        kmerMapToHapIds[4L] = setOf("hap2","hap4")
        kmerMapToHapIds[5L] = setOf("hap1","hap2","hap3","hap4")
        kmerMapToHapIds[10L] = setOf("hap10","hap11","hap12")
        kmerMapToHapIds[11L] = setOf("hap10","hap11")
        kmerMapToHapIds[12L] = setOf("hap12")


        //create a hapId ToRangeMap
        hapIdToRefRangeMap["hap1"] = listOf(ReferenceRange("chr1",10,50))
        hapIdToRefRangeMap["hap2"] = listOf(ReferenceRange("chr1",10,50))
        hapIdToRefRangeMap["hap3"] = listOf(ReferenceRange("chr1",10,50))
        hapIdToRefRangeMap["hap4"] = listOf(ReferenceRange("chr1",10,50))
        hapIdToRefRangeMap["hap10"] = listOf(ReferenceRange("chr1",100,150))
        hapIdToRefRangeMap["hap11"] = listOf(ReferenceRange("chr1",100,150))
        hapIdToRefRangeMap["hap12"] = listOf(ReferenceRange("chr1",100,150))

        //setup the truth
        val truth = mapOf(
            ReferenceRange("chr1",10,50) to setOf(1L,2L,3L,4L,5L),
            ReferenceRange("chr1",100,150) to setOf(10L,11L,12L)
        )

        assertEquals(truth,buildKmerIndex.getRefRangeToKmerSetMap(kmerMapToHapIds,hapIdToRefRangeMap))

    }


    @Test
    fun testCreateHapIdToKmerMap() {

        val buildKmerIndex = BuildKmerIndex()

        //create a kmerMapToHapIds
        val kmerMapToHapIds = mutableMapOf<Long,Set<String>>()

        //add in some KmerToHapIdSets
        kmerMapToHapIds[1L] = setOf("hap1","hap2")
        kmerMapToHapIds[2L] = setOf("hap3","hap4")
        kmerMapToHapIds[3L] = setOf("hap1","hap3")
        kmerMapToHapIds[4L] = setOf("hap2","hap4")
        kmerMapToHapIds[5L] = setOf("hap1","hap2","hap3","hap4")
        kmerMapToHapIds[10L] = setOf("hap10","hap11","hap12")
        kmerMapToHapIds[11L] = setOf("hap10","hap11")
        kmerMapToHapIds[12L] = setOf("hap12")

        val truth = mapOf(setOf("hap1", "hap2") to listOf(1L),
            setOf("hap3", "hap4") to listOf(2L),
            setOf("hap1", "hap2", "hap3", "hap4") to listOf(5L),
            setOf("hap1", "hap3") to listOf(3L),
            setOf("hap2", "hap4") to listOf(4L),
            setOf("hap10", "hap11", "hap12") to listOf(10L),
            setOf("hap10", "hap11") to listOf(11L),
            setOf("hap12") to listOf(12L))

        assertEquals(truth,buildKmerIndex.createHapIdToKmerMap(kmerMapToHapIds))

    }

    @Test
    fun testBuildEncodedHapSetsAndHashOffsets() {
        val buildKmerIndex = BuildKmerIndex()
        //make a range
        val refRange = ReferenceRange("chr1",10,50)
        //There are 5 haplotypes in this refRange
        val numberOfHaplotypesInRange = 5
        //There are 3 hapSets in this refRange
        val numberOfHapSets = 3
        //make a refRangeToHapIndexMap
        val refRangeToHapIndexMap = mapOf(
            ReferenceRange("chr1",10,50) to mapOf("hap1" to 0,"hap2" to 1,"hap3" to 2,"hap4" to 3,"hap5" to 4)
        )
        //make a hapidKmerHashMap
        val hapidKmerHashMap = mapOf<Set<String>,List<Long>>(
            setOf("hap1","hap2") to listOf(1L),
            setOf("hap3","hap4") to listOf(2L),
            setOf("hap1","hap2","hap3","hap4") to listOf(5L),
        )

        //Create the truth which should have this type
        //Pair<BitSet, MutableList<Pair<Long, Int>>>
        val truthBitset = BitSet()
        truthBitset.set(0)
        truthBitset.set(1)
        truthBitset.set(7)
        truthBitset.set(8)
        truthBitset.set(10)
        truthBitset.set(11)
        truthBitset.set(12)
        truthBitset.set(13)

        //TODO Check this, I am not fully understanding the kmer hash long here.  The offsets make sense, but the kmer hashes need another look.
        val truth = Pair(truthBitset, mutableListOf(Pair(1L,0),Pair(2L,5),Pair(5L,10)))

        assertEquals(truth,buildKmerIndex.buildEncodedHapSetsAndHashOffsets(numberOfHaplotypesInRange,numberOfHapSets,refRangeToHapIndexMap,refRange,hapidKmerHashMap))

    }
    @Test
    fun testExportRefRangeKmerIndex() {
        //Very simple unit test to make sure that we are exporting what is expected.
        //create a refRange
        val refRange = ReferenceRange("chr1",10,50)

        //make an encoded HapBitSet
        val encodedHapSets = BitSet()
        encodedHapSets.set(1)
        encodedHapSets.set(2)
        encodedHapSets.set(3)

        //make a kmerHashOffsets
        val kmerHashOffsets = listOf(Pair(1L,1),Pair(2L,2),Pair(3L,3))

        getBufferedWriter(TestExtension.testKmerIndex).use { writer ->
            BuildKmerIndex().exportRefRangeKmerIndex(writer,refRange,encodedHapSets,kmerHashOffsets)
        }


        //read in the file and make sure it is correct
        val truth = listOf(">chr1:10-50","14","1@1,2@2,3@3")
        //compare truth with the file
        val inputFile = getBufferedReader(TestExtension.testKmerIndex).readLines()
        assertEquals(truth,inputFile)
    }

    @Test
    fun testCreateAgcCommandLists() {
        val graph = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFile,
            TestExtension.smallseqLineBHvcfFile, TestExtension.smallseqRefHvcfFile))
        //chr 1, lineB regions in the alt header are the same as the reference range
        val commandLists = BuildKmerIndex.rangeListsForAgcCommand(graph, graph.rangesByContig(), "1")
        println("---------sampleContigList---------")
        println(commandLists.sampleContigList)
        println("---------otherRegionsList---------")
        println(commandLists.otherRegionsList)
        assertTrue(commandLists.sampleContigList.contains("1@LineA"), "command list does not contain 1@LineA")
        assertTrue ( commandLists.sampleContigList.contains("1@LineB"), "command list does not contain 1@LineB")
        assertTrue ( commandLists.sampleContigList.contains("1@Ref"), "command list does not contain 1@Ref" )
        assertTrue(commandLists.otherRegionsList.isEmpty())

        val graph2 = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFile,
            "${TestExtension.smallSeqInputDir}LineB_kmer_index_test.h.vcf", TestExtension.smallseqRefHvcfFile))
        //all lineB regions for reference chr1 are changed to chr2 in LineB_kmer_index_test.h.vcf. Also, chr2 reference range regions are change to chr 1
        //as a result graph chr1 sequence should come from LineB chr2 (2@LineB)
        val commandLists2 = BuildKmerIndex.rangeListsForAgcCommand(graph2, graph2.rangesByContig(), "1")
        println("---------sampleContigList---------")
        println(commandLists2.sampleContigList)
        println("---------otherRegionsList---------")
        println(commandLists2.otherRegionsList)
        assertTrue(commandLists2.sampleContigList.contains("1@LineA"), "command list does not contain 1@LineA")
        assertTrue ( commandLists2.sampleContigList.contains("2@LineB"), "command list does not contain 2@LineB")
        assertTrue ( commandLists2.sampleContigList.contains("1@Ref"), "command list does not contain 1@Ref" )
        assertTrue(commandLists2.otherRegionsList.isEmpty())

        val graph3 = HaplotypeGraph(listOf(TestExtension.smallseqLineAHvcfFile,
            "${TestExtension.smallSeqInputDir}LineB_agc_command_test.h.vcf", TestExtension.smallseqRefHvcfFile))
        //LineB_kmer_index_test.h.vcf but with two regions in each chromosome set back to same chr as reference
        val commandLists3 = BuildKmerIndex.rangeListsForAgcCommand(graph3, graph3.rangesByContig(), "1")
        println("---------sampleContigList---------")
        println(commandLists3.sampleContigList)
        println("---------otherRegionsList---------")
        println(commandLists3.otherRegionsList)
        assertTrue(commandLists3.sampleContigList.contains("1@LineA"), "command list does not contain 1@LineA")
        assertTrue ( commandLists3.sampleContigList.contains("2@LineB"), "command list does not contain 2@LineB")
        assertTrue ( commandLists3.sampleContigList.contains("1@Ref"), "command list does not contain 1@Ref" )
        assertEquals(listOf("1@LineB:5501-6500", "1@LineB:34001-38500"), commandLists3.otherRegionsList, "other regions list not as expected for commandLists3")
    }

}