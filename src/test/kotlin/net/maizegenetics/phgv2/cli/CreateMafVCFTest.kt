package net.maizegenetics.phgv2.cli

import biokotlin.genome.AssemblyVariantInfo
import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.testing.test
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.utils.*
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.*

@ExtendWith(TestExtension::class)
class CreateMafVCFTest {
    companion object {
        @JvmStatic
        @AfterAll
        fun tearDown() {
            File(TestExtension.testVCFDir).deleteRecursively()
            File(TestExtension.testTileDBURI).deleteRecursively()
        }
    }

    @Test
    fun testCliktParams() {
        val createMAFVCF = CreateMafVcf()

        val resultMissingBed =
            createMAFVCF.test("--db-path ${TestExtension.testTileDBURI} --maf-dir ${TestExtension.testMafDir} --reference-file ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingBed.statusCode, 1)
        assertEquals(
            "Usage: create-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --bed\n", resultMissingBed.output
        )
        val resultMissingRef =
            createMAFVCF.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} --maf-dir ${TestExtension.testMafDir} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingRef.statusCode, 1)
        assertEquals(
            "Usage: create-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --reference-file\n", resultMissingRef.output
        )

        val resultMissingOutput =
            createMAFVCF.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} --maf-dir ${TestExtension.testMafDir} --reference-file ${TestExtension.testRefFasta}")
        assertEquals(resultMissingOutput.statusCode, 1)
        assertEquals(
            "Usage: create-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --output-dir\n", resultMissingOutput.output
        )

        val resultMissingMafDir =
            createMAFVCF.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} --reference-file ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")
        assertEquals(resultMissingMafDir.statusCode, 1)
        assertEquals(
            "Usage: create-maf-vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --maf-dir\n", resultMissingMafDir.output
        )

    }


    @Test
    fun testLoadingBEDFile() {
        val bedFile = "data/test/buildMAFVCF/B73_Test.bed"

        val ranges = loadRanges(bedFile)

        assertEquals(4, ranges.size)

        assertEquals(Position("chr1", 1), ranges[0].first)
        assertEquals(Position("chr1", 40), ranges[0].second)

        assertEquals(Position("chr7", 15), ranges[1].first)
        assertEquals(Position("chr7", 48), ranges[1].second)
        assertEquals(Position("chr7", 451), ranges[2].first)
        assertEquals(Position("chr7", 456), ranges[2].second)

        assertEquals(Position("chr10", 1), ranges[3].first)
        assertEquals(Position("chr10", 40), ranges[3].second)

    }

    @Test
    fun testBuildRefSeq() {
        val testRef = "data/test/buildMAFVCF/B73_Test.fa"

        val createMafVcf = CreateMafVcf()
        val refSeq = createMafVcf.buildRefGenomeSeq(testRef)

        assertEquals(3, refSeq.size)
        //check seqSizes
        assertEquals(461, refSeq["chr7"]!!.size())
        assertEquals(40, refSeq["chr1"]!!.size())
        assertEquals(40, refSeq["chr10"]!!.size())

        //check seqs
        val truthChr7Seq = "AAAAAAAAAAAAAAAGGGAATGTTAACCAAATGAATTGTCTCTTACGGTGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAGATGGGT"

        assertEquals(truthChr7Seq, refSeq["chr7"]!!.seq())
        assertEquals("GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT", refSeq["chr1"]!!.seq())
        assertEquals("GCAGCTGAAAACAGTCAATCTTACACACTTGGGGCCTACT", refSeq["chr10"]!!.seq())
    }

    @Test
    fun testSimpleBuildMAFVCF() {
        //Need to create the agc record before we run this:
        val fastaCreateFileNamesFile = "data/test/buildMAFVCF/fastaCreateFileNames.txt"
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/buildMAFVCF/B73_Test.fa"

        //Create the tileDB datasets - these are verified in AgcCompress()
        Initdb().createDataSets(TestExtension.testTileDBURI,"")

        val agcCompress = AgcCompress()
        // Create the initial compressed file
        val agcResult = agcCompress.test("--fasta-list ${fastaCreateFileNamesFile} --db-path ${dbPath} --reference-file ${refFasta}")
        println(agcResult.output)

        val createMAFVCF = CreateMafVcf()
        val result = createMAFVCF.test("--db-path ${dbPath} --bed data/test/buildMAFVCF/B73_Test.bed --reference-file ${refFasta} --maf-dir data/test/buildMAFVCF/mafs/ -o ${TestExtension.testVCFDir}")
        println(result.output)
        //compare the contents of the output gVCF files to the expected output
        compareTwoGVCFFiles("data/test/buildMAFVCF/truthGVCFs/B97_truth.g.vcf", "${TestExtension.testVCFDir}/B97_ASM_Test.g.vcf.gz")

        //test that metrics file was created
        assertTrue(File("${TestExtension.testVCFDir}/VCFMetrics.tsv").exists())

        //Now we need to compare the hVCF's sequence with the sequence coming from the MAF files to make sure things match correctly as well as the boundaries
        val outputHVCF = "${TestExtension.testVCFDir}/B97_ASM_Test.h.vcf.gz"
        val outputHVCFReader = VCFFileReader(File(outputHVCF),false)

        val outputHeader = outputHVCFReader.header.metaDataInInputOrder.filter { it.key == "ALT" }

        //pull out the hashes from the output header.
        val outputHashes = outputHeader.map { it.toString().split("<")[1].split(",").first().split("=").last() }.toSet()

        //Manually get the sequences out of the reference file
        val seqsManuallyExtractedFromMAF = setOf("AAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT",
            "AGGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATG",
            "TAAGGA"
            )

        val mafHashes = seqsManuallyExtractedFromMAF.map { getChecksumForString(it, "Md5") }.toSet()
        //compare the hashsets
        assertEquals(mafHashes,outputHashes)

        //Check the coordinates of the GVCF and make sure these match the hvcfs.
        //These were manually found by comparing BED and GVCF so they are right
        //Ref coord chr1	0	40 cooresponds to chr6 98 .. 142
        //chr7	14	48 corresponds to chr4 4245 .. 4281
        //chr7	450	456 corresponds to chr4 5247 .. 5252
        //chr10	0	40 corresponds to chr6 1098 .. 1142
        val truthCoords = setOf(Pair(Position("chr6", 98), Position("chr6", 142)),
            Pair(Position("chr4", 4245), Position("chr4", 4281)),
            Pair(Position("chr4", 5247), Position("chr4", 5252)),
            )

        val seqCoords = outputHeader.map {
            val tokens = it.toString().split("<")[1]
                .split(",")
                .map { token -> token.split("=") }.associate { token -> Pair(token[0], token[1]) }

            val region = tokens["Regions"]!!
            val regionSplit = region.split(":")
            val chr = regionSplit[0]
            val bounds = regionSplit[1].split("-")
            val start = bounds[0].toInt()
            val end = bounds[1].toInt()
            Pair(Position(chr, start), Position(chr, end))
        }.toSet()

        assertEquals(truthCoords, seqCoords)

        //Check the refCoords of the variants
        val bedCoords = setOf(Pair(Position("chr1", 1), Position("chr1", 40)),
            Pair(Position("chr7", 15), Position("chr7", 48)),
            Pair(Position("chr7", 451), Position("chr7", 456)),
            Pair(Position("chr10", 1), Position("chr10", 40)))
        val variants = outputHVCFReader.iterator().toList()
        val refCoords = variants.map { Pair(Position(it.contig,it.start),Position(it.contig, it.getAttributeAsInt("END",0))) }.toSet()

        assertEquals(bedCoords, refCoords)
    }

    @Test
    fun verifyRefMatchesASM() {
        val createMAFVCF = CreateMafVcf()
        val testFasta = "data/test/buildMAFVCF/B97_ASM_Test.fa"
        val testSeqs = createMAFVCF.buildRefGenomeSeq(testFasta)
        //check first record
        //NucSeq is zerobased inclusive inclusive, VCF is 1-based inclusive inclusive.
        val firstSeq = testSeqs["chr6"]!![97..141].seq()
        assertEquals("AAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT",firstSeq)


        //Looking for gvcfCoords [1098,1142] which should be [1097,1141] in the NucSeq
        val secondSeq = testSeqs["chr6"]!![1097..1141].seq()
        assertEquals("AAAAAGACAGCTGAAAATATCAATCTTACACACTTGGGGCCTACT",secondSeq)

        //Looking for gvcfCoords [4243,4283]
        val firstPartOfChr7_asm4 = testSeqs["chr4"]!![4242 .. 4282].seq()
        assertEquals("AAAGGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG", firstPartOfChr7_asm4) //Need to add an A at start here because MAF block begins with deletion so we need to start at position 12

        //Looking for gvcfCoords[5247,5257]
        val secondPartOfChr7_asm4 = testSeqs["chr4"]!![5246 .. 5256].seq()
        assertEquals("TAAGGATCCCT",secondPartOfChr7_asm4)

    }

    @Test
    fun testCreateOutputVariantMetadata() {
        val createMafVcf = CreateMafVcf()
        // bedRanges: List<Pair<Position, Position>>,
        //        refGenomeSequence: Map<String, NucSeq>,
        //        variantInfos: List<AssemblyVariantInfo>,
        //        sampleName: String

        val bedRanges = listOf(Pair(Position("chr1",1),Position("chr1",10)),
            Pair(Position("chr1",11),Position("chr1",30)),
            Pair(Position("chr1",31),Position("chr1",40)),
            Pair(Position("chr1",41),Position("chr1",60))
        )

        //make a ref genome sequence with chr1 of length 60
        val refGenomeSequence = mapOf("chr1" to NucSeq("A".repeat(60)))

        val sampleName = "test"
        //make simple variant infos
        //data class AssemblyVariantInfo(
        //    var chr: String, var startPos: Int, var endPos: Int, var genotype: String, var refAllele: String,
        //    var altAllele: String, var isVariant: Boolean, var alleleDepths: IntArray = intArrayOf(),
        //    var asmChrom: String = "", var asmStart: Int = -1, var asmEnd: Int = -1, var asmStrand: String = "",
        //    var isMissing: Boolean = false
        //)
        val variantInfos = listOf(
            AssemblyVariantInfo("chr1", 1,10,"A","A",".",false, intArrayOf(10,0),"asm1",100,109,"+",false),
            AssemblyVariantInfo("chr1", 11,15,"A","A",".",false, intArrayOf(5,0),"asm1",110,114,"+",false),
            AssemblyVariantInfo("chr1", 16,16,"G","A","G",true, intArrayOf(3,7),"asm1",115,115,"+",false),
            AssemblyVariantInfo("chr1", 17,35,"A","A",".",false, intArrayOf(10,0),"asm1",116,134,"+",false),
            AssemblyVariantInfo("chr1", 36,36,"TTT","A","TTT",true, intArrayOf(4,6),"asm1",135,137,"+",false),
            AssemblyVariantInfo("chr1", 37,50,"A","A",".",false, intArrayOf(15,0),"asm1",138,149,"+",false),
            AssemblyVariantInfo("chr1", 51,53,"G","AAA","G",true, intArrayOf(10,0),"asm1",150,152,"+",false),
            AssemblyVariantInfo("chr1", 54,60,"A","A",".",false, intArrayOf(7,0),"asm1",153,159,"+",false)
        )

        val outputMetadata = createMafVcf.createOutputVariantMetadata(
            bedRanges,
            refGenomeSequence,
            variantInfos,
            sampleName
        )
        //data class HVCFRecordMetadata(val sampleName: String, val refSeq : String = "", val asmSeq : String = "",
        //                              val refContig : String, val refStart: Int, val refEnd: Int,
        //                              val asmRegions: List<Pair<Position,Position>>)
        assertEquals(4, outputMetadata.size)
        //check first metadata
        val firstMeta = outputMetadata[0]
        assertEquals("test", firstMeta.sampleName)
        assertEquals("A".repeat(10), firstMeta.refSeq)
        assertEquals("", firstMeta.asmSeq) //We put a "" for asmSeq in this function
        assertEquals("chr1", firstMeta.refContig)
        assertEquals(1, firstMeta.refStart)
        assertEquals(10, firstMeta.refEnd)
        assertEquals(1, firstMeta.asmRegions.size)
        assertEquals(Position("asm1",100), firstMeta.asmRegions[0].first)
        assertEquals(Position("asm1",109), firstMeta.asmRegions[0].second)

        val secondMeta = outputMetadata[1]
        assertEquals("test", secondMeta.sampleName)
        assertEquals("A".repeat(20), secondMeta.refSeq)
        assertEquals("", secondMeta.asmSeq)
        assertEquals("chr1", secondMeta.refContig)
        assertEquals(11, secondMeta.refStart)
        assertEquals(30, secondMeta.refEnd)
        assertEquals(1, secondMeta.asmRegions.size) //This merges the SNP into the refRange as it is consecutive ASM position wise
        assertEquals(Position("asm1",110), secondMeta.asmRegions[0].first)
        assertEquals(Position("asm1",129), secondMeta.asmRegions[0].second)


        val thirdMeta = outputMetadata[2]
        assertEquals("test", thirdMeta.sampleName)
        assertEquals("A".repeat(10), thirdMeta.refSeq)
        assertEquals("", thirdMeta.asmSeq)
        assertEquals("chr1", thirdMeta.refContig)
        assertEquals(31, thirdMeta.refStart)
        assertEquals(40, thirdMeta.refEnd)
        assertEquals(1, thirdMeta.asmRegions.size)
        assertEquals(Position("asm1",130), thirdMeta.asmRegions[0].first)
        assertEquals(Position("asm1",141), thirdMeta.asmRegions[0].second)

        val fourthMeta = outputMetadata[3]
        assertEquals("test", fourthMeta.sampleName)
        assertEquals("A".repeat(20), fourthMeta.refSeq)
        assertEquals("", fourthMeta.asmSeq)
        assertEquals("chr1", fourthMeta.refContig)
        assertEquals(41, fourthMeta.refStart)
        assertEquals(60, fourthMeta.refEnd)
        assertEquals(1, fourthMeta.asmRegions.size)
        assertEquals(Position("asm1",142), fourthMeta.asmRegions[0].first)
        assertEquals(Position("asm1",159), fourthMeta.asmRegions[0].second)


        //Test some more difficult indels
        //Overlapping the bed range start and end and a BED fully contained within a deletion
        val bedRanges2 = listOf(Pair(Position("chr1",5),Position("chr1",25)),
            Pair(Position("chr1",26),Position("chr1",40)),
            Pair(Position("chr1",41),Position("chr1",45)),
            Pair(Position("chr1",46),Position("chr1",60)))

        val variantInfos2 = listOf(
            AssemblyVariantInfo("chr1", 1,10,"AAAAAAAAAA","AAAAAAAAAA","A",true, intArrayOf(10,0),"asm1",200,200,"+",false),
            AssemblyVariantInfo("chr1",24,27,"A","AAAA","A",true, intArrayOf(8,2),"asm1",210,210,"+",false),
            AssemblyVariantInfo("chr1", 38,47,"A","AAAAAAAAAA","A",true, intArrayOf(5,5),"asm1",220,220,"+",false),

//            AssemblyVariantInfo("chr2", 11,20,"A","A",".",false, intArrayOf(5,0),"asm2",210,219,"+",false),
//            AssemblyVariantInfo("chr2", 21,21,"G","A","G",true, intArrayOf(3,7),"asm2",220,220,"+",false),
//            AssemblyVariantInfo("chr2", 22,40,"A","A",".",false, intArrayOf(10,0),"asm2",221,239,"+",false),
//            AssemblyVariantInfo("chr2", 41,60,"TTT","A","TTT",true, intArrayOf(4,6),"asm2",240,242,"+",false),
//            AssemblyVariantInfo("chr2", 61,70,"A","A",".",false, intArrayOf(15,0),"asm2",243,252,"+",false),
//            AssemblyVariantInfo("chr2", 71,80,"G","AAA","G",true, intArrayOf(10,0),"asm2",253,262,"+",false),
//            AssemblyVariantInfo("chr2", 81,90,"A","A",".",false, intArrayOf(7,0),"asm2",263,272,"+",false)
        )

        val outputMetadata2 = createMafVcf.createOutputVariantMetadata(
            bedRanges2,
            refGenomeSequence,
            variantInfos2,
            sampleName
        )
        assertEquals(2, outputMetadata2.size)
        //Even though this should span 2 ref ranges, the left alignment means that it is assigned to 5-25 not 26-40 and definitely not both
        val metaData2_0 = outputMetadata2[0]
        assertEquals("test", metaData2_0.sampleName)
        assertEquals("A".repeat(21), metaData2_0.refSeq) //5 to 25 inclusive  full ref seq from BED
        assertEquals("", metaData2_0.asmSeq)
        assertEquals("chr1", metaData2_0.refContig)
        assertEquals(5, metaData2_0.refStart)
        assertEquals(25, metaData2_0.refEnd)
        assertEquals(1, metaData2_0.asmRegions.size)
        assertEquals(Position("asm1",210), metaData2_0.asmRegions[0].first)
        assertEquals(Position("asm1",210), metaData2_0.asmRegions[0].second)

        //This one should be assigned to 26-40 even though it fully overlaps 41-45 and partially overlaps 46-60 because of left alignment
        val metaData2_1 = outputMetadata2[1]
        assertEquals("test", metaData2_1.sampleName)
        assertEquals("A".repeat(15), metaData2_1.refSeq) //
        assertEquals("", metaData2_1.asmSeq)
        assertEquals("chr1", metaData2_1.refContig)
        assertEquals(26, metaData2_1.refStart)
        assertEquals(40, metaData2_1.refEnd)
        assertEquals(1, metaData2_1.asmRegions.size)
        assertEquals(Position("asm1",220), metaData2_1.asmRegions[0].first)
        assertEquals(Position("asm1",220), metaData2_1.asmRegions[0].second)

    }

    /**
     * Function to compare the output gVCF file with the expected gVCF
     * It compares the alleles, depths, genotypes and the ASM metadata.
     */
    fun compareTwoGVCFFiles(truthGVCFFile : String, generatedFile: String) {
        //Load in the output GVCF  and the truth GVCF and verify that the output is correct
        val truthVariantIterator = VCFFileReader(File(truthGVCFFile),false).iterator()
        val truthVariants = mutableListOf<VariantContext>()
        while(truthVariantIterator.hasNext()) {
            truthVariants.add(truthVariantIterator.next())
        }
        val truthMap = truthVariants.associateBy { Position(it.contig, it.start) }

        val outputVariantIterator = VCFFileReader(File(generatedFile), false).iterator()
        val outputVariants = mutableListOf<VariantContext>()
        while(outputVariantIterator.hasNext()) {
            outputVariants.add(outputVariantIterator.next())
        }

        assertEquals(truthVariants.size, outputVariants.size,"Number of Variants does not match:")

        for(variant in outputVariants) {
            if(!truthMap.containsKey(Position(variant.contig, variant.start))) {
                fail("No matching variant found: ${variant.contig}:${variant.start}")
            }
            val matchingTruth = truthMap[Position(variant.contig, variant.start)]!!

            //Check END
            assertEquals(matchingTruth.end, variant.end,
                "End position does not match: outputFile: ${variant.contig}:${variant.start}-${variant.end}, " +
                        "Truth: ${matchingTruth.contig}:${matchingTruth.start}-${matchingTruth.end}")

            //Check alleles
            for(element in matchingTruth.alleles.zip(variant.alleles)) {
                assertEquals(element.first, element.second, "Alleles do not match: ${variant.contig}:${variant.start}-${variant.end}")
            }
            //Check GT
            assertEquals( matchingTruth.getGenotype(0).genotypeString,variant.getGenotype(0).genotypeString,
                "GT Fields do not match: ${variant.contig}:${variant.start}-${variant.end}")
            //Check AD
            for(element in matchingTruth.getGenotype(0).ad.zip(variant.getGenotype(0).ad)) {
                assertEquals(element.first, element.second, "AD fields do not match: ${variant.contig}:${variant.start}-${variant.end}")
            }
            //Check ASM Contig
            assertEquals(matchingTruth.getAttribute("ASM_Chr"), variant.getAttribute("ASM_Chr"),
                "ASM_Contig does not match: ${variant.contig}:${variant.start}-${variant.end}")
            //Check ASM Start
            assertEquals(matchingTruth.getAttribute("ASM_Start"), variant.getAttribute("ASM_Start"),
                "ASM_Start does not match: ${variant.contig}:${variant.start}-${variant.end}")
            //Check ASM END
            assertEquals(matchingTruth.getAttribute("ASM_End"), variant.getAttribute("ASM_End"),
                "ASM_End does not match: ${variant.contig}:${variant.start}-${variant.end}")
            //Check ASM Strand
            assertEquals(matchingTruth.getAttribute("ASM_Strand"), variant.getAttribute("ASM_Strand"),
                "ASM_Strand does not match: ${variant.contig}:${variant.start}-${variant.end}")

        }
    }

}