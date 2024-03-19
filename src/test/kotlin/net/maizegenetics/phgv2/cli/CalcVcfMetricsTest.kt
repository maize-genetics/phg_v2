package net.maizegenetics.phgv2.cli

import htsjdk.variant.vcf.VCFHeader
import net.maizegenetics.phgv2.utils.getBufferedWriter
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.junit.jupiter.api.*
import java.io.File
import kotlin.test.assertEquals

class CalcVcfMetricsTest {

    @Test
    fun testSingleGVCFFromFile() {
        val inFile = "data/test/buildMAFVCF/truthGVCFs/B97_truth.g.vcf"

        val stats = CalcVcfMetrics().getGVCFStats(File(inFile))!!

        assertEquals(stats.size, 4)

        val chr1Stats = stats.filter{it.chrom == "chr1"}[0]
        val chr7Stats = stats.filter{it.chrom == "chr7"}[0]
        val chr10Stats = stats.filter{it.chrom == "chr10"}[0]
        val allStats = stats.filter{it.chrom == "ALL"}[0]

        // chr1
        // 29 ref bp, 11 SNPs, 1 ins (5 bp) 0 del, 0 Ns
        // length 100: %ref is 0.29, %mapped is 0.4
        assertEquals(0.29, chr1Stats.percentIdentityWithRef)
        assertEquals(0.4, chr1Stats.percentMappedToRef)
        assertEquals(11, chr1Stats.numSNPs)
        assertEquals(0, chr1Stats.numDel)
        assertEquals(1, chr1Stats.numIns)
        assertEquals(0, chr1Stats.numNs)
        assertEquals(0, chr1Stats.numBasesDeleted)
        assertEquals(5, chr1Stats.numBasesInserted)

        // chr7
        // 37 ref bp, 12 SNPs, 1 ins (3 bp) 1 del (1 bp), 0 Ns
        // length 461: %ref is 0.08, %mapped is 0.108
        assertEquals(37.0/461, chr7Stats.percentIdentityWithRef)
        assertEquals((37.0 + 12.0 + 1.0) / 461, chr7Stats.percentMappedToRef)
        assertEquals(12, chr7Stats.numSNPs)
        assertEquals(1, chr7Stats.numDel)
        assertEquals(1, chr7Stats.numIns)
        assertEquals(0, chr7Stats.numNs)
        assertEquals(1, chr7Stats.numBasesDeleted)
        assertEquals(3, chr7Stats.numBasesInserted)

        // chr10
        //  ref bp, 11 SNPs, 1 ins (5 bp) 0 del, 0 Ns
        // length 100: %ref is 0.29, %mapped is 0.4
        assertEquals(0.29, chr10Stats.percentIdentityWithRef)
        assertEquals(0.4, chr10Stats.percentMappedToRef)
        assertEquals(11, chr10Stats.numSNPs)
        assertEquals(0, chr10Stats.numDel)
        assertEquals(1, chr10Stats.numIns)
        assertEquals(0, chr10Stats.numNs)
        assertEquals(0, chr10Stats.numBasesDeleted)
        assertEquals(5, chr10Stats.numBasesInserted)

        // ALL
        // 95 ref bp, 34 SNPs, 1 ins (5 bp) 0 del, 0 Ns
        // length 661: %ref is 0.14, %mapped is 0.4
        assertEquals(95.0/661, allStats.percentIdentityWithRef)
        assertEquals((95.0 + 34.0 + 1)/661, allStats.percentMappedToRef)
        assertEquals(34, allStats.numSNPs)
        assertEquals(1, allStats.numDel)
        assertEquals(3, allStats.numIns)
        assertEquals(0, allStats.numNs)
        assertEquals(1, allStats.numBasesDeleted)
        assertEquals(13, allStats.numBasesInserted)


    }

    @Test
    fun testHVCFStats() {

        val rangeMap = mapOf(
            Pair(RefRangeInfo("a", "chr1", 1, 10), HVCFStatsByRange("a", 10, true)),
            Pair(RefRangeInfo("b", "chr1", 23, 66), HVCFStatsByRange("a", 23, false)),
            Pair(RefRangeInfo("c", "chr1", 524, 568), HVCFStatsByRange("a", 18, true)),
            Pair(RefRangeInfo("d", "chr1", 3658, 3685), HVCFStatsByRange("a", 879, false)),
            Pair(RefRangeInfo("e", "chr2", 12, 16), HVCFStatsByRange("a", 654, false)),
            Pair(RefRangeInfo("f", "chr2", 154, 521), HVCFStatsByRange("a", 136, false)),
            Pair(RefRangeInfo("g", "chr2", 542, 589), HVCFStatsByRange("a", 253, true)),
            Pair(RefRangeInfo("h", "chr2", 4678, 5809), HVCFStatsByRange("a", 534, false)),
        )

        val stats = CalcVcfMetrics().getHVCFStats("sample1", rangeMap)

        stats.forEach{println(it)}

        assertEquals(stats.size, 3)

        assertEquals(stats[0].refRangesWithHap, 8)
        assertEquals(stats[0].hapsIdenticalToRef, 3)

        assertEquals(stats[1].refRangesWithHap, 4)
        assertEquals(stats[1].hapsIdenticalToRef, 2)

        assertEquals(stats[2].refRangesWithHap, 4)
        assertEquals(stats[2].hapsIdenticalToRef, 1)

    }


    @Test
    fun testNCounting() {
        val metrics = CalcVcfMetrics().getGVCFStats(File("$testingDir/$NGvcfFile"))

        val chr1metrics = metrics!![0]

        assertEquals(61, chr1metrics.numNs)
    }

    @Test
    fun testWriteOutput() {
        val outGVCFMetrics = "$testingDir/gvcf_metrics.tsv"

        CalcVcfMetrics().calculateVcfMetrics(gvcfDir, outGVCFMetrics)
        val gvcfLines = File(outGVCFMetrics).readLines()

        assertEquals(getTruthTSV(), gvcfLines)
    }

    companion object {
        val testingDir = System.getProperty("user.home") + "/temp/VCFMetricsTests/"
        val gvcfDir = testingDir + "/gvcfs/"
        val gvcfFile = "sampleName.gvcf"
        val refGvcfFile = "Ref.g.vcf"
        val NGvcfFile = "n.g.vcf"

        @JvmStatic
        @BeforeAll
        fun setup() {
            setupDebugLogging()

            File(testingDir).deleteRecursively()
            File(testingDir).mkdirs()
            File(gvcfDir).mkdirs()

            createDummyGVCF("$gvcfDir/$gvcfFile")
            createRefGVCF("$gvcfDir/$refGvcfFile")
            createNGVCF("$testingDir/$NGvcfFile")
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
        }

        /**
         * Helper function creates a dummy reference gvcf. Formatted correctly, but with nonsense data.
         */
        private fun createRefGVCF(outputFile: String) {

            getBufferedWriter(outputFile).use { output ->
                output.write(
                    "##fileformat=VCFv4.2\n" +
                            "##FORMAT=<ID=AD,Number=3,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n" +
                            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n" +
                            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" +
                            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                            "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n" +
                            "##INFO=<ID=AF,Number=3,Type=Integer,Description=\"Allele Frequency\">\n" +
                            "##INFO=<ID=ASM_Chr,Number=1,Type=String,Description=\"Assembly chromosome\">\n" +
                            "##INFO=<ID=ASM_End,Number=1,Type=Integer,Description=\"Assembly end position\">\n" +
                            "##INFO=<ID=ASM_Start,Number=1,Type=Integer,Description=\"Assembly start position\">\n" +
                            "##INFO=<ID=ASM_Strand,Number=1,Type=String,Description=\"Assembly strand\">\n" +
                            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n" +
                            "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n" +
                            "##contig=<ID=1,length=40308>\n" +
                            "##contig=<ID=2,length=32822>\n" +
                            "##contig=<ID=3,length=10000>\n" +
                            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tRef\n" +
                            "1\t1\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5000;ASM_Start=1;ASM_Strand=+;END=5000\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "1\t5001\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5000;ASM_Start=1;ASM_Strand=+;END=10000\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "1\t10001\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5000;ASM_Start=1;ASM_Strand=+;END=15000\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "1\t15001\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5000;ASM_Start=1;ASM_Strand=+;END=20000\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "1\t20001\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5000;ASM_Start=1;ASM_Strand=+;END=25000\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "1\t25001\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5000;ASM_Start=1;ASM_Strand=+;END=30000\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "1\t30001\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5000;ASM_Start=1;ASM_Strand=+;END=35000\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "1\t35001\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5000;ASM_Start=1;ASM_Strand=+;END=40000\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "1\t40001\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=308;ASM_Start=1;ASM_Strand=+;END=40308\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "2\t1\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5000;ASM_Start=1;ASM_Strand=+;END=5000\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "2\t5001\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5000;ASM_Start=1;ASM_Strand=+;END=10000\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "2\t10001\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5000;ASM_Start=1;ASM_Strand=+;END=15000\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "2\t15001\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5000;ASM_Start=1;ASM_Strand=+;END=20000\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "2\t20001\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5000;ASM_Start=1;ASM_Strand=+;END=25000\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "2\t25001\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5000;ASM_Start=1;ASM_Strand=+;END=30000\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "2\t30001\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2822;ASM_Start=1;ASM_Strand=+;END=32822\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"
                )
            }

            val builder = ProcessBuilder("bgzip", "$outputFile")
            val process = builder.start()
            var error = process.waitFor()

            val builder2 = ProcessBuilder("tabix", "$outputFile.gz")
            val process2 = builder2.start()
            error += process2.waitFor()

            if (error != 0) {
                println("Something went wrong creating the reference gvcf. Check that bgzip and tabix are installed on this machine.")
            }
        }


        /**
         * Helper function creates a dummy gvcf. Formatted correctly, but with nonsense data
         * Contains two chromosomes, snps, indels, and unaligned segments
         * The ASM start and end positions are nonsense
         */
        private fun createDummyGVCF(outputFile: String) {
            getBufferedWriter(outputFile).use { output ->
                output.write("##fileformat=VCFv4.2\n" +
                        "##FORMAT=<ID=AD,Number=3,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n" +
                        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n" +
                        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" +
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                        "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n" +
                        "##INFO=<ID=AF,Number=3,Type=Integer,Description=\"Allele Frequency\">\n" +
                        "##INFO=<ID=ASM_Chr,Number=1,Type=String,Description=\"Assembly chromosome\">\n" +
                        "##INFO=<ID=ASM_End,Number=1,Type=Integer,Description=\"Assembly end position\">\n" +
                        "##INFO=<ID=ASM_Start,Number=1,Type=Integer,Description=\"Assembly start position\">\n" +
                        "##INFO=<ID=ASM_Strand,Number=1,Type=String,Description=\"Assembly strand\">\n" +
                        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n" +
                        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n" +
                        "##contig=<ID=1,length=40308>\n" +
                        "##contig=<ID=2,length=32822>\n" +
                        "##contig=<ID=3,length=10000>\n" +
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleName\n" +
                        "1\t1\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2864;ASM_Start=2857;ASM_Strand=+;END=8\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t9\t.\tA\tTATTCTGCGGAGAGCCTTAGTGAGTATAATGAAGTCACATAACTGTTTTGCTACCCTTCTCTCGGACTCCCCCGGATGAT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7181;ASM_Start=7103;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t10\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8036;ASM_Start=7987;ASM_Strand=+;END=59\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t60\t.\tGGCCCGTGTGAGATACTCTTATAGCGCAGGCTGGACAATACGGAGC\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8829;ASM_Start=8829;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t106\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1166;ASM_Start=697;ASM_Strand=+;END=575\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t576\t.\tG\tACCGATCAATCCTAATGTGTCTACGTGCCAGCCTGGACGTCTTAATAGTCGTATTAAGCGACGTAG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=702;ASM_Start=638;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t577\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2375;ASM_Start=2366;ASM_Strand=+;END=586\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t587\t.\tT\tCTCGCCGCCGAAGAGCTACCGTCTCGCCGGTATT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2792;ASM_Start=2760;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t588\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9649;ASM_Start=9266;ASM_Strand=+;END=971\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t972\t.\tCTTTCTGCCGTCTCAGCTGGAAGAAGAATCTCGCAGAAAGAA\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2654;ASM_Start=2654;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t1014\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8709;ASM_Start=8476;ASM_Strand=+;END=1247\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t1248\t.\tC\tGCTTTTTTGGCGCTCCGACATGAACAACTAACGAAACCCCGCTCGAGGCGCTGTTAGTAGGCCAGGCGTCCGTCTTGCTTCCGATCAGGAAAAA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5225;ASM_Start=5133;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t1249\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=10570;ASM_Start=9751;ASM_Strand=+;END=2068\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t2069\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=188;ASM_Start=188;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t2070\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3212;ASM_Start=2849;ASM_Strand=+;END=2433\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t2434\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4300;ASM_Start=4300;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t2435\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2070;ASM_Start=1655;ASM_Strand=+;END=2850\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t2851\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2706;ASM_Start=2706;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t2852\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5374;ASM_Start=4438;ASM_Strand=+;END=3788\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t3789\t.\tT\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2588;ASM_Start=2588;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t3790\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6652;ASM_Start=6501;ASM_Strand=+;END=3941\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t4283\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2656;ASM_Start=2521;ASM_Strand=+;END=4418\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t4419\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3201;ASM_Start=3201;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t4420\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9408;ASM_Start=8690;ASM_Strand=+;END=5138\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t5139\t.\tT\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1584;ASM_Start=1584;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t5140\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3636;ASM_Start=2814;ASM_Strand=+;END=5962\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t5963\t.\tC\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4713;ASM_Start=4713;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t5964\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3686;ASM_Start=3319;ASM_Strand=+;END=6331\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t6332\t.\tC\tTGTGTGCAGGAGCCTGGTAGACGAGGCAATGTCGGTGAACAAGACCCTTGATTAGGCATT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7640;ASM_Start=7582;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t6333\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8891;ASM_Start=8732;ASM_Strand=+;END=6492\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t6493\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3195;ASM_Start=3195;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t6494\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=856;ASM_Start=777;ASM_Strand=+;END=6573\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t6574\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7018;ASM_Start=7018;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t6575\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8262;ASM_Start=7884;ASM_Strand=+;END=6953\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t6954\t.\tT\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8536;ASM_Start=8536;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t6955\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9468;ASM_Start=8647;ASM_Strand=+;END=7776\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t7777\t.\tCCTCCCGCAAAGTCGCGAGCGTTAGGTGGCCCGACGCCACTTCGGTTGTTAG\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=684;ASM_Start=684;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t7829\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=10116;ASM_Start=9541;ASM_Strand=+;END=8404\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t8405\t.\tA\tCACTCTAAGGTGCGTTCCCCACTATGGCGAATTCTACTGTCCA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4431;ASM_Start=4390;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t8406\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2996;ASM_Start=2847;ASM_Strand=+;END=8555\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t8556\t.\tAATTTCGAATAGTCAATGCAGCCGAATTAGGG\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5409;ASM_Start=5409;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t8588\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5036;ASM_Start=5030;ASM_Strand=+;END=8594\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t9582\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1644;ASM_Start=801;ASM_Strand=+;END=10425\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t10426\t.\tC\tGTTTGGCGGAAAAGGTGACGGGGCTAACAGAGAGACAAGGAAAGGAACTACCTAAA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8844;ASM_Start=8790;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t10427\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4003;ASM_Start=3448;ASM_Strand=+;END=10982\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t10983\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5517;ASM_Start=5517;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t10984\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=10200;ASM_Start=9653;ASM_Strand=+;END=11531\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t11532\t.\tT\tGGATCTAATCCTTCCGATTATCTTGGTCCGTGCTTTGAAAGTATGTCCCGCCTAATATATGACACTTACTCGCCTT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8835;ASM_Start=8761;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t11533\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7480;ASM_Start=6763;ASM_Strand=+;END=12250\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t12251\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2469;ASM_Start=2469;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t12252\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7925;ASM_Start=7686;ASM_Strand=+;END=12491\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t12492\t.\tT\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8829;ASM_Start=8829;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t12493\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2996;ASM_Start=2854;ASM_Strand=+;END=12635\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t12636\t.\tA\tAAAACTCT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3543;ASM_Start=3537;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t12637\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1752;ASM_Start=963;ASM_Strand=+;END=13426\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t13427\t.\tGGTTACTAATATCCGCCCCTTTACAAAGCACCTGGATTTACTGAGCAGTT\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5770;ASM_Start=5770;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t13477\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=734;ASM_Start=698;ASM_Strand=+;END=13513\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t13514\t.\tCCTCGCCGCTGCCATGGGATTCGCCCGGGGAGCTGGTCCCAGCTTCTTGTTGATACTCAAGCCGTAATGGTACGTACTCCA\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6822;ASM_Start=6822;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t13595\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6172;ASM_Start=5250;ASM_Strand=+;END=14517\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t14518\t.\tGGACAACAGCAGGACGTAAATGAGTTGTTCTTTCGTAGGTGGCTGTAACCGTCCTGGCCAGAGACACAATCAACACCACAATCTTA\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9078;ASM_Start=9078;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t14604\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8716;ASM_Start=8071;ASM_Strand=+;END=15249\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t15250\t.\tC\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1862;ASM_Start=1862;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t15251\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8123;ASM_Start=7344;ASM_Strand=+;END=16030\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t16031\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4628;ASM_Start=4628;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t16032\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3400;ASM_Start=2766;ASM_Strand=+;END=16666\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t16667\t.\tT\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2513;ASM_Start=2513;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t16668\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5627;ASM_Start=4922;ASM_Strand=+;END=17373\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t17374\t.\tA\tCGTGATCGTGGTATGCCTTCAACGTCAGCGCCATAAGGAGCACCATTCGGCGCGTGGCA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2798;ASM_Start=2741;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t17375\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5258;ASM_Start=4951;ASM_Strand=+;END=17682\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t17683\t.\tC\tAACTGCTCATTTGACTCTTAGGTTGTATTCAAGCCGATATTTAAACA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1633;ASM_Start=1588;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t17684\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8839;ASM_Start=8586;ASM_Strand=+;END=17937\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t17938\t.\tTAGACAACTGTATGGCGGTTAAGAGCGTGCAAACCTCAGATGGCACGAGGGCCAGGGTACCGGCACACTT\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5497;ASM_Start=5497;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t18008\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3842;ASM_Start=3646;ASM_Strand=+;END=18204\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t18205\t.\tG\tATTTGCCTGATTGAGTGGTCGCGCCTACGGGTTTCTTGGCGGAATGCAGCGGAACATTACCTTAAGT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1177;ASM_Start=1112;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t18206\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1873;ASM_Start=1552;ASM_Strand=+;END=18527\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t18528\t.\tT\tGCCTCGACCCTGT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8127;ASM_Start=8116;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t18529\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9252;ASM_Start=8468;ASM_Strand=+;END=19313\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t19314\t.\tT\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=913;ASM_Start=913;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t19315\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6523;ASM_Start=6148;ASM_Strand=+;END=19690\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t19691\t.\tG\tGTGTTACGAGCCGCACT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=338;ASM_Start=323;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t19692\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1333;ASM_Start=905;ASM_Strand=+;END=20120\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t20121\t.\tG\tAGCGGCATT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1896;ASM_Start=1889;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t20122\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2066;ASM_Start=1469;ASM_Strand=+;END=20719\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t20720\t.\tT\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6198;ASM_Start=6198;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t20721\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6225;ASM_Start=5418;ASM_Strand=+;END=21528\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t21529\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3233;ASM_Start=3233;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t21530\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4886;ASM_Start=4314;ASM_Strand=+;END=22102\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t22103\t.\tC\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9127;ASM_Start=9127;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t22104\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8475;ASM_Start=8312;ASM_Strand=+;END=22267\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t22268\t.\tT\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6000;ASM_Start=6000;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t22269\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3270;ASM_Start=2410;ASM_Strand=+;END=23129\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t23130\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7864;ASM_Start=7864;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t23131\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1080;ASM_Start=107;ASM_Strand=+;END=24104\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t24105\t.\tGGCATTACCGCTGCTCCTAATATCAAAGTCCTGCAACTGTTAGCGGACTAGCAGACTCGGGTTGAACTCGGTTTTCTCCGTTCTGG\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=687;ASM_Start=687;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t24191\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7430;ASM_Start=7158;ASM_Strand=+;END=24463\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t24840\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9843;ASM_Start=9092;ASM_Strand=+;END=25591\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t25592\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2445;ASM_Start=2445;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t25593\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2922;ASM_Start=2250;ASM_Strand=+;END=26265\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t26266\t.\tAACGGGTGAGGACCCTGATATCCATTCGATTGCCACCGATCAGTCGTC\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7342;ASM_Start=7342;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t26314\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3091;ASM_Start=2984;ASM_Strand=+;END=26421\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t27037\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5131;ASM_Start=4971;ASM_Strand=+;END=27197\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t27198\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5066;ASM_Start=5066;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t27199\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=10216;ASM_Start=9581;ASM_Strand=+;END=27834\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t27835\t.\tG\tATCGACAAGATAGCCCTGCTCGTAACGGCGAGAGGACGTGTAGTGGATTTTTGATGACTAACA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8925;ASM_Start=8864;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t27836\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4678;ASM_Start=4251;ASM_Strand=+;END=28263\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t28264\t.\tG\tATCGGAGGCAGGTGAATCGATTCGGTTCTCTGAAACCTTCTTAGCCATCCAGGATCCCTGTAGCTAACGAGTGCTGCCCGCATTAGCACATCCATA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5290;ASM_Start=5196;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t28265\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1420;ASM_Start=1025;ASM_Strand=+;END=28660\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t28661\t.\tA\tGATCA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3097;ASM_Start=3094;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t28662\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=10022;ASM_Start=9491;ASM_Strand=+;END=29193\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t29194\t.\tT\tAGTTGCACCAATGTACGGAT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3107;ASM_Start=3089;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t29195\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9182;ASM_Start=9108;ASM_Strand=+;END=29269\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t29270\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3657;ASM_Start=3657;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t29271\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1370;ASM_Start=1325;ASM_Strand=+;END=29316\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t29317\t.\tTCCTCCTCGCGAAGTTGAACCACCGTATCTCCTACCCCCAGGAAGGTCGGCTAGATGTGCTGACAACTA\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9266;ASM_Start=9266;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t29386\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4078;ASM_Start=3917;ASM_Strand=+;END=29547\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t29548\t.\tTTGCAACACGAAAGGGCTATCATATGTCTCTGCTCGAGTGCA\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6311;ASM_Start=6311;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t29590\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=10362;ASM_Start=9528;ASM_Strand=+;END=30424\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t30425\t.\tAGGCTCTCA\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8477;ASM_Start=8477;ASM_Strand=+\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t30434\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4241;ASM_Start=3777;ASM_Strand=+;END=30898\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t30899\t.\tGAGTAAGATCCAAAGTAGACTAACCCACAACAAGTCAAAGCATGATGTCATCCGTGTTAACTCGGAAGAAGGTTGGAGTCTTCG\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=852;ASM_Start=852;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t30983\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4805;ASM_Start=3911;ASM_Strand=+;END=31877\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t31878\t.\tGTATTCCTTGAGCGTGCGCCGTCCC\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7425;ASM_Start=7425;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t31903\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7923;ASM_Start=7327;ASM_Strand=+;END=32499\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t32500\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8856;ASM_Start=8856;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t32501\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2209;ASM_Start=2105;ASM_Strand=+;END=32605\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t32606\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9943;ASM_Start=9943;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t32607\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3513;ASM_Start=3116;ASM_Strand=+;END=33004\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t33005\t.\tTCACACTATTCGGCGGGCCGAAAGCTTTACCCGCACCATCCTATACTCCTAAAGGTTTCCGAAGCCAACAGATATAAACGATATTACA\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4399;ASM_Start=4399;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t33093\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1722;ASM_Start=812;ASM_Strand=+;END=34003\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t34004\t.\tCAACCGCAGGCGAATGAAGACGACCAGGACCATGCAATTACAACATTGGCGGCAGCCCTAAGTAGTTAGAATTCTATCCTTCGCCTGAAAGACCGCCTGC\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3781;ASM_Start=3781;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t34104\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5960;ASM_Start=5828;ASM_Strand=+;END=34236\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t34237\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8297;ASM_Start=8297;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t34238\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8983;ASM_Start=8102;ASM_Strand=+;END=35119\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t35120\t.\tA\tTATAGCACTGGTCTGGCATCAGCGGTTGAGGCTCTACACGGTTAATTCTTCCGATTCGTAATTACCCGCGAGTGTTAGGGGCCAGTTGTCCGCAATTCG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2116;ASM_Start=2019;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t35121\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2876;ASM_Start=2792;ASM_Strand=+;END=35205\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t35206\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9002;ASM_Start=9002;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t35207\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8816;ASM_Start=8756;ASM_Strand=+;END=35267\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t35268\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5159;ASM_Start=5159;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t35269\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9708;ASM_Start=9044;ASM_Strand=+;END=35933\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t35934\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1371;ASM_Start=1371;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t35935\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1895;ASM_Start=1737;ASM_Strand=+;END=36093\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t36094\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2386;ASM_Start=2386;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t36095\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4777;ASM_Start=3829;ASM_Strand=+;END=37043\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t37044\t.\tT\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4497;ASM_Start=4497;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t37045\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=10111;ASM_Start=9363;ASM_Strand=+;END=37793\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t37794\t.\tG\tTATCAGAGTCTATCGACATTGAACACAGTGCAAAGTACGGCCTGTCCTTTTGTTCACGCGCCATAAGAGGGACACA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2967;ASM_Start=2893;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t37795\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2213;ASM_Start=2136;ASM_Strand=+;END=37872\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t37873\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=845;ASM_Start=845;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t37874\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7181;ASM_Start=6617;ASM_Strand=+;END=38438\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t38439\t.\tTACACGCTTTCGTGTTTCACGCCCAAGATATCCTACCGTCCTCGTTTCAACAGAACGGACGCGGCGGCTACCAGGGAGTTGCGGGACGACATGTGGCCA\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3471;ASM_Start=3471;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t38538\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1084;ASM_Start=257;ASM_Strand=+;END=39365\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t39366\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4250;ASM_Start=4250;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t39367\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5185;ASM_Start=4270;ASM_Strand=+;END=40282\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t40283\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1397;ASM_Start=1397;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "1\t40284\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8836;ASM_Start=8813;ASM_Strand=+;END=40307\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "1\t40308\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4006;ASM_Start=4006;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t1\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6936;ASM_Start=6271;ASM_Strand=+;END=666\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t667\t.\tTCCTACATTGGTGCAT\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1447;ASM_Start=1447;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t683\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7899;ASM_Start=7590;ASM_Strand=+;END=992\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t993\t.\tTCCTGTTTCTCGTGTATTGTGTTGGCTATAACAACGATGGA\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6936;ASM_Start=6936;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t1034\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8346;ASM_Start=8197;ASM_Strand=+;END=1183\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t1184\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4522;ASM_Start=4522;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t1185\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3524;ASM_Start=2857;ASM_Strand=+;END=1852\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t1853\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1148;ASM_Start=1148;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t1854\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1995;ASM_Start=1314;ASM_Strand=+;END=2535\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t2718\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9814;ASM_Start=9114;ASM_Strand=+;END=3418\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t3419\t.\tT\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7435;ASM_Start=7435;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t3420\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8750;ASM_Start=8266;ASM_Strand=+;END=3904\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t3905\t.\tT\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9012;ASM_Start=9012;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t3906\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1046;ASM_Start=421;ASM_Strand=+;END=4531\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t4532\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9334;ASM_Start=9334;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t4533\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=820;ASM_Start=573;ASM_Strand=+;END=4780\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t4781\t.\tG\tCTGGGAAGACTTCAGA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6774;ASM_Start=6760;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t4782\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=10079;ASM_Start=9896;ASM_Strand=+;END=4965\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t4966\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7302;ASM_Start=7302;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t4967\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=10586;ASM_Start=9661;ASM_Strand=+;END=5892\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t5893\t.\tT\tAAATCTAACCTACCCCACTTTCAAGCGCAAAAGCTCGCTACGGTTGCCATTGCCATATGATGTTCATCGTCACCACCGCGAGAACTCGGTTTTCGGCGT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4438;ASM_Start=4341;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t5894\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5004;ASM_Start=4389;ASM_Strand=+;END=6509\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t6510\t.\tGGATTCGGCGTTGTGCTAAGCCGCAGCGACGAGACTCAAAGT\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6002;ASM_Start=6002;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t6552\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9793;ASM_Start=9115;ASM_Strand=+;END=7230\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t7231\t.\tGGGCGATCCTGCCCACGGTGCTCTGCGTATGATAAGCGCATGGCAGGGCATT\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9232;ASM_Start=9232;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t7283\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6038;ASM_Start=5685;ASM_Strand=+;END=7636\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t7637\t.\tT\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1550;ASM_Start=1550;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t7638\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2509;ASM_Start=2377;ASM_Strand=+;END=7770\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t7771\t.\tC\tGAGGATGTTTCCGCACGAGCTGAGTGGCATGAGTC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7403;ASM_Start=7370;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t7772\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6716;ASM_Start=5837;ASM_Strand=+;END=8651\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t8652\t.\tC\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2628;ASM_Start=2628;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t8653\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7145;ASM_Start=7131;ASM_Strand=+;END=8667\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t8668\t.\tC\tTCGGTGCG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8285;ASM_Start=8279;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t8669\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3169;ASM_Start=2987;ASM_Strand=+;END=8851\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t8852\t.\tACAAAGCATGGGCGGCATGATACCTCGCTGCATGGGGTCATGAAGGTCGAACGCCGTGGCTGTAAACCGGCGTTAGTCACTATTATTC\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=560;ASM_Start=560;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t8940\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5672;ASM_Start=5049;ASM_Strand=+;END=9563\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t9978\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7841;ASM_Start=6853;ASM_Strand=+;END=10966\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t10967\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1415;ASM_Start=1415;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t10968\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2313;ASM_Start=1433;ASM_Strand=+;END=11848\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t11849\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5379;ASM_Start=5379;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t11850\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9861;ASM_Start=9387;ASM_Strand=+;END=12324\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t12325\t.\tT\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2405;ASM_Start=2405;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t12326\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1596;ASM_Start=1044;ASM_Strand=+;END=12878\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t12879\t.\tTCCCAAGAATGGACTGCTG\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8753;ASM_Start=8753;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t12898\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7269;ASM_Start=6839;ASM_Strand=+;END=13328\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t13329\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9452;ASM_Start=9452;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t13330\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7149;ASM_Start=6217;ASM_Strand=+;END=14262\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t14263\t.\tA\tATCCCAAGCCACGGATAAGGAACGTTGCAAGCTATAGTAAAGAGAGTAACGTCTAAGCGGAAGGTTGAAGTCCCGGGTGGAG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6125;ASM_Start=6045;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t14264\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5551;ASM_Start=5486;ASM_Strand=+;END=14329\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t14330\t.\tA\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3559;ASM_Start=3559;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t14331\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5683;ASM_Start=5179;ASM_Strand=+;END=14835\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t14836\t.\tT\tACAATTCGCATGATCTTGGGATTAGTTGCACCCGCGTACGAGCCCGGTTTTGCTACGTTCTCTGTTCTTAATGGCG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3416;ASM_Start=3342;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t14837\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5964;ASM_Start=5926;ASM_Strand=+;END=14875\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t14876\t.\tC\tAACCGACCTTGGA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4460;ASM_Start=4449;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t14877\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9078;ASM_Start=8138;ASM_Strand=+;END=15817\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t15818\t.\tT\tCAATATACACGACGCTCACTAGGGCCCAACTTGAAGGGTTACAGTACGTCTTGGGTGTGAAGACAGAACTTCAAAGGAAAGGCCGG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5454;ASM_Start=5370;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t15819\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2172;ASM_Start=1603;ASM_Strand=+;END=16388\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t16389\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8893;ASM_Start=8893;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t16390\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4675;ASM_Start=3900;ASM_Strand=+;END=17165\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t17166\t.\tT\tCGGTTTGGCAGGGATCTGACGCCTATCTAGACGGTGGTGAGAG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9087;ASM_Start=9046;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t17167\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=10124;ASM_Start=9569;ASM_Strand=+;END=17722\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t17723\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9990;ASM_Start=9990;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t17724\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6892;ASM_Start=6020;ASM_Strand=+;END=18596\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t18597\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4545;ASM_Start=4545;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t18598\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2685;ASM_Start=1985;ASM_Strand=+;END=19298\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t19299\t.\tA\tGTGCCAATTCAC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6328;ASM_Start=6318;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t19300\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3293;ASM_Start=2897;ASM_Strand=+;END=19696\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t19697\t.\tG\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4229;ASM_Start=4229;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t19698\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2671;ASM_Start=2234;ASM_Strand=+;END=20135\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t20136\t.\tG\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7936;ASM_Start=7936;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t20137\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7556;ASM_Start=7553;ASM_Strand=+;END=20140\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t20141\t.\tG\tGGCATACACGACAGGCATCTCTTGCACCCGTGGGATCATTAAGCTGTCGGCAAAGCACCCAGTACGCATTGAATGCCCGAGACAGTCTGCAGCCAG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6161;ASM_Start=6067;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t20142\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2009;ASM_Start=1247;ASM_Strand=+;END=20904\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t20905\t.\tT\tGGTGATATTAATCTTCAGTGGAGC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=565;ASM_Start=543;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t20906\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=756;ASM_Start=258;ASM_Strand=+;END=21404\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t21405\t.\tC\tAACAGATGTTACAGTGAGTACGCCATGCCATGGCTTTCGTGGACTGGCTCAGGCCGTTGACACCTCGAACAGTTCTCTCATGATGGAGAGGA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5900;ASM_Start=5810;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t21406\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6459;ASM_Start=5814;ASM_Strand=+;END=22051\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t22052\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9014;ASM_Start=9014;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t22053\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2340;ASM_Start=1801;ASM_Strand=+;END=22592\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t22593\t.\tA\tGTGAACTTACAACGCCAGTATCACTTTAAGTTCCTATTATGCTGAAGGGCA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5234;ASM_Start=5185;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t22594\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5158;ASM_Start=4186;ASM_Strand=+;END=23566\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t23567\t.\tGGGTACCAGTGTTTGTGGACCTAAACCCCTTGTCCCAACTCATAAGGGGGTCAG\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4486;ASM_Start=4486;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t23621\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5290;ASM_Start=4746;ASM_Strand=+;END=24165\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t24166\t.\tA\tCGGTCGGCAAACTCGAGGGCATCTCGTCAAAGCCAAAGGGATCCAATTCGGGT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5673;ASM_Start=5622;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t24167\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8762;ASM_Start=8463;ASM_Strand=+;END=24466\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t24467\t.\tGTTGTGCCGTAGACGCCATCCGCGCTCTCTTACGCCTAGGATGGCGTGGCAGGCGT\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=227;ASM_Start=227;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t24523\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=654;ASM_Start=259;ASM_Strand=+;END=24918\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t24919\t.\tA\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5046;ASM_Start=5046;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t24920\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9306;ASM_Start=8387;ASM_Strand=+;END=25839\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t25840\t.\tGGTATTGACGGTGGTTCCTCAGTTTGGCAAACCT\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1888;ASM_Start=1888;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t25874\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6044;ASM_Start=5645;ASM_Strand=+;END=26273\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t26274\t.\tCGGAAGATGAAGTCACCC\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1417;ASM_Start=1417;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t26292\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7589;ASM_Start=6775;ASM_Strand=+;END=27106\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t27107\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4189;ASM_Start=4189;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t27108\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3753;ASM_Start=3444;ASM_Strand=+;END=27417\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t27418\t.\tCCCCATGGAAGTGACACCACCAGGGTCTGCGCTTGGGTATTC\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=816;ASM_Start=816;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t27460\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9785;ASM_Start=9485;ASM_Strand=+;END=27760\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t27761\t.\tC\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2749;ASM_Start=2749;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t27762\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8427;ASM_Start=7898;ASM_Strand=+;END=28291\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t28292\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2316;ASM_Start=2316;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t28293\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2764;ASM_Start=1941;ASM_Strand=+;END=29116\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t29117\t.\tA\tC,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1352;ASM_Start=1352;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t29118\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2533;ASM_Start=2180;ASM_Strand=+;END=29471\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t29472\t.\tAATCGGCTTGTTGGGTCTGCTTCGATCCAGCCGGGCAGACGAGGAGGCCGATTACTTATTTAGAGACTTGACGG\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7399;ASM_Start=7399;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t29546\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6869;ASM_Start=6839;ASM_Strand=+;END=29576\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t29577\t.\tC\tA,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2947;ASM_Start=2947;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t29578\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9254;ASM_Start=9209;ASM_Strand=+;END=29623\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t29624\t.\tC\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7699;ASM_Start=7699;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t29625\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2369;ASM_Start=1977;ASM_Strand=+;END=30017\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t30556\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2796;ASM_Start=2504;ASM_Strand=+;END=30848\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t30849\t.\tCTGTCTTTGGCATTTCGCACTGGCTAGTCACCGG\tG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6910;ASM_Start=6910;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t30883\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=9069;ASM_Start=8986;ASM_Strand=+;END=30966\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t30967\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1578;ASM_Start=1578;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t30968\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2301;ASM_Start=1811;ASM_Strand=+;END=31458\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t31459\t.\tCAACATCGGCGGTCATATAGCTTATCAGCACAACAGAAACAGCAATGCTTTGA\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4281;ASM_Start=4281;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t31512\t.\tT\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=3439;ASM_Start=2697;ASM_Strand=+;END=32254\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t32255\t.\tG\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6146;ASM_Start=6146;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n"  +
                        "2\t32256\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=4342;ASM_Start=3777;ASM_Strand=+;END=32821\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"  +
                        "2\t32822\t.\tA\tGTAGTTCGCTCTCAGGCACCGGGCGCGCTATCGTGGCCTGTAGAGGCTGCGGCGGTTACGAATTGCTGG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5983;ASM_Start=5916;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n")
            }


            val builder = ProcessBuilder("bgzip", "$outputFile")
            val process = builder.start()
            var error = process.waitFor()

            val builder2 = ProcessBuilder("tabix", "$outputFile.gz")
            val process2 = builder2.start()
            error += process2.waitFor()

            if (error != 0) {
                println("Something went wrong creating the testing gvcf. Check that bgzip and tabix are installed on this machine.")
            }
        }

        /**
         * the truth output for dummyGVCF and RefGVCF
         */
        private fun getTruthTSV(): List<String> {
            return listOf("taxa\tchrom\trefLength\tnumSNPs\tnumIns\tnumDel\tnumNs\tnumBasesInserted\tnumBasesDeleted\tpercentIdentityWithRef\tpercentMappedToRef\tmeanInsertionSize\tmedianInsertionSize\tlargestInsertion\tmeanDeletionSize\tmedianDeletionSize\tlargestDeletion",
                    "sampleName\tALL\t73130\t120\t37\t31\t0\t1906\t1692\t0.9280049227403254\t0.9527827157117462\t51.513513513513516\t55.0\t98\t54.58064516129032\t51.0\t99",
                    "sampleName\t1\t40308\t67\t21\t17\t0\t1067\t1083\t0.9139376798650392\t0.9424679964275082\t50.80952380952381\t58.0\t98\t63.705882352941174\t68.0\t99",
                    "sampleName\t2\t32822\t53\t16\t14\t0\t839\t609\t0.9452806044726099\t0.9654500030467369\t52.4375\t51.0\t98\t43.5\t41.0\t87",
                    "Ref\tALL\t73130\t0\t0\t0\t0\t0\t0\t1.0\t1.0\t0.0\t0.0\t0\t0.0\t0.0\t0",
                    "Ref\t1\t40308\t0\t0\t0\t0\t0\t0\t1.0\t1.0\t0.0\t0.0\t0\t0.0\t0.0\t0",
                    "Ref\t2\t32822\t0\t0\t0\t0\t0\t0\t1.0\t1.0\t0.0\t0.0\t0\t0.0\t0.0\t0")
        }

        private fun createNGVCF(outputFile: String) {
            getBufferedWriter(outputFile).use { output ->
                output.write(
                    "##fileformat=VCFv4.2\n" +
                            "##FORMAT=<ID=AD,Number=3,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n" +
                            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n" +
                            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" +
                            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                            "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n" +
                            "##INFO=<ID=AF,Number=3,Type=Integer,Description=\"Allele Frequency\">\n" +
                            "##INFO=<ID=ASM_Chr,Number=1,Type=String,Description=\"Assembly chromosome\">\n" +
                            "##INFO=<ID=ASM_End,Number=1,Type=Integer,Description=\"Assembly end position\">\n" +
                            "##INFO=<ID=ASM_Start,Number=1,Type=Integer,Description=\"Assembly start position\">\n" +
                            "##INFO=<ID=ASM_Strand,Number=1,Type=String,Description=\"Assembly strand\">\n" +
                            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n" +
                            "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n" +
                            "##contig=<ID=1,length=40308>\n" +
                            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleName\n" +
                            "1\t1\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2864;ASM_Start=2857;ASM_Strand=+;END=8\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "1\t9\t.\tA\tTATTCTGCGGAGAGCCTTAGTGAGTATAATGAAGTCACATAACTGTTTTGCTACCCTTCTCTCGGACTCCCCCGGATGAT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=7181;ASM_Start=7103;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                            "1\t10\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8036;ASM_Start=7987;ASM_Strand=+;END=59\tGT:AD:DP:PL\t.:0,30,0:30:90,0,90\n" +
                            "1\t60\t.\tGGCCCGTGTGAGATACTCTTATAGCGCAGGCTGGACAATACGGAGC\tT,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=8829;ASM_Start=8829;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                            "1\t106\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1166;ASM_Start=697;ASM_Strand=+;END=575\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n" +
                            "1\t576\t.\tG\tACCGATCAATCCTAATGTGTCTANNNNNNNNNNNTTAAGCGACGTAG,<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=702;ASM_Start=638;ASM_Strand=+\tGT:AD:DP:PL\t1:0,30,0:30:90,0,90\n" +
                            "1\t577\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=2375;ASM_Start=2366;ASM_Strand=+;END=586\tGT:AD:DP:PL\t0:0,30,0:30:90,0,90\n"
                )
            }
        }

    }

}