package net.maizegenetics.phgv2.cli

import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.testing.test
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.vcf.VCFFileReader
import io.kotest.common.runBlocking
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.utils.*
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertNotNull
import kotlin.test.assertNull
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
class Hvcf2GvcfTest {
    companion object {

        val diploidFixtureDir = "${TestExtension.testVCFDir}/diploidFixtures"

        val prebuiltDataDir = "data/test/hvcf2vcf"
        val prebuiltAsmHvcfDir = "$prebuiltDataDir/asmHvcfs"
        val prebuiltGvcfDir = "$prebuiltDataDir/temp"
        val prebuiltImputeHvcfDir = "$prebuiltDataDir/imputeHvcfs"

        @JvmStatic
        @BeforeAll
        fun setup() {
            val dbPath = TestExtension.testTileDBURI
            val refName = "Ref"
            val refUrl = TestExtension.refURL

            val ranges = "data/test/smallseq/anchors.bed"
            val refFasta = "data/test/smallseq/Ref.fa"

            File(TestExtension.testVCFDir).mkdirs()
            File(TestExtension.testTileDBURI).mkdirs()
            Initdb().createDataSets(dbPath, "")

            // Create the agc compressed file
            println("testSimpleHvcf2Gvcf:running agcCompress")
            val agcCompress = AgcCompress()
            var agcResult =
                agcCompress.test("--fasta-list ${TestExtension.smallseqAssembliesListFile} --db-path ${dbPath} --reference-file ${TestExtension.smallseqRefFile}")
            println(agcResult.output)

            println("testSimpleHvcf2Gvcf:running CreateRefVcf")
            var result =
                CreateRefVcf().test("--bed $ranges --reference-name $refName --reference-file $refFasta --reference-url ${refUrl} --db-path $dbPath")
            assertEquals(0, result.statusCode)

            // Run alignAssemblies test to get MAF files.
            // Run createMAFVCf on the assemblies LineA and LIneB to get
            // data into the db.

            println("testSimpleHvcf2Gvcf: running AlignAssemblies")
            val alignAssemblies = AlignAssemblies()

            result = alignAssemblies.test(
                "--gff ${TestExtension.smallseqAnchorsGffFile} --reference-file ${TestExtension.smallseqRefFile} " +
                        "--assembly-file-list ${TestExtension.smallseqAssembliesListFile} -o ${TestExtension.tempDir} --total-threads 1 --in-parallel 1"
            )

            println("testSimpleHvcf2Gvcf: result output: ${result.output}")
            assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

            // Load assemblies using CreateMafVcf - creates and loads gvcf and hvcf
            println("testSimpleHvcf2Gvcf: running CreateMafVcf")
            val createMafVcf = CreateMafVcf()
            result =
                createMafVcf.test("--db-path ${dbPath} --bed ${ranges} --reference-file ${refFasta} --maf-dir ${TestExtension.tempDir} -o ${TestExtension.testVCFDir}")
            println(result.output)

            // Need to load the vcf now!
            // Load the vcf files into the tiledb dataset
            println("testSimpleHvcf2Gvcf: running LoadVcf")
            val loadVcf = LoadVcf()
            result = loadVcf.test("--db-path ${dbPath} --vcf-dir ${TestExtension.testVCFDir}")

            generateDiploidHvcfFixtures(dbPath, refFasta)
        }

        /**
         * Generate diploid hVCF fixture files that can be loaded from disk by subsequent tests.
         * Creates heterozygous (LineA/LineB), homozygous (LineA/LineA), and
         * ref-gamete (Ref/LineA) diploid hVCF files.
         */
        private fun generateDiploidHvcfFixtures(dbPath: String, refFasta: String) {
            val refSeq = CreateMafVcf().buildRefGenomeSeq(refFasta)
            File(diploidFixtureDir).mkdirs()

            val hvcfFiles = mutableListOf<String>()
            File(TestExtension.testVCFDir).listFiles()!!
                .filter { it.name.endsWith(".h.vcf") || it.name.endsWith(".h.vcf.gz") }
                .forEach { hvcfFiles.add(it.path) }
            File("${dbPath}/hvcf_files").listFiles()!!
                .filter { it.name.endsWith(".h.vcf") || it.name.endsWith(".h.vcf.gz") }
                .forEach { hvcfFiles.add(it.path) }
            val graph = HaplotypeGraph(hvcfFiles.distinct())

            val sampleConfigs = listOf(
                Triple("DiploidHet_LineA_LineB", SampleGamete("LineA"), SampleGamete("LineB")),
                Triple("DiploidHom_LineA", SampleGamete("LineA"), SampleGamete("LineA")),
                Triple("DiploidRef_Ref_LineA", SampleGamete("Ref"), SampleGamete("LineA"))
            )

            for ((sampleName, gamete1, gamete2) in sampleConfigs) {
                buildDiploidHvcfFile(graph, sampleName, gamete1, gamete2, refSeq, "$diploidFixtureDir/$sampleName.h.vcf")
            }

            println("Generated diploid hVCF fixtures in $diploidFixtureDir")
        }

        private fun buildDiploidHvcfFile(
            graph: HaplotypeGraph,
            sampleName: String,
            gamete1: SampleGamete,
            gamete2: SampleGamete,
            refSeq: Map<String, NucSeq>,
            outputFile: String
        ) {
            val altHeaders = mutableListOf<AltHeaderMetaData>()
            val variantContexts = mutableListOf<VariantContext>()

            for (range in graph.ranges()) {
                val hapId1 = graph.sampleToHapId(range, gamete1)
                val hapId2 = graph.sampleToHapId(range, gamete2)
                val refAllele = refSeq[range.contig]!![range.start - 1].toString()
                val startPos = Position(range.contig, range.start)
                val endPos = Position(range.contig, range.end)

                variantContexts.add(createDiploidHVCFRecord(sampleName, startPos, endPos, listOf(hapId1, hapId2), refAllele))
                listOfNotNull(hapId1, hapId2).distinct().mapNotNull { graph.altHeader(it) }
                    .forEach { altHeaders.add(it) }
            }

            val headerSet = altHeaders.distinctBy { it.id }.map { altHeaderMetadataToVCFHeaderLine(it) }.toSet()
            exportVariantContext(sampleName, variantContexts, outputFile, refSeq, headerSet)
        }

        @JvmStatic
        @AfterAll
        fun tearDown() {
            File(TestExtension.tempDir).deleteRecursively()
            File(TestExtension.testVCFDir).deleteRecursively()
            File(TestExtension.testTileDBURI).deleteRecursively()

        }
    }

    @Test
    fun testCliktParams() {
        val hvcf2gvcf = Hvcf2Gvcf()

        // There are only 3 required parameters - test for missing each one
        assertThrows<IllegalArgumentException> {
            //Check that an error is thrown when the dbPath folder does not exist
            hvcf2gvcf.test("--hvcf-dir ${TestExtension.testVCFDir} --reference-file ${TestExtension.testRefFasta} ")
        }

        val resultMissingRef =
            hvcf2gvcf.test("--db-path ${TestExtension.testTileDBURI}  --hvcf-dir ${TestExtension.testMafDir}")
        assertEquals(resultMissingRef.statusCode, 1)
        assertEquals(
            "Usage: hvcf2gvcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --reference-file\n", resultMissingRef.output
        )


        val resultMissingHvcfDir =
            hvcf2gvcf.test("--db-path ${TestExtension.testTileDBURI}  --reference-file ${TestExtension.testRefFasta}")
        assertEquals(resultMissingHvcfDir.statusCode, 1)
        assertEquals(
            "Usage: hvcf2gvcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --hvcf-dir\n", resultMissingHvcfDir.output
        )

    }

    @Test
    fun testSimpleHvcf2Gvcf() {
        // This is a basic test.  We copy an hvcf file created from CreateMafVcf to a new location
        // and run that through hvcf2gvcf.  We then verify that the output file exists.

        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/smallseq/Ref.fa"

        // Copy the $TestExtension.testVCFDir/LineB.h.vcf.gz to $TestExtension.testVCFDir/LineBPath.h.vcf.gz
        // Make directory ${TestExtension.testVCFDir}/testOutputGVCFDir

        println("\nNOW ... running hvcf2gvcf")
        val testGVCFdir = "${TestExtension.testVCFDir}/testOutputGVCFDir"
        File(testGVCFdir).mkdirs()
        val lineBPathHvcf = "${testGVCFdir}/LineBPath.h.vcf.gz"
        val lineBHvcf = "${TestExtension.testVCFDir}/LineB.h.vcf.gz"
        File(lineBHvcf).copyTo(File(lineBPathHvcf))

        // Run hvcf2gvcf on the copied file
        val hvcf2gvcf = Hvcf2Gvcf()
        val result =
            hvcf2gvcf.test("--db-path ${dbPath} --hvcf-dir $testGVCFdir --output-dir ${testGVCFdir} --reference-file ${refFasta}")
        // verify the output file exists, which will be LineBPath.g.vcf
        assertTrue(File("${testGVCFdir}/LineBPath.g.vcf").exists())
        // Verify header lines
        val lines = File("${testGVCFdir}/LineBPath.g.vcf").readLines()
        assertEquals(17, lines.filter { it.startsWith("#") }.size)
        assertEquals("##fileformat=VCFv4.2", lines[0])
        assertTrue(`lines`.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tLineBPath"))

        // The LineBPath.g.vcf file will contain lines representing the beginning of the ref ranges, e.g.1\t1001 1\t5501 1\t6501 2\t1001 2\t5501 2\t6501
        // These are not in the original LineB.vcf file, but are added by hvcf2gvcf to represent the beginning of the ref ranges as the
        // h.vcf files from which it is made shows those positions.  These entries are created as refBlocks from the beginning of the refRange
        assertTrue(lines.contains("1\t1001\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=1001;ASM_Start=1001;ASM_Strand=+;END=1001\tGT:AD:DP:PL\t0:30,0:30:0,90,90"))
        assertTrue(lines.contains("1\t5501\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5561;ASM_Start=5501;ASM_Strand=+;END=5561\tGT:AD:DP:PL\t0:30,0:30:0,90,90"))
        assertTrue(lines.contains("1\t6501\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6531;ASM_Start=6501;ASM_Strand=+;END=6531\tGT:AD:DP:PL\t0:30,0:30:0,90,90"))
        assertTrue(lines.contains("2\t1001\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=2;ASM_End=1064;ASM_Start=1001;ASM_Strand=+;END=1064\tGT:AD:DP:PL\t0:30,0:30:0,90,90"))
        assertTrue(lines.contains("2\t5501\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=2;ASM_End=5508;ASM_Start=5501;ASM_Strand=+;END=5508\tGT:AD:DP:PL\t0:30,0:30:0,90,90"))
        assertTrue(lines.contains("2\t6501\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=2;ASM_End=6501;ASM_Start=6501;ASM_Strand=+;END=6501\tGT:AD:DP:PL\t0:30,0:30:0,90,90"))

        println("done !!")


    }

    @Test
    fun testPathHvcf2Gvcf() {
        // This is a test of the Hvcf2Gvcf:run function.  An hvcf file created from the
        // smallSeq FindPathsTest is copied to a new location and run through hvcf2gvcf.

        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/smallseq/Ref.fa"

        // The tiledb datasets have been created and populated in the setup functions.
        // Using hvcf file TestLine2.h.vcf, created from the FindPathsTest junit
        // This has been stored to test/smallseq folder.  Its contents represent sequence
        // from 2 assemblies, LineA and LineB.


        println("\ntestPathHvcf2Gvcf... running hvcf2gvcf")
        val testGVCFdir = "${TestExtension.testVCFDir}/testOutputGVCFDir"
        File(testGVCFdir).deleteRecursively() // make sure folder is clean
        File(testGVCFdir).mkdirs() // start fresh
        val testLine2Hvcf = "${testGVCFdir}/TestLine2.h.vcf.gz"
        val line2Hvcf = "data/test/smallseq/TestLine2.h.vcf.gz"
        File(line2Hvcf).copyTo(File(testLine2Hvcf))

        // Initial verification that neither LineA.vcf nor LineB.vcf exist in the testGVCFdir
        // THese files will be exported and added after the hvcf2gvcf run
        assertTrue(!File("${testGVCFdir}/LineA.vcf").exists())
        assertTrue(!File("${testGVCFdir}/LineB.vcf").exists())

        // Run hvcf2gvcf on the copied file
        val hvcf2gvcf = Hvcf2Gvcf()
        val result =
            hvcf2gvcf.test("--db-path ${dbPath} --hvcf-dir $testGVCFdir --output-dir ${testGVCFdir} --reference-file ${refFasta}")
        // verify the output file exists, which will be LineBPath.g.vcf
        assertTrue(File("${testGVCFdir}/TestLine2.g.vcf").exists())

        // verify the vcfs for LineA and LineB were exported
        assertTrue(File("${testGVCFdir}/LineA.vcf").exists())
        assertTrue(File("${testGVCFdir}/LineB.vcf").exists())

        // Verify header lines
        val lines = File("${testGVCFdir}/TestLine2.g.vcf").readLines()
        assertEquals(17, lines.filter { it.startsWith("#") }.size)
        assertEquals("##fileformat=VCFv4.2", lines[0])
        assertTrue(`lines`.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTestLine2"))

        // compare to the "truth" file
        // Regarding the truth file: Because the LineA.vcf and LineB.vcf files have no data past position 50300
        // in them, the TestLine2.g.vcf file will have no data past position 50300 in it.
        // The entries for positions between 1-27500 on both chrosomes 1 and 2 all come from LineA.vcf
        // The entries for positions between 27501-50300 on both chromosomes 1 and 2 all come from LineB.vcf.
        // In additiona, the TestLine2.g.vcf file has entries at the beginning of each reference range as the
        // hvcf file lists the region beginning with the refRange beginning. RefBLocks were created for these
        // positions from refBLock beginning to first variant in the vcf file.
        CreateMafVCFTest().compareTwoGVCFFiles(
            "data/test/hvcfToGvcf/TestLine2_truth.g.vcf",
            "${testGVCFdir}/TestLine2.g.vcf"
        )

        println("testPathHvcf2Gvcf done !!")
    }

    @Test
    fun testPathHvcf2Gvcf_withExports() {
        // This is a test of the Hvcf2Gvcf:run function.  An hvcf file created from the
        // smallSeq FindPathsTest is copied to a new location and run through hvcf2gvcf.

        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/smallseq/Ref.fa"

        // The tiledb datasets have been created and populated in the setup functions.
        // Using hvcf file TestLine2.h.vcf, created from the FindPathsTest junit
        // This has been stored to test/smallseq folder.  Its contents represent sequence
        // from 2 assemblies, LineA and LineB.


        println("\ntestPathHvcf2Gvcf... running hvcf2gvcf")
        val testGVCFdir = "${TestExtension.testVCFDir}/testOutputGVCFDir"
        File(testGVCFdir).deleteRecursively() // make sure folder is clean
        File(testGVCFdir).mkdirs() // start fresh
        val testLine2Hvcf = "${testGVCFdir}/TestLine2.h.vcf.gz"
        val line2Hvcf = "data/test/smallseq/TestLine2.h.vcf.gz"
        File(line2Hvcf).copyTo(File(testLine2Hvcf))

        // export LineA.vcf.  This is to test the code does not export
        // when the file exists.  If the -O v option is used in the call to
        //  tiledb export, the extension from tiledb is .vcf.
        // If it os NOT used, the extension seems to be
        // what was used in the original file.  In this case, the extension is .g.vcf
        // In the hvcf2gvcf code I call tiledb with the -O v option.  But the
        // ExportVCF() clikt command does not, so we get extension .g.vcf
        // THis ends up being useful for determining at which point the files
        // were exported.
        val resultExport = ExportVcf().test(
            "--db-path ${dbPath} --sample-names LineA  --dataset-type gvcf -o ${testGVCFdir}"
        )
        println("Export LineA.vcf result: ${resultExport.output}")

        // Initial verification after export that LineA.vcf exist but LineB.vcf does not exist in the testGVCFdir
        // LineB should be exported
        assertTrue(File("${testGVCFdir}/LineA.g.vcf").exists())
        assertTrue(!File("${testGVCFdir}/LineB.vcf").exists())

        // Run hvcf2gvcf on the copied file
        val hvcf2gvcf = Hvcf2Gvcf()
        val result =
            hvcf2gvcf.test("--db-path ${dbPath} --hvcf-dir $testGVCFdir --output-dir ${testGVCFdir} --reference-file ${refFasta}")
        // verify the output file exists, which will be LineBPath.g.vcf
        assertTrue(File("${testGVCFdir}/TestLine2.g.vcf").exists())

        // verify the vcfs for LineA and LineB were exported
        assertTrue(File("${testGVCFdir}/LineA.g.vcf").exists())
        assertTrue(File("${testGVCFdir}/LineB.vcf").exists())

        // Verify header lines
        val lines = File("${testGVCFdir}/TestLine2.g.vcf").readLines()
        assertEquals(17, lines.filter { it.startsWith("#") }.size)
        assertEquals("##fileformat=VCFv4.2", lines[0])
        assertTrue(`lines`.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTestLine2"))

        // compare to the "truth" file
        // Regarding the truth file: Because the LineA.vcf and LineB.vcf files have no data past position 50300
        // in them, the TestLine2.g.vcf file will have no data past position 50300 in it.
        // The entries for positions between 1-27500 on both chrosomes 1 and 2 all come from LineA.vcf
        // The entries for positions between 27501-50300 on both chromosomes 1 and 2 all come from LineB.vcf.
        // In additiona, the TestLine2.g.vcf file has entries at the beginning of each reference range as the
        // hvcf file lists the region beginning with the refRange beginning. RefBLocks were created for these
        // positions from refBLock beginning to first variant in the vcf file.
        CreateMafVCFTest().compareTwoGVCFFiles(
            "data/test/hvcfToGvcf/TestLine2_truth.g.vcf",
            "${testGVCFdir}/TestLine2.g.vcf"
        )

        println("testPathHvcf2Gvcf done !!")
    }

    @Test
    fun testExportGvcfs() {
        // THis test verifies batching of the export files from tiledb
        // val exportSuccess = exportGvcfFiles(missingSamples, outputDir, dbPath, condaEnvPrefix, 5)
        val dbPath = TestExtension.testTileDBURI
        val outputDir = TestExtension.tempDir + "/vcfDir/exportedFiles"
        // make the exporedFiles folder
        File(outputDir).mkdirs()
        val condaEnvPrefix = ""
        val missingSamples = setOf("LineA", "LineB")
        val batchSize = 1 // using 1 as there are only 2 samples in the test tiledbURI
        runBlocking {
            val exportSuccess =
                Hvcf2Gvcf().exportSamplesGvcfFiles(missingSamples, outputDir, dbPath, condaEnvPrefix, batchSize)
        }
        // verify the gvcf files exist in the outputDir
        assertTrue(File("$outputDir/LineA.vcf").exists())
        assertTrue(File("$outputDir/LineB.vcf").exists())

        //remove the files written to the outputDir
        File("$outputDir/LineA.vcf").delete()
        File("$outputDir/LineB.vcf").delete()

        runBlocking {
            //  export with all in the same batch, batch size = 5
            val exportSuccess2 =
                Hvcf2Gvcf().exportSamplesGvcfFiles(missingSamples, outputDir, dbPath, condaEnvPrefix, 5)
        }

        // verify the gvcf files exist in the outputDir
        assertTrue(File("$outputDir/LineA.vcf").exists())
        assertTrue(File("$outputDir/LineB.vcf").exists())


    }

    @Test
    fun testCreateRefGvcf() {
        // This is a test of the Hvcf2Gvcf:createRefGvcf function.  We create a reference gvcf file
        // and verify it has data for all reference ranges in the bed file.

        val outputDir = TestExtension.tempDir
        val refName = "Ref"

        val ranges = "data/test/smallseq/anchors.bed"
        val refFasta = "data/test/smallseq/Ref.fa"
        val refSeq = CreateMafVcf().buildRefGenomeSeq(refFasta)

        Hvcf2Gvcf().createRefGvcf(refName, ranges, refSeq, outputDir)
        assertTrue(File("$outputDir/$refName.vcf").exists())

        // read the file and verify that it has data for all reference ranges in the bed file
        val lines = File("$outputDir/$refName.vcf").readLines()
        val refRanges = File(ranges).readLines()


        // Verify there are 17 lines in the file that begin with "#" (header lines)
        assertEquals(17, lines.filter { it.startsWith("#") }.size)

        // Verify the first line is the file is "##fileformat=VCFv4.2"
        assertEquals("##fileformat=VCFv4.2", lines[0])

        // verify there are 20 lines in the created file that begin with "1" and 20 that begin with "2"
        assertEquals(20, lines.filter { it.startsWith("1") }.size)
        assertEquals(20, lines.filter { it.startsWith("2") }.size)

        val chrom1entry =
            "1\t1001\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5500;ASM_Start=1001;ASM_Strand=+;END=5500\tGT:AD:DP:PL\t0:30,0:30:0,90,90"
        val chom2entry =
            "2\t1001\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=2;ASM_End=5500;ASM_Start=1001;ASM_Strand=+;END=5500\tGT:AD:DP:PL\t0:30,0:30:0,90,90"

        // verify the file contains both the chrom1entry and the chrom2entry
        assertTrue(lines.contains(chrom1entry))
        assertTrue(lines.contains(chom2entry))
        println("ref gvcf file created: $outputDir/$refName.vcf")
    }

    @Test
    fun testResizeVcandASMPositions_ForwardStrand() {
        // Build a couple variant context objects
        // Make the simple, we don't care about depth or pl at this point,
        // we just want to see if we resize correctly.  But we do need genotypes
        val alleles = listOf(Allele.create("G", true), Allele.NON_REF_ALLELE)
        val genotype = listOf(alleles[0])

        val vc1 = VariantContextBuilder(".", "chr1", 1001, 1010, alleles)
        val currentGenotype = GenotypeBuilder("LineA", genotype).make()
        vc1.attribute("ASM_Chr", "1").attribute("ASM_Start", 1003).attribute("ASM_End", 1012)
            .attribute("ASM_Strand", "+")
        val firstVariant = vc1.genotypes(currentGenotype).make()

        // First test sending in single variant context
        // consider the ref range is position 1005-1015, but the variant is at 1001-1010
        // The starts are adjusted, the ends remain the same
        val resizedVc1 =
            Hvcf2Gvcf().resizeVCandASMpositions(Pair(firstVariant, firstVariant), Pair(1005, 1015), Pair("+", "+"))
        println("resizedVc1 ref/asm start positions: ${resizedVc1[0].first} ${resizedVc1[0].second}")
        println("resizedVc1 ref/asm end positions: ${resizedVc1[1].first} ${resizedVc1[1].second}")
        assertEquals(1005, resizedVc1[0].first) // ref start
        assertEquals(1010, resizedVc1[1].first) // ref end
        assertEquals(1007, resizedVc1[0].second) // asm start
        assertEquals(1012, resizedVc1[1].second) // asm end

        // Now test sending in two variants (ie last is different than first)
        // consider the ref range is position 1005-1015, but the first variant is 1001-1010 and
        // second variant is 1011-1020
        // Both start and end are adjusted
        val vc2 = VariantContextBuilder(".", "chr1", 1011, 1020, alleles)
        val currentGenotype2 = GenotypeBuilder("LineA", genotype).make()
        vc2.attribute("ASM_Chr", "1").attribute("ASM_Start", 1013).attribute("ASM_End", 1022)
            .attribute("ASM_Strand", "+")
        val secondVariant = vc2.genotypes(currentGenotype2).make()
        val resizedVc2 =
            Hvcf2Gvcf().resizeVCandASMpositions(Pair(firstVariant, secondVariant), Pair(1005, 1015), Pair("+", "+"))

        println("resizedVc2 ref/asm start positions: ${resizedVc2[0].first} ${resizedVc2[0].second}")
        println("resizedVc2 ref/asm end positions: ${resizedVc2[1].first} ${resizedVc2[1].second}")
        assertEquals(1005, resizedVc2[0].first) // variant start
        assertEquals(1015, resizedVc2[1].first) // variant end
        assertEquals(1007, resizedVc2[0].second) // asm start
        assertEquals(1017, resizedVc2[1].second) // asm end

    }

    @Test
    fun testResizeVcandASMPositions_reverseStrand() {
        // Note the BioKotlin created gvcfs have the ASM_Start value greater
        // than the ASM_End value when the strand is reverse.

        // Build a couple variant context objects
        // Make it simple, we don't care about depth or pl at this point,
        // we just want to see if we resize correctly.  But we do need genotypes
        val alleles = listOf(Allele.create("G", true), Allele.NON_REF_ALLELE)
        val genotype = listOf(alleles[0])

        val vc1 = VariantContextBuilder(".", "chr1", 1001, 1010, alleles)
        val currentGenotype = GenotypeBuilder("LineA", genotype).make()
        // WHen the strand is reverse, the start and end are reversed
        vc1.attribute("ASM_Chr", "1").attribute("ASM_Start", 1022).attribute("ASM_End", 1013)
            .attribute("ASM_Strand", "-")
        val firstVariant = vc1.genotypes(currentGenotype).make()

        // Test a single variant context
        // consider the ref range is position 1005-1010, but the variant is at 1001-1010
        // The starts remain, the ends are adjusted
        println("\nReverse strand test")
        val resizedVc1 =
            Hvcf2Gvcf().resizeVCandASMpositions(Pair(firstVariant, firstVariant), Pair(1005, 1010), Pair("-", "-"))

        println("resizedVc1 ref/asm start positions: ${resizedVc1[0].first} ${resizedVc1[0].second}")
        println("resizedVc1 ref/asm end positions: ${resizedVc1[1].first} ${resizedVc1[1].second}")

        assertEquals(1005, resizedVc1[0].first) // variant start changes
        assertEquals(1010, resizedVc1[1].first) // variant end
        assertEquals(1018, resizedVc1[0].second) // asm start changes
        assertEquals(1013, resizedVc1[1].second) // asm end

        val vc2 = VariantContextBuilder(".", "chr1", 1011, 1020, alleles)
        val currentGenotype2 = GenotypeBuilder("LineA", genotype).make()
        // NOTE: The ASM_Start for the second variant is less than the ASM_Start for the first variant.
        // This is how these sequences would appear in a real file, and is necessary for the
        // string of reverse strand entries to have their positions adjusted correctly.
        vc2.attribute("ASM_Chr", "1").attribute("ASM_Start", 1012).attribute("ASM_End", 1003)
            .attribute("ASM_Strand", "-")
        val secondVariant = vc2.genotypes(currentGenotype2).make()
        val resizedVc2 =
            Hvcf2Gvcf().resizeVCandASMpositions(Pair(firstVariant, secondVariant), Pair(1005, 1015), Pair("-", "-"))

        println("\nresizedVc2 ref/asm start positions: ${resizedVc2[0].first} ${resizedVc2[0].second}")
        println("resizedVc2 ref/asm end positions: ${resizedVc2[1].first} ${resizedVc2[1].second}")
        assertEquals(1005, resizedVc2[0].first) // variant start
        assertEquals(1015, resizedVc2[1].first) // variant end
        assertEquals(1018, resizedVc2[0].second) // asm start
        assertEquals(1008, resizedVc2[1].second) // asm end

    }
    
    @Test
    fun testProcessCurrentHVCFVariant_diploid() {
        val hvcf2gvcf = Hvcf2Gvcf()
        val refRange = ReferenceRange("1", 1001, 5500)
        val hapIdA = "hapA_checksum123"
        val hapIdB = "hapB_checksum456"

        val hapIdAndRangeToSampleMap = mapOf(
            Pair(refRange, hapIdA) to listOf(HvcfVariant(refRange, "LineA", hapIdA)),
            Pair(refRange, hapIdB) to listOf(HvcfVariant(refRange, "LineB", hapIdB))
        )

        val refAllele = Allele.create("A", true)
        val alt1 = Allele.create("<$hapIdA>", false)
        val alt2 = Allele.create("<$hapIdB>", false)
        val genotype = GenotypeBuilder("TestSample", listOf(alt1, alt2)).phased(true).make()
        val vc = VariantContextBuilder()
            .chr("1")
            .start(1001L)
            .stop(5500L)
            .alleles(listOf(refAllele, alt1, alt2))
            .genotypes(genotype)
            .attribute("END", 5500)
            .make()

        val (isHaploid, result) = hvcf2gvcf.processCurrentHVCFVariant(vc, hapIdAndRangeToSampleMap)

        assertFalse(isHaploid, "Diploid variant should not be haploid")
        assertNotNull(result.first, "First gamete variant should not be null")
        assertNotNull(result.second, "Second gamete variant should not be null")
        assertEquals("LineA", result.first!!.sampleName)
        assertEquals("LineB", result.second!!.sampleName)
        assertEquals(hapIdA, result.first!!.hapId)
        assertEquals(hapIdB, result.second!!.hapId)
    }

    @Test
    fun testProcessCurrentHVCFVariant_haploid() {
        val hvcf2gvcf = Hvcf2Gvcf()
        val refRange = ReferenceRange("1", 1001, 5500)
        val hapIdA = "hapA_checksum123"

        val hapIdAndRangeToSampleMap = mapOf(
            Pair(refRange, hapIdA) to listOf(HvcfVariant(refRange, "LineA", hapIdA))
        )

        val refAllele = Allele.create("A", true)
        val alt1 = Allele.create("<$hapIdA>", false)
        val genotype = GenotypeBuilder("TestSample", listOf(alt1)).make()
        val vc = VariantContextBuilder()
            .chr("1")
            .start(1001L)
            .stop(5500L)
            .alleles(listOf(refAllele, alt1))
            .genotypes(genotype)
            .attribute("END", 5500)
            .make()

        val (isHaploid, result) = hvcf2gvcf.processCurrentHVCFVariant(vc, hapIdAndRangeToSampleMap)

        assertTrue(isHaploid, "Haploid variant should be haploid")
        assertNotNull(result.first, "First gamete variant should not be null")
        assertNotNull(result.second, "Second gamete should also resolve (same as first for haploid)")
        assertEquals("LineA", result.first!!.sampleName)
        assertEquals("LineA", result.second!!.sampleName)
    }

    @Test
    fun testProcessCurrentHVCFVariant_diploidHomozygous() {
        val hvcf2gvcf = Hvcf2Gvcf()
        val refRange = ReferenceRange("1", 1001, 5500)
        val hapIdA = "hapA_checksum123"

        val hapIdAndRangeToSampleMap = mapOf(
            Pair(refRange, hapIdA) to listOf(HvcfVariant(refRange, "LineA", hapIdA))
        )

        val refAllele = Allele.create("A", true)
        val alt1 = Allele.create("<$hapIdA>", false)
        val genotype = GenotypeBuilder("TestSample", listOf(alt1, alt1)).phased(true).make()
        val vc = VariantContextBuilder()
            .chr("1")
            .start(1001L)
            .stop(5500L)
            .alleles(listOf(refAllele, alt1))
            .genotypes(genotype)
            .attribute("END", 5500)
            .make()

        val (isHaploid, result) = hvcf2gvcf.processCurrentHVCFVariant(vc, hapIdAndRangeToSampleMap)

        assertFalse(isHaploid, "Homozygous diploid should still be diploid")
        assertNotNull(result.first)
        assertNotNull(result.second)
        assertEquals("LineA", result.first!!.sampleName)
        assertEquals("LineA", result.second!!.sampleName)
        assertEquals(result.first!!.hapId, result.second!!.hapId)
    }

    @Test
    fun testProcessCurrentHVCFVariant_diploidMissingHapId() {
        val hvcf2gvcf = Hvcf2Gvcf()
        val refRange = ReferenceRange("1", 1001, 5500)
        val hapIdA = "hapA_checksum123"
        val hapIdB = "hapB_checksum456"

        val hapIdAndRangeToSampleMap = mapOf(
            Pair(refRange, hapIdA) to listOf(HvcfVariant(refRange, "LineA", hapIdA))
        )

        val refAllele = Allele.create("A", true)
        val alt1 = Allele.create("<$hapIdA>", false)
        val alt2 = Allele.create("<$hapIdB>", false)
        val genotype = GenotypeBuilder("TestSample", listOf(alt1, alt2)).phased(true).make()
        val vc = VariantContextBuilder()
            .chr("1")
            .start(1001L)
            .stop(5500L)
            .alleles(listOf(refAllele, alt1, alt2))
            .genotypes(genotype)
            .attribute("END", 5500)
            .make()

        val (isHaploid, result) = hvcf2gvcf.processCurrentHVCFVariant(vc, hapIdAndRangeToSampleMap)

        assertFalse(isHaploid)
        assertNotNull(result.first, "First gamete should resolve")
        assertNull(result.second, "Second gamete should be null when hapId not in map")
        assertEquals("LineA", result.first!!.sampleName)
    }

    /**
     * Tests the full hvcf2gvcf pipeline using a heterozygous diploid hVCF file loaded from disk.
     * Verifies that:
     * - The input hVCF file has proper diploid structure (phased genotypes with 2 alleles)
     * - Two separate gamete gVCF files are produced
     * - Both gamete gVCFs have proper headers and data records
     * - The gamete gVCFs contain variants from the correct chromosomes
     * - The gamete gVCFs have expected ASM attributes
     * - Gamete 1 (LineA) and Gamete 2 (LineB) have different variant content
     */
    @Test
    fun testHeterozygousDiploidFromFile() {
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/smallseq/Ref.fa"
        val sampleName = "DiploidHet_LineA_LineB"

        val fixtureFile = File("$diploidFixtureDir/$sampleName.h.vcf")
        assertTrue(fixtureFile.exists(), "Heterozygous diploid hVCF fixture should exist")

        verifyDiploidHvcfStructure(fixtureFile, sampleName, isHeterozygous = true)

        val testDir = "${TestExtension.testVCFDir}/testHetDiploidFromFile"
        File(testDir).deleteRecursively()
        File(testDir).mkdirs()
        fixtureFile.copyTo(File("$testDir/$sampleName.h.vcf"))

        Hvcf2Gvcf().test(
            "--db-path $dbPath --hvcf-dir $testDir --output-dir $testDir --reference-file $refFasta"
        )

        val gamete1File = File("$testDir/${sampleName}_1.g.vcf")
        val gamete2File = File("$testDir/${sampleName}_2.g.vcf")
        assertTrue(gamete1File.exists(), "Gamete 1 gVCF should exist")
        assertTrue(gamete2File.exists(), "Gamete 2 gVCF should exist")

        val gamete1Variants = readGvcfVariants(gamete1File)
        val gamete2Variants = readGvcfVariants(gamete2File)

        assertTrue(gamete1Variants.isNotEmpty(), "Gamete 1 should have variant records")
        assertTrue(gamete2Variants.isNotEmpty(), "Gamete 2 should have variant records")

        verifyGvcfHeaders(gamete1File, sampleName)
        verifyGvcfHeaders(gamete2File, sampleName)

        val gamete1Chroms = gamete1Variants.map { it.contig }.toSet()
        val gamete2Chroms = gamete2Variants.map { it.contig }.toSet()
        assertTrue(gamete1Chroms.contains("1"), "Gamete 1 should have chr1 variants")
        assertTrue(gamete1Chroms.contains("2"), "Gamete 1 should have chr2 variants")
        assertTrue(gamete2Chroms.contains("1"), "Gamete 2 should have chr1 variants")
        assertTrue(gamete2Chroms.contains("2"), "Gamete 2 should have chr2 variants")

        verifyAsmAttributes(gamete1Variants)
        verifyAsmAttributes(gamete2Variants)

        val gamete1Positions = gamete1Variants.map { "${it.contig}:${it.start}" }.toSet()
        val gamete2Positions = gamete2Variants.map { "${it.contig}:${it.start}" }.toSet()
        assertTrue(
            gamete1Positions != gamete2Positions,
            "Heterozygous gametes should have different variant positions"
        )
    }

    /**
     * Tests the full hvcf2gvcf pipeline using a homozygous diploid hVCF file loaded from disk.
     * Both gametes come from LineA, so the output gVCFs should have identical variant counts
     * and matching variant positions.
     */
    @Test
    fun testHomozygousDiploidFromFile() {
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/smallseq/Ref.fa"
        val sampleName = "DiploidHom_LineA"

        val fixtureFile = File("$diploidFixtureDir/$sampleName.h.vcf")
        assertTrue(fixtureFile.exists(), "Homozygous diploid hVCF fixture should exist")

        verifyDiploidHvcfStructure(fixtureFile, sampleName, isHeterozygous = false)

        val testDir = "${TestExtension.testVCFDir}/testHomDiploidFromFile"
        File(testDir).deleteRecursively()
        File(testDir).mkdirs()
        fixtureFile.copyTo(File("$testDir/$sampleName.h.vcf"))

        Hvcf2Gvcf().test(
            "--db-path $dbPath --hvcf-dir $testDir --output-dir $testDir --reference-file $refFasta"
        )

        val gamete1File = File("$testDir/${sampleName}_1.g.vcf")
        val gamete2File = File("$testDir/${sampleName}_2.g.vcf")
        assertTrue(gamete1File.exists(), "Gamete 1 gVCF should exist")
        assertTrue(gamete2File.exists(), "Gamete 2 gVCF should exist")

        val gamete1Variants = readGvcfVariants(gamete1File)
        val gamete2Variants = readGvcfVariants(gamete2File)

        assertTrue(gamete1Variants.isNotEmpty(), "Gamete 1 should have variant records")
        assertTrue(gamete2Variants.isNotEmpty(), "Gamete 2 should have variant records")
        assertEquals(
            gamete1Variants.size, gamete2Variants.size,
            "Homozygous gametes should have the same number of variant records"
        )

        verifyGvcfHeaders(gamete1File, sampleName)
        verifyGvcfHeaders(gamete2File, sampleName)
        verifyAsmAttributes(gamete1Variants)
        verifyAsmAttributes(gamete2Variants)

        val gamete1Positions = gamete1Variants.map { "${it.contig}:${it.start}-${it.end}" }
        val gamete2Positions = gamete2Variants.map { "${it.contig}:${it.start}-${it.end}" }
        assertEquals(
            gamete1Positions, gamete2Positions,
            "Homozygous gametes should have identical variant positions"
        )
    }

    /**
     * Tests the full hvcf2gvcf pipeline using a diploid hVCF file where one gamete
     * is the reference (Ref). Verifies that the ref gamete gVCF is composed entirely
     * of ref-block entries (NON_REF alleles with genotype 0) while the non-ref gamete
     * contains actual variants.
     */
    @Test
    fun testDiploidWithRefGameteFromFile() {
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/smallseq/Ref.fa"
        val sampleName = "DiploidRef_Ref_LineA"

        val fixtureFile = File("$diploidFixtureDir/$sampleName.h.vcf")
        assertTrue(fixtureFile.exists(), "Ref-gamete diploid hVCF fixture should exist")

        verifyDiploidHvcfStructure(fixtureFile, sampleName, isHeterozygous = true)

        val testDir = "${TestExtension.testVCFDir}/testRefDiploidFromFile"
        File(testDir).deleteRecursively()
        File(testDir).mkdirs()
        fixtureFile.copyTo(File("$testDir/$sampleName.h.vcf"))

        Hvcf2Gvcf().test(
            "--db-path $dbPath --hvcf-dir $testDir --output-dir $testDir --reference-file $refFasta"
        )

        val gamete1File = File("$testDir/${sampleName}_1.g.vcf")
        val gamete2File = File("$testDir/${sampleName}_2.g.vcf")
        assertTrue(gamete1File.exists(), "Gamete 1 (Ref) gVCF should exist")
        assertTrue(gamete2File.exists(), "Gamete 2 (LineA) gVCF should exist")

        val gamete1Variants = readGvcfVariants(gamete1File)
        val gamete2Variants = readGvcfVariants(gamete2File)

        assertTrue(gamete1Variants.isNotEmpty(), "Ref gamete gVCF should have data records")
        assertTrue(gamete2Variants.isNotEmpty(), "LineA gamete gVCF should have data records")

        verifyGvcfHeaders(gamete1File, sampleName)
        verifyGvcfHeaders(gamete2File, sampleName)
        verifyAsmAttributes(gamete1Variants)
        verifyAsmAttributes(gamete2Variants)

        val refGameteAllRefBlocks = gamete1Variants.all { vc ->
            vc.genotypes.first().isHomRef || vc.alternateAlleles.all { it.displayString == "<NON_REF>" }
        }
        assertTrue(refGameteAllRefBlocks, "Ref gamete gVCF should consist entirely of ref-block entries")

        val lineAHasSnps = gamete2Variants.any { vc ->
            vc.alternateAlleles.any { !it.isSymbolic && !it.isNoCall }
        }
        assertTrue(lineAHasSnps, "LineA gamete gVCF should contain actual variant calls (SNPs)")
    }

    /**
     * Verify that an hVCF file has proper diploid structure:
     * phased genotypes with exactly 2 alleles per record.
     * For heterozygous files, checks that at least some records have distinct ALT alleles.
     */
    private fun verifyDiploidHvcfStructure(hvcfFile: File, expectedSampleName: String, isHeterozygous: Boolean) {
        val reader = VCFFileReader(hvcfFile, false)
        val header = reader.fileHeader

        assertTrue(
            header.sampleNamesInOrder.contains(expectedSampleName),
            "hVCF header should contain sample name $expectedSampleName"
        )

        val altHeaderCount = header.metaDataInInputOrder.count { it.key == "ALT" }
        assertTrue(altHeaderCount > 0, "Diploid hVCF should have ALT header lines")

        var recordCount = 0
        var hetRecordCount = 0
        for (vc in reader) {
            recordCount++
            val genotype = vc.getGenotype(expectedSampleName)
            assertNotNull(genotype, "Each record should have a genotype for $expectedSampleName")
            assertEquals(2, genotype.alleles.size, "Diploid genotype should have exactly 2 alleles at ${vc.contig}:${vc.start}")

            if (vc.alternateAlleles.size >= 2) {
                hetRecordCount++
            }
        }
        reader.close()

        assertTrue(recordCount > 0, "hVCF should have at least one variant record")
        if (isHeterozygous) {
            assertTrue(
                hetRecordCount > 0,
                "Heterozygous diploid hVCF should have at least some records with 2 distinct ALT alleles"
            )
        }
    }

    private fun readGvcfVariants(gvcfFile: File): List<VariantContext> {
        val reader = VCFFileReader(gvcfFile, false)
        val variants = reader.toList()
        reader.close()
        return variants
    }

    private fun verifyGvcfHeaders(gvcfFile: File, expectedSampleName: String) {
        val lines = gvcfFile.readLines()
        assertEquals("##fileformat=VCFv4.2", lines[0], "gVCF should start with VCF 4.2 format header")

        val chromHeader = lines.first { it.startsWith("#CHROM") }
        assertTrue(
            chromHeader.endsWith(expectedSampleName),
            "gVCF #CHROM header should end with sample name $expectedSampleName"
        )

        val headerLines = lines.filter { it.startsWith("#") }
        assertTrue(headerLines.size >= 17, "gVCF should have standard header lines")
    }

    /**
     * Validate that a gamete gVCF produced from diploid splitting preserves the same
     * SNP calls as the original per-sample gVCF. Matches records by assembly
     * coordinates (ASM_Chr, ASM_Start) rather than reference position, because
     * the hVCF round-trip may re-place variants inside large indels at different
     * reference positions while preserving the assembly mapping.
     *
     * Only compares simple SNPs (single-base REF) since complex structural variants
     * (large deletions/insertions) are re-encoded differently after the round trip.
     * Requires at least [minMatchRate] (default 95%) of SNPs to match, since SNPs
     * near large indel boundaries may be re-encoded differently.
     */
    private fun verifyGameteMatchesOriginalGvcf(
        originalGvcfPath: String,
        generatedGvcfFile: File,
        minMatchRate: Double = 0.95
    ) {
        val originalVariants = readGvcfVariants(File(originalGvcfPath))
        val generatedVariants = readGvcfVariants(generatedGvcfFile)

        data class AsmKey(val chr: String, val start: String)

        val generatedByAsm = generatedVariants
            .filter { it.hasAttribute("ASM_Chr") && it.hasAttribute("ASM_Start") }
            .associateBy { AsmKey(it.getAttributeAsString("ASM_Chr", ""), it.getAttributeAsString("ASM_Start", "")) }

        val originalSnps = originalVariants.filter { vc ->
            vc.reference.length() == 1 &&
                    vc.alternateAlleles.any { !it.isSymbolic && !it.isNoCall && it.length() == 1 }
        }
        assertTrue(originalSnps.isNotEmpty(), "Original gVCF should have simple SNP records")

        var matchedSnps = 0
        var mismatchedSnps = 0
        for (original in originalSnps) {
            val asmChr = original.getAttributeAsString("ASM_Chr", "")
            val asmStart = original.getAttributeAsString("ASM_Start", "")
            val key = AsmKey(asmChr, asmStart)
            val generated = generatedByAsm[key]

            if (generated == null) {
                mismatchedSnps++
                continue
            }

            val originalAltBases = original.alternateAlleles
                .filter { !it.isSymbolic && !it.isNoCall && it.length() == 1 }
                .map { it.baseString }.toSet()
            val generatedAltBases = generated.alternateAlleles
                .filter { !it.isSymbolic && !it.isNoCall && it.length() == 1 }
                .map { it.baseString }.toSet()

            if (original.reference != generated.reference || originalAltBases != generatedAltBases) {
                mismatchedSnps++
                continue
            }

            val generatedGt = generated.getGenotype(0)
            if (!generatedGt.isNoCall) {
                val originalGt = original.getGenotype(0)
                if (originalGt.genotypeString != generatedGt.genotypeString) {
                    mismatchedSnps++
                    continue
                }

                val origAd = originalGt.ad
                val genAd = generatedGt.ad
                if (origAd != null && genAd != null && !origAd.contentEquals(genAd)) {
                    mismatchedSnps++
                    continue
                }
            }

            matchedSnps++
        }

        val total = matchedSnps + mismatchedSnps
        val matchRate = matchedSnps.toDouble() / total
        println("Round-trip validation: matched $matchedSnps/$total SNPs (${String.format("%.1f", matchRate * 100)}%) " +
                "from ${File(originalGvcfPath).name} against ${generatedGvcfFile.name}")

        assertTrue(matchRate >= minMatchRate,
            "SNP match rate ${String.format("%.1f", matchRate * 100)}% is below ${minMatchRate * 100}% threshold. " +
                    "Matched $matchedSnps/$total, missing $mismatchedSnps SNPs (near indel boundaries).")
    }

    /**
     * Verifies that SNPs in the generated gVCF are present in the donor gVCF.
     * Unlike [verifyGameteMatchesOriginalGvcf] which checks what fraction of donor
     * SNPs appear in the output, this checks the forward direction: what fraction of
     * output SNPs are found in the donor. This is appropriate when the output only
     * covers a subset of reference ranges (e.g., assembly hVCF round-trips) so that
     * not all donor SNPs are expected in the output.
     */
    private fun verifyOutputSnpsFoundInDonor(
        donorGvcfPath: String,
        generatedGvcfFile: File,
        minMatchRate: Double = 0.95
    ) {
        val donorVariants = readGvcfVariants(File(donorGvcfPath))
        val generatedVariants = readGvcfVariants(generatedGvcfFile)

        data class AsmKey(val chr: String, val start: String)

        val donorByAsm = donorVariants
            .filter { it.hasAttribute("ASM_Chr") && it.hasAttribute("ASM_Start") }
            .associateBy { AsmKey(it.getAttributeAsString("ASM_Chr", ""), it.getAttributeAsString("ASM_Start", "")) }

        val generatedSnps = generatedVariants.filter { vc ->
            vc.reference.length() == 1 &&
                    vc.alternateAlleles.any { !it.isSymbolic && !it.isNoCall && it.length() == 1 }
        }
        assertTrue(generatedSnps.isNotEmpty(), "Generated gVCF should have simple SNP records")

        var matchedSnps = 0
        var mismatchedSnps = 0
        for (generated in generatedSnps) {
            val asmChr = generated.getAttributeAsString("ASM_Chr", "")
            val asmStart = generated.getAttributeAsString("ASM_Start", "")
            val donor = donorByAsm[AsmKey(asmChr, asmStart)]

            if (donor == null) {
                mismatchedSnps++
                continue
            }

            val generatedAltBases = generated.alternateAlleles
                .filter { !it.isSymbolic && !it.isNoCall && it.length() == 1 }
                .map { it.baseString }.toSet()
            val donorAltBases = donor.alternateAlleles
                .filter { !it.isSymbolic && !it.isNoCall && it.length() == 1 }
                .map { it.baseString }.toSet()

            if (generated.reference != donor.reference || generatedAltBases != donorAltBases) {
                mismatchedSnps++
                continue
            }

            matchedSnps++
        }

        val total = matchedSnps + mismatchedSnps
        val matchRate = matchedSnps.toDouble() / total
        println("Forward validation: matched $matchedSnps/$total output SNPs (${String.format("%.1f", matchRate * 100)}%) " +
                "from ${generatedGvcfFile.name} against donor ${File(donorGvcfPath).name}")

        assertTrue(matchRate >= minMatchRate,
            "Output SNP match rate ${String.format("%.1f", matchRate * 100)}% is below ${minMatchRate * 100}% threshold. " +
                    "Matched $matchedSnps/$total, $mismatchedSnps output SNPs not found in donor.")
    }

    /**
     * Verifies that files with the .hvcf extension (as opposed to .h.vcf) are
     * correctly handled by buildGvcfFromHvcf: the sample name is extracted via
     * substringBeforeLast(".hvcf") and the pipeline produces the expected gVCF
     * output. Covers lines 129-130 in Hvcf2Gvcf.kt.
     */
    @Test
    fun testBuildGvcfFromHvcf_hvcfExtension() {
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/smallseq/Ref.fa"

        val testDir = "${TestExtension.testVCFDir}/testHvcfExtension"
        File(testDir).deleteRecursively()
        File(testDir).mkdirs()

        val sourceHvcf = "$diploidFixtureDir/DiploidHom_LineA.h.vcf"
        val destHvcf = "$testDir/TestHvcfExt.hvcf"
        File(sourceHvcf).copyTo(File(destHvcf))

        val hvcf2gvcf = Hvcf2Gvcf()
        hvcf2gvcf.test("--db-path $dbPath --hvcf-dir $testDir --output-dir $testDir --reference-file $refFasta")

        // "TestHvcfExt.hvcf".substringBeforeLast(".hvcf") == "TestHvcfExt"
        // Diploid input produces two gamete output files
        val gamete1 = File("$testDir/TestHvcfExt_1.g.vcf")
        val gamete2 = File("$testDir/TestHvcfExt_2.g.vcf")
        assertTrue(gamete1.exists(), "Gamete 1 gVCF should be created from .hvcf extension file")
        assertTrue(gamete2.exists(), "Gamete 2 gVCF should be created from .hvcf extension file")

        val g1Variants = readGvcfVariants(gamete1)
        val g2Variants = readGvcfVariants(gamete2)
        assertTrue(g1Variants.isNotEmpty(), "Gamete 1 should have variant records")
        assertTrue(g2Variants.isNotEmpty(), "Gamete 2 should have variant records")

        File(testDir).deleteRecursively()
    }

    /**
     * Same as above but for .hvcf.gz extension. Copies an existing bgzipped
     * .h.vcf.gz file and renames it to .hvcf.gz to verify the hvcf.gz branch.
     */
    @Test
    fun testBuildGvcfFromHvcf_hvcfGzExtension() {
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/smallseq/Ref.fa"

        val testDir = "${TestExtension.testVCFDir}/testHvcfGzExtension"
        File(testDir).deleteRecursively()
        File(testDir).mkdirs()

        val sourceHvcf = "${TestExtension.testVCFDir}/LineB.h.vcf.gz"
        val destHvcf = "$testDir/LineBGzExt.hvcf.gz"
        File(sourceHvcf).copyTo(File(destHvcf))

        val hvcf2gvcf = Hvcf2Gvcf()
        hvcf2gvcf.test("--db-path $dbPath --hvcf-dir $testDir --output-dir $testDir --reference-file $refFasta")

        // "LineBGzExt.hvcf.gz".substringBeforeLast(".hvcf") == "LineBGzExt"
        // LineB is haploid so a single gVCF file is produced
        val outputFile = File("$testDir/LineBGzExt.g.vcf")
        assertTrue(outputFile.exists(), "gVCF should be created from .hvcf.gz extension file")

        val variants = readGvcfVariants(outputFile)
        assertTrue(variants.isNotEmpty(), "gVCF from .hvcf.gz should have variant records")

        File(testDir).deleteRecursively()
    }

    /**
     * Verifies that exportSamplesGvcfFiles returns false when the underlying
     * tiledbvcf export command fails (e.g. invalid dbPath). This covers the
     * exception path at lines 558 and 615 in Hvcf2Gvcf.kt where the export
     * process returns a non-zero exit code and the error is caught.
     */
    @Test
    fun testExportSamplesGvcfFiles_returnsFalseOnBadDbPath() {
        val badDbPath = "${TestExtension.tempDir}/nonexistent_db_for_export"
        val outputDir = "${TestExtension.tempDir}/exportTestBadDb"
        File(outputDir).mkdirs()

        val result = runBlocking {
            Hvcf2Gvcf().exportSamplesGvcfFiles(
                setOf("LineA"),
                outputDir,
                badDbPath,
                "",
                1
            )
        }
        assertFalse(result, "exportSamplesGvcfFiles should return false when tiledbvcf export fails")

        File(outputDir).deleteRecursively()
    }

    /**
     * Verifies that the hvcf2gvcf pipeline silently ignores files that do not
     * have a recognized hvcf extension (.h.vcf, .h.vcf.gz, .hvcf, .hvcf.gz).
     * This tests the filter logic that guards the unreachable else-error at
     * line 132 in Hvcf2Gvcf.kt.
     */
    @Test
    fun testHvcf2Gvcf_nonHvcfFilesIgnored() {
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/smallseq/Ref.fa"

        val testDir = "${TestExtension.testVCFDir}/testNonHvcfFilesIgnored"
        File(testDir).deleteRecursively()
        File(testDir).mkdirs()

        File("$testDir/somefile.txt").writeText("not an hvcf file")
        File("$testDir/data.vcf").writeText("regular vcf, not hvcf")
        File("$testDir/data.csv").writeText("comma separated, not hvcf")

        val hvcf2gvcf = Hvcf2Gvcf()
        hvcf2gvcf.test("--db-path $dbPath --hvcf-dir $testDir --output-dir $testDir --reference-file $refFasta")

        val pipelineGvcfFiles = File(testDir).listFiles()?.filter { it.name.endsWith(".g.vcf") }
        assertTrue(pipelineGvcfFiles.isNullOrEmpty(), "No gvcf files should be produced from non-hvcf input files")

        File(testDir).deleteRecursively()
    }

    /**
     * Verifies that extractHVCFVariantsFromSample skips variant records that
     * have no alternate alleles (ALT="."). This covers the skip logic at
     * lines 340-342 in Hvcf2Gvcf.kt.
     */
    @Test
    fun testExtractHVCFVariants_skipsRecordsWithMissingAlts() {
        val testDir = "${TestExtension.testVCFDir}/testSkipMissingAlts"
        File(testDir).deleteRecursively()
        File(testDir).mkdirs()

        val vcfContent = listOf(
            "##fileformat=VCFv4.2",
            "##ALT=<ID=hapA_checksum123,Description=\"haplotype data\">",
            "##contig=<ID=1,length=100000>",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTestSample",
            "1\t1001\t.\tA\t.\t.\t.\tEND=5500\tGT\t.",
            "1\t6001\t.\tA\t<hapA_checksum123>\t.\t.\tEND=10000\tGT\t1"
        ).joinToString("\n")

        val vcfFile = File("$testDir/test_skip.h.vcf")
        vcfFile.writeText(vcfContent)

        val reader = VCFFileReader(vcfFile, false)
        val hapIdAndRangeToSampleMap = mapOf(
            Pair(ReferenceRange("1", 6001, 10000), "hapA_checksum123") to
                    listOf(HvcfVariant(ReferenceRange("1", 6001, 10000), "LineA", "hapA_checksum123"))
        )

        val (isHaploid, results) = Hvcf2Gvcf().extractHVCFVariantsFromSample(reader, hapIdAndRangeToSampleMap)
        reader.close()

        assertEquals(1, results.size, "Only the record with a valid ALT allele should be processed")
        assertTrue(isHaploid, "Single-allele genotypes should be treated as haploid")

        File(testDir).deleteRecursively()
    }

    /**
     * Verifies that buildRefBlockVariantInfo throws a NullPointerException
     * when the chromosome is not present in the reference sequence map.
     * This covers the !! operator at line 522 in Hvcf2Gvcf.kt.
     */
    @Test
    fun testBuildRefBlockVariantInfo_throwsOnMissingChrom() {
        val refSeq = mapOf("1" to NucSeq("ACGTACGTACGT"))

        assertThrows<NullPointerException> {
            Hvcf2Gvcf().buildRefBlockVariantInfo(
                refSeq,
                "nonexistent_chrom",
                Pair(1, 10),
                "nonexistent_chrom",
                Pair(1, 10),
                "+"
            )
        }
    }

    /**
     * Verifies that processCurrentHVCFVariant returns null for both gametes
     * when neither hapId is found in the hapIdAndRangeToSampleMap. This covers
     * the map lookup at lines 383-384 in Hvcf2Gvcf.kt.
     */
    @Test
    fun testProcessCurrentHVCFVariant_bothHapIdsMissingFromMap() {
        val hvcf2gvcf = Hvcf2Gvcf()
        val hapIdA = "hapA_unknown"
        val hapIdB = "hapB_unknown"

        val hapIdAndRangeToSampleMap = emptyMap<Pair<ReferenceRange, String>, List<HvcfVariant>>()

        val refAllele = Allele.create("A", true)
        val alt1 = Allele.create("<$hapIdA>", false)
        val alt2 = Allele.create("<$hapIdB>", false)
        val genotype = GenotypeBuilder("TestSample", listOf(alt1, alt2)).phased(true).make()
        val vc = VariantContextBuilder()
            .chr("1")
            .start(1001L)
            .stop(5500L)
            .alleles(listOf(refAllele, alt1, alt2))
            .genotypes(genotype)
            .attribute("END", 5500)
            .make()

        val (isHaploid, result) = hvcf2gvcf.processCurrentHVCFVariant(vc, hapIdAndRangeToSampleMap)

        assertFalse(isHaploid, "Diploid variant should not be haploid")
        assertNull(result.first, "First gamete should be null when hapId not in map")
        assertNull(result.second, "Second gamete should be null when hapId not in map")
    }

    /**
     * Tests the hvcf2gvcf pipeline using pre-built test data files. Uses an imputed hVCF
     * (LineAB.h.vcf) that references haplotypes from both LineA and LineB, along with
     * pre-built gVCF donor files. Verifies that the output gVCF is correctly composed
     * from variants of both donors across both chromosomes.
     *
     * LineAB.h.vcf composition:
     *   Chr1: ranges 1-6500 from LineA, 6501-16500 from LineB
     *   Chr2: ranges 1-6500 from LineB, 6501-12000 from LineA
     */
    @Test
    fun testImputedHvcf2GvcfWithPrebuiltData() {
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/smallseq/Ref.fa"

        val testDir = "${TestExtension.testVCFDir}/testImputedPrebuilt"
        File(testDir).deleteRecursively()
        File(testDir).mkdirs()

        File("$prebuiltGvcfDir/LineA.g.vcf").copyTo(File("$testDir/LineA.g.vcf"))
        File("$prebuiltGvcfDir/LineB.g.vcf").copyTo(File("$testDir/LineB.g.vcf"))
        File("$prebuiltImputeHvcfDir/LineAB.h.vcf").copyTo(File("$testDir/LineAB.h.vcf"))

        val hvcf2gvcf = Hvcf2Gvcf()
        hvcf2gvcf.test(
            "--db-path $dbPath --hvcf-dir $testDir --output-dir $testDir --reference-file $refFasta"
        )

        val outputFile = File("$testDir/LineAB.g.vcf")
        assertTrue(outputFile.exists(), "LineAB.g.vcf should be created from imputed hVCF")

        val variants = readGvcfVariants(outputFile)
        assertTrue(variants.isNotEmpty(), "Output gVCF should have variant records")

        verifyGvcfHeaders(outputFile, "LineAB")
        verifyAsmAttributes(variants)

        val chroms = variants.map { it.contig }.toSet()
        assertTrue(chroms.contains("1"), "Output should have chromosome 1 records")
        assertTrue(chroms.contains("2"), "Output should have chromosome 2 records")

        val chr1Variants = variants.filter { it.contig == "1" }
        val chr2Variants = variants.filter { it.contig == "2" }
        assertTrue(chr1Variants.isNotEmpty(), "Should have chr1 variants")
        assertTrue(chr2Variants.isNotEmpty(), "Should have chr2 variants")

        val chr1HasSnps = chr1Variants.any { vc ->
            vc.alternateAlleles.any { !it.isSymbolic && !it.isNoCall }
        }
        val chr2HasSnps = chr2Variants.any { vc ->
            vc.alternateAlleles.any { !it.isSymbolic && !it.isNoCall }
        }
        assertTrue(chr1HasSnps, "Chr1 output should contain SNP calls from donors")
        assertTrue(chr2HasSnps, "Chr2 output should contain SNP calls from donors")
    }

    /**
     * Tests the hvcf2gvcf pipeline using pre-built assembly hVCFs and gVCFs.
     * Runs hvcf2gvcf on renamed copies of the assembly hVCFs, using the pre-built
     * gVCFs as donor data, and verifies the output matches the original gVCFs
     * via round-trip SNP comparison.
     *
     * Assembly hVCFs are copied with alternate names (LineA_prebuilt, LineB_prebuilt)
     * to avoid output filename conflicts with the donor gVCFs.
     */
    @Test
    fun testAsmHvcf2GvcfMatchesPrebuiltGvcf() {
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/smallseq/Ref.fa"

        val testDir = "${TestExtension.testVCFDir}/testAsmPrebuiltGvcf"
        File(testDir).deleteRecursively()
        File(testDir).mkdirs()

        File("$prebuiltAsmHvcfDir/LineA.h.vcf").copyTo(File("$testDir/LineA_prebuilt.h.vcf"))
        File("$prebuiltAsmHvcfDir/LineB.h.vcf").copyTo(File("$testDir/LineB_prebuilt.h.vcf"))

        File("$prebuiltGvcfDir/LineA.g.vcf").copyTo(File("$testDir/LineA.g.vcf"))
        File("$prebuiltGvcfDir/LineB.g.vcf").copyTo(File("$testDir/LineB.g.vcf"))

        val hvcf2gvcf = Hvcf2Gvcf()
        hvcf2gvcf.test(
            "--db-path $dbPath --hvcf-dir $testDir --output-dir $testDir --reference-file $refFasta"
        )

        val lineAOutput = File("$testDir/LineA_prebuilt.g.vcf")
        val lineBOutput = File("$testDir/LineB_prebuilt.g.vcf")
        assertTrue(lineAOutput.exists(), "LineA_prebuilt.g.vcf should be created")
        assertTrue(lineBOutput.exists(), "LineB_prebuilt.g.vcf should be created")

        val lineAVariants = readGvcfVariants(lineAOutput)
        val lineBVariants = readGvcfVariants(lineBOutput)
        assertTrue(lineAVariants.isNotEmpty(), "LineA output should have variants")
        assertTrue(lineBVariants.isNotEmpty(), "LineB output should have variants")

        verifyGvcfHeaders(lineAOutput, "LineA_prebuilt")
        verifyGvcfHeaders(lineBOutput, "LineB_prebuilt")
        verifyAsmAttributes(lineAVariants)
        verifyAsmAttributes(lineBVariants)

        val lineAChroms = lineAVariants.map { it.contig }.toSet()
        val lineBChroms = lineBVariants.map { it.contig }.toSet()
        assertTrue(lineAChroms.contains("1") && lineAChroms.contains("2"),
            "LineA output should have variants on both chromosomes")
        assertTrue(lineBChroms.contains("1") && lineBChroms.contains("2"),
            "LineB output should have variants on both chromosomes")

        verifyOutputSnpsFoundInDonor("$prebuiltGvcfDir/LineA.g.vcf", lineAOutput)
        verifyOutputSnpsFoundInDonor("$prebuiltGvcfDir/LineB.g.vcf", lineBOutput)
    }

    private fun verifyAsmAttributes(variants: List<VariantContext>) {
        for (vc in variants) {
            assertNotNull(
                vc.getAttribute("ASM_Chr"),
                "Variant at ${vc.contig}:${vc.start} should have ASM_Chr attribute"
            )
            assertNotNull(
                vc.getAttribute("ASM_Start"),
                "Variant at ${vc.contig}:${vc.start} should have ASM_Start attribute"
            )
            assertNotNull(
                vc.getAttribute("ASM_End"),
                "Variant at ${vc.contig}:${vc.start} should have ASM_End attribute"
            )
            assertNotNull(
                vc.getAttribute("ASM_Strand"),
                "Variant at ${vc.contig}:${vc.start} should have ASM_Strand attribute"
            )

            val asmStrand = vc.getAttributeAsString("ASM_Strand", "")
            assertTrue(
                asmStrand == "+" || asmStrand == "-",
                "ASM_Strand should be + or - at ${vc.contig}:${vc.start}"
            )
        }
    }
}