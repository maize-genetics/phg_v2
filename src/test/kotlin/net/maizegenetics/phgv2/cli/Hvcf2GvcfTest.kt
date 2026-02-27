package net.maizegenetics.phgv2.cli

import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.testing.test
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContextBuilder
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

    @Test
    fun testDiploidHvcf2Gvcf() {
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/smallseq/Ref.fa"
        val refSeq = CreateMafVcf().buildRefGenomeSeq(refFasta)

        val testDir = "${TestExtension.testVCFDir}/testDiploidDir"
        File(testDir).deleteRecursively()
        File(testDir).mkdirs()

        val hvcfFiles = File(TestExtension.testVCFDir).listFiles()!!
            .filter { it.name.endsWith(".h.vcf") || it.name.endsWith(".h.vcf.gz") }
            .map { it.path }
        val graph = HaplotypeGraph(hvcfFiles)

        buildDiploidHvcf(
            graph, "TestDiploid",
            SampleGamete("LineA"), SampleGamete("LineB"),
            refSeq, "$testDir/TestDiploid.h.vcf"
        )

        Hvcf2Gvcf().test(
            "--db-path $dbPath --hvcf-dir $testDir --output-dir $testDir --reference-file $refFasta"
        )

        assertTrue(File("$testDir/TestDiploid_1.g.vcf").exists(), "Gamete 1 gVCF should exist")
        assertTrue(File("$testDir/TestDiploid_2.g.vcf").exists(), "Gamete 2 gVCF should exist")

        val lines1 = File("$testDir/TestDiploid_1.g.vcf").readLines()
        val lines2 = File("$testDir/TestDiploid_2.g.vcf").readLines()

        assertEquals("##fileformat=VCFv4.2", lines1[0])
        assertEquals("##fileformat=VCFv4.2", lines2[0])
        val headerLine1 = lines1.first { it.startsWith("#CHROM") }
        val headerLine2 = lines2.first { it.startsWith("#CHROM") }
        assertTrue(headerLine1.endsWith("TestDiploid"), "Gamete 1 gVCF should have correct sample name")
        assertTrue(headerLine2.endsWith("TestDiploid"), "Gamete 2 gVCF should have correct sample name")

        val dataLines1 = lines1.filter { !it.startsWith("#") }
        val dataLines2 = lines2.filter { !it.startsWith("#") }
        assertTrue(dataLines1.isNotEmpty(), "Gamete 1 gVCF should have data records")
        assertTrue(dataLines2.isNotEmpty(), "Gamete 2 gVCF should have data records")
    }

    @Test
    fun testDiploidHvcf2Gvcf_homozygous() {
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/smallseq/Ref.fa"
        val refSeq = CreateMafVcf().buildRefGenomeSeq(refFasta)

        val testDir = "${TestExtension.testVCFDir}/testDiploidHomozygousDir"
        File(testDir).deleteRecursively()
        File(testDir).mkdirs()

        val hvcfFiles = File(TestExtension.testVCFDir).listFiles()!!
            .filter { it.name.endsWith(".h.vcf") || it.name.endsWith(".h.vcf.gz") }
            .map { it.path }
        val graph = HaplotypeGraph(hvcfFiles)

        buildDiploidHvcf(
            graph, "TestHomozygous",
            SampleGamete("LineA"), SampleGamete("LineA"),
            refSeq, "$testDir/TestHomozygous.h.vcf"
        )

        Hvcf2Gvcf().test(
            "--db-path $dbPath --hvcf-dir $testDir --output-dir $testDir --reference-file $refFasta"
        )

        assertTrue(File("$testDir/TestHomozygous_1.g.vcf").exists(), "Gamete 1 gVCF should exist")
        assertTrue(File("$testDir/TestHomozygous_2.g.vcf").exists(), "Gamete 2 gVCF should exist")

        val dataLines1 = File("$testDir/TestHomozygous_1.g.vcf").readLines().filter { !it.startsWith("#") }
        val dataLines2 = File("$testDir/TestHomozygous_2.g.vcf").readLines().filter { !it.startsWith("#") }
        assertTrue(dataLines1.isNotEmpty(), "Gamete 1 gVCF should have data records")
        assertTrue(dataLines2.isNotEmpty(), "Gamete 2 gVCF should have data records")
        assertEquals(dataLines1.size, dataLines2.size, "Homozygous diploid should have equal records per gamete")
    }

    @Test
    fun testDiploidHvcf2Gvcf_withRefGamete() {
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/smallseq/Ref.fa"
        val refSeq = CreateMafVcf().buildRefGenomeSeq(refFasta)

        val testDir = "${TestExtension.testVCFDir}/testDiploidRefDir"
        File(testDir).deleteRecursively()
        File(testDir).mkdirs()

        val hvcfFiles = mutableListOf<String>()
        File(TestExtension.testVCFDir).listFiles()!!
            .filter { it.name.endsWith(".h.vcf") || it.name.endsWith(".h.vcf.gz") }
            .forEach { hvcfFiles.add(it.path) }
        File("${dbPath}/hvcf_files").listFiles()!!
            .filter { it.name.endsWith(".h.vcf") || it.name.endsWith(".h.vcf.gz") }
            .forEach { hvcfFiles.add(it.path) }
        val graph = HaplotypeGraph(hvcfFiles.distinct())

        buildDiploidHvcf(
            graph, "TestRefGamete",
            SampleGamete("Ref"), SampleGamete("LineA"),
            refSeq, "$testDir/TestRefGamete.h.vcf"
        )

        Hvcf2Gvcf().test(
            "--db-path $dbPath --hvcf-dir $testDir --output-dir $testDir --reference-file $refFasta"
        )

        assertTrue(File("$testDir/TestRefGamete_1.g.vcf").exists(), "Gamete 1 (Ref) gVCF should exist")
        assertTrue(File("$testDir/TestRefGamete_2.g.vcf").exists(), "Gamete 2 (LineA) gVCF should exist")

        val dataLines1 = File("$testDir/TestRefGamete_1.g.vcf").readLines().filter { !it.startsWith("#") }
        val dataLines2 = File("$testDir/TestRefGamete_2.g.vcf").readLines().filter { !it.startsWith("#") }
        assertTrue(dataLines1.isNotEmpty(), "Ref gamete gVCF should have data records")
        assertTrue(dataLines2.isNotEmpty(), "LineA gamete gVCF should have data records")
    }

    private fun buildDiploidHvcf(
        graph: HaplotypeGraph,
        sampleName: String,
        gamete1: SampleGamete,
        gamete2: SampleGamete,
        refSeq: Map<String, NucSeq>,
        outputFile: String
    ) {
        val altHeaders = mutableListOf<AltHeaderMetaData>()
        val variantContexts = mutableListOf<htsjdk.variant.variantcontext.VariantContext>()

        for (range in graph.ranges()) {
            val hapId1 = graph.sampleToHapId(range, gamete1)
            val hapId2 = graph.sampleToHapId(range, gamete2)
            val refAllele = refSeq[range.contig]!!.get(range.start - 1).toString()
            val startPos = Position(range.contig, range.start)
            val endPos = Position(range.contig, range.end)

            variantContexts.add(createDiploidHVCFRecord(sampleName, startPos, endPos, listOf(hapId1, hapId2), refAllele))
            listOfNotNull(hapId1, hapId2).distinct().mapNotNull { graph.altHeader(it) }
                .forEach { altHeaders.add(it) }
        }

        val headerSet = altHeaders.distinctBy { it.id }.map { altHeaderMetadataToVCFHeaderLine(it) }.toSet()
        exportVariantContext(sampleName, variantContexts, outputFile, refSeq, headerSet)
    }
}