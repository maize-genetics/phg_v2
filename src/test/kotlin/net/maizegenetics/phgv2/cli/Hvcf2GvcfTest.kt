package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.utils.verifyURI
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals
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
            Initdb().createDataSets(dbPath,"")

            // Create the agc compressed file
            println("testSimpleHvcf2Gvcf:running agcCompress")
            val agcCompress = AgcCompress()
            var agcResult = agcCompress.test("--fasta-list ${TestExtension.smallseqAssembliesListFile} --db-path ${dbPath} --reference-file ${TestExtension.smallseqRefFile}")
            println(agcResult.output)

            println("testSimpleHvcf2Gvcf:running CreateRefVcf")
            var result = CreateRefVcf().test("--bed $ranges --reference-name $refName --reference-file $refFasta --reference-url ${refUrl} --db-path $dbPath")
            assertEquals(0, result.statusCode )

            // RUn alignAssemblies test to get MAF files.
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
            // remember - there is no gvcf for the ref
            println("testSimpleHvcf2Gvcf: running CreateMafVcf")
            val createMafVcf = CreateMafVcf()
            result = createMafVcf.test("--db-path ${dbPath} --bed ${ranges} --reference-file ${refFasta} --maf-dir ${TestExtension.tempDir} -o ${TestExtension.testVCFDir}")
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
                    "Error: invalid value for --reference-file: --reference-file must not be blank\n", resultMissingRef.output
        )


        val resultMissingHvcfDir =
            hvcf2gvcf.test("--db-path ${TestExtension.testTileDBURI}  --reference-file ${TestExtension.testRefFasta}")
        assertEquals(resultMissingHvcfDir.statusCode, 1)
        assertEquals(
            "Usage: hvcf2gvcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --hvcf-dir: --hvcf-dir must not be blank\n", resultMissingHvcfDir.output
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
        val result = hvcf2gvcf.test("--db-path ${dbPath} --hvcf-dir $testGVCFdir --output-dir ${testGVCFdir} --reference-file ${refFasta}")
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
        assertTrue(lines.contains("1\t5501\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5561;ASM_Start=5501;ASM_Strand=+;END=5561\tGT:AD:DP:PL\t0:30,0:30:0,90,90"))
        assertTrue(lines.contains("1\t6501\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=6531;ASM_Start=6501;ASM_Strand=+;END=6531\tGT:AD:DP:PL\t0:30,0:30:0,90,90"))
        assertTrue(lines.contains("2\t1001\t.\tC\t<NON_REF>\t.\t.\tASM_Chr=2;ASM_End=1064;ASM_Start=1001;ASM_Strand=+;END=1064\tGT:AD:DP:PL\t0:30,0:30:0,90,90"))
        assertTrue(lines.contains("2\t5501\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=2;ASM_End=5508;ASM_Start=5501;ASM_Strand=+;END=5508\tGT:AD:DP:PL\t0:30,0:30:0,90,90"))
        assertTrue(lines.contains("2\t6501\t.\tG\t<NON_REF>\t.\t.\tASM_Chr=2;ASM_End=6501;ASM_Start=6501;ASM_Strand=+;END=6501\tGT:AD:DP:PL\t0:30,0:30:0,90,90"))

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

        // Run hvcf2gvcf on the copied file
        val hvcf2gvcf = Hvcf2Gvcf()
        val result = hvcf2gvcf.test("--db-path ${dbPath} --hvcf-dir $testGVCFdir --output-dir ${testGVCFdir} --reference-file ${refFasta}")
        // verify the output file exists, which will be LineBPath.g.vcf
        assertTrue(File("${testGVCFdir}/TestLine2.g.vcf").exists())
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
        CreateMafVCFTest().compareTwoGVCFFiles("data/test/hvcfToGvcf/TestLine2_truth.g.vcf", "${testGVCFdir}/TestLine2.g.vcf")

        println("testPathHvcf2Gvcf done !!")
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

        val chrom1entry = "1\t1001\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=1;ASM_End=5500;ASM_Start=1001;ASM_Strand=+;END=5500\tGT:AD:DP:PL\t0:30,0:30:0,90,90"
        val chom2entry = "2\t1001\t.\tA\t<NON_REF>\t.\t.\tASM_Chr=2;ASM_End=5500;ASM_Start=1001;ASM_Strand=+;END=5500\tGT:AD:DP:PL\t0:30,0:30:0,90,90"

        // verify the file contains both the chrom1entry and the chrom2entry
        assertTrue(lines.contains(chrom1entry))
        assertTrue(lines.contains(chom2entry))
        println("ref gvcf file created: $outputDir/$refName.vcf")
    }

}