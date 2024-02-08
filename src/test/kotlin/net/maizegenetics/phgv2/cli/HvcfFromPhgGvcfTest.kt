package net.maizegenetics.phgv2.cli

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.testing.test
import com.google.common.collect.Range
import com.google.common.collect.RangeSet
import com.google.common.collect.TreeRangeSet
import net.maizegenetics.phgv2.utils.Position
import org.junit.jupiter.api.*
import java.io.File
import java.nio.file.Paths
import kotlin.test.assertEquals

class HvcfFromPhgGvcfTest {
    companion object {

        @JvmStatic
        @BeforeAll
        //@BeforeEach
        fun setup() {
            File(TestExtension.testVCFDir).mkdirs()
            File(TestExtension.testTileDBURI).mkdirs()
        }

        @JvmStatic
        @AfterAll
        //@AfterEach
        fun tearDown() {
            File(TestExtension.testVCFDir).deleteRecursively()
            File(TestExtension.testTileDBURI).deleteRecursively()
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testCliktParams() {
        val hvcfFromPhgGvcf = HvcfFromPhgGvcf()

        // There are only 3 required parameters - test for missing each one
        val resultMissingBed =
            hvcfFromPhgGvcf.test("--db-path ${TestExtension.testTileDBURI} --gvcf-dir ${TestExtension.testVCFDir} --reference-file ${TestExtension.testRefFasta} ")
        assertEquals(resultMissingBed.statusCode, 1)
        assertEquals(
            "Usage: hvcf-from-phg-gvcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --bed: --bed must not be blank\n", resultMissingBed.output
        )
        val resultMissingRef =
            hvcfFromPhgGvcf.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} --gvcf-dir ${TestExtension.testMafDir}")
        assertEquals(resultMissingRef.statusCode, 1)
        assertEquals(
            "Usage: hvcf-from-phg-gvcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --reference-file: --reference-file must not be blank\n", resultMissingRef.output
        )


        val resultMissingGvcfDir =
            hvcfFromPhgGvcf.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} --reference-file ${TestExtension.testRefFasta}")
        assertEquals(resultMissingGvcfDir.statusCode, 1)
        assertEquals(
            "Usage: hvcf-from-phg-gvcf [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --gvcf-dir: --gvcf-dir must not be blank\n", resultMissingGvcfDir.output
        )

    }

    @Test
    fun testHvcfFromGvcf() {
        // Copy the gvcf files from data/test/smallseq to the testVCFDir
        val gvcfDir = TestExtension.testVCFDir
        val gvcfFile = File(TestExtension.smallSeqLineAGvcfFile)
        gvcfFile.copyTo(File(gvcfDir, gvcfFile.name))

        //Need to create the agc record before we run this:
        // We will not be loading to tiledb, but we will be pulling sequence from AGC.
        // AGC lives in the same directory as the tiledb datasets.  Need all the fasta files in the agc record
        // for each sample in the gvcf files

        // copy the fasta files
        val fastaInputDir = "data/test/smallseq"
        val fastaOutputDir = TestExtension.testOutputFastaDir
        val fastaFiles = File(fastaInputDir).listFiles { file -> file.extension == "fa" }
        fastaFiles.forEach { file -> file.copyTo(File(fastaOutputDir, file.name) ) }

        val dbPath = TestExtension.testTileDBURI
        val refFasta = TestExtension.smallseqRefFile
        println("refFasta: $refFasta")

        // get the full path fasta file names  from the fastaInput, write to fileList
        val fileList = mutableListOf<String>()
        File(fastaOutputDir).walk(FileWalkDirection.TOP_DOWN).filter{it.name.endsWith(".fa")}.forEach {
            fileList.add(it.toString())
        }
        // write the full path fasta file names  from the fastaOutputDir to a single file, one per line, in tempDir
        // This file will be used as input to the agc-compress command
        val fastaCreateFileNamesFile = File(dbPath, "agcFileList.txt")
        fastaCreateFileNamesFile.writeText(fileList.joinToString("\n"))

        val agcCompress = AgcCompress()
        // Create the compressed file
        val agcResult =
            agcCompress.test("--fasta-list ${fastaCreateFileNamesFile} --db-path ${dbPath} --reference-file ${refFasta}")
        println(agcResult.output)

        val bedFile = TestExtension.smallseqAnchorsBedFile
        val hvcfFromGvcf = HvcfFromPhgGvcf()
        val result =
            hvcfFromGvcf.test("--db-path ${dbPath} --bed ${bedFile} --reference-file ${refFasta} --gvcf-dir ${gvcfDir} ")
        println(result.output)

        // read the created hvcf file, which lives in the gvcfDir and has extension .h.vcf.gz
        val hvcfFile = File(gvcfDir).listFiles { _, name -> name.endsWith(".h.vcf.gz") }[0].toString()
        val hvcfLines = bufferedReader(hvcfFile).readLines()
        println("hvcfLines.size: ${hvcfLines.size}")

        var altLines = 0
        var variantLines = 0
        for (line in hvcfLines) {
            if (line.startsWith("#")) {
                if (line.startsWith("##ALT")) altLines++
            } else {
                variantLines++ // why isn't this count 33 ??
            }
        }
        println("altLines: $altLines, variantLines: $variantLines")

        // read number of lines in the bedFile
        val bedLines = bufferedReader(bedFile).readLines()

        // The bed file has 40 anchors
        // The gvcf file may not have data that represents all 40 anchors
        // Can do these verifications:
        // 1.  Verify the number of variant lines matches the number of ALT header lines in the hvcf file
        //     This is adjusted when there are hash collisions
        // 2.  Manually check the gvcf to see how many ref ranges are represented.  This should match
        //     the number of ALT header lines in the hvcf file (and the number of variant lines)

        val refRangesMissing = findMissingRefRanges(gvcfFile,File(bedFile))
        val expectedVariantLines = bedLines.size - refRangesMissing
        // the number of variantLines should match the expectedVariantLines
        assertEquals(expectedVariantLines,variantLines)

        // In this gvcf file, there is 1 hash collision. There will only be 1 ALT line per seq hash, so the number of ALT lines
        // should be 1 less than the number of variant lines
        assertEquals(variantLines-1,altLines)

    }

    @Test
    fun compareToCreateMafVcfOutput() {
        // This test will compare the output of HvcfFromPhgGvcf to the hvcf output of CreateMafVcf
        // the gvcf used to create the hvcf via the new HvcfFromPhgGvcf code is the gvcf
        // created by CreateMafVcf.

        // Setup and run CreateMafVcf to get the gvcf file for testing
        val fastaCreateFileNamesFile = "data/test/buildMAFVCF/fastaCreateFileNames.txt"
        val dbPath = TestExtension.testTileDBURI
        val refFasta = "data/test/buildMAFVCF/B73_Test.fa"

        val agcCompress = AgcCompress()
        // Create the initial compressed file
        val agcResult = agcCompress.test("--fasta-list ${fastaCreateFileNamesFile} --db-path ${dbPath} --reference-file ${refFasta}")
        println(agcResult.output)

        val createMAFVCF = CreateMafVcf()
        val result = createMAFVCF.test("--db-path ${dbPath} --bed data/test/buildMAFVCF/B73_Test.bed --reference-file ${refFasta} --maf-dir data/test/buildMAFVCF/mafs/ -o ${TestExtension.testVCFDir}")
        println(result.output)

        // Copy the gvcf files created by CreateMafVcf to a new folder
        // Run the HvcfFromPhgGvcf code on the gvcf files in the new folder
        // Verify the h.vcf files are identical to the hvcf file created by CreateMafVcf
        val newHvcfDir = "${TestExtension.tempDir}newHvcf/"
        File(newHvcfDir).mkdirs()
        val gvcfFiles = File(TestExtension.testVCFDir).listFiles { _, name -> name.endsWith(".g.vcf.gz") }
        gvcfFiles.forEach { file -> file.copyTo(File(newHvcfDir, file.name)) }

        // Run the HvcfFromPhgGvcf code on the gvcf files in the testVCFDir
        val hvcfFromGvcf = HvcfFromPhgGvcf()
        val hvcfResult = hvcfFromGvcf.test("--db-path ${dbPath} --bed data/test/buildMAFVCF/B73_Test.bed --reference-file ${refFasta} --gvcf-dir ${newHvcfDir} ")

        // Verify the h.vcf.gz file created by HvcfFromPhgGvcf is identical to the hvcf file created by CreateMafVcf
        // This will be done by comparing the number of lines in the files and the first 10 lines of the files
        val hvcfFile = File(newHvcfDir).listFiles { _, name -> name.endsWith(".h.vcf.gz") }[0]
        val createMafVcfHvcfFile = File(TestExtension.testVCFDir).listFiles { _, name -> name.endsWith(".h.vcf.gz") }[0]
        // Compare the contents of the two files
        val hvcfLines = hvcfFile.readLines()
        val createMafVcfHvcfLines = createMafVcfHvcfFile.readLines()
        println("hvcfLines.size: ${hvcfLines.size}, createMafVcfHvcfLines.size: ${createMafVcfHvcfLines.size}")

        assertEquals(hvcfLines.size, createMafVcfHvcfLines.size)
        for (idx in 0..hvcfLines.size-1) {
            assertEquals(hvcfLines[idx], createMafVcfHvcfLines[idx])
        }
    }


    @Test
    fun testFindMissingRefRanges() {
        val gvcfFile = File(TestExtension.smallSeqLineAGvcfFile)
        val bedFile = File(TestExtension.smallseqAnchorsBedFile)
        val missing = findMissingRefRanges(gvcfFile,bedFile)
        // There are 40 anchors in the bed file,the gvcf file has 34 variant records,
        // so there are 6 ref ranges missing from the gvcf file
        assertEquals(6,missing)
    }

    fun findMissingRefRanges(gvcfFile:File,bedFile:File):Int {
        // This test takes a gvcf file and a bed file and finds the ref ranges in the bed file that are not
        // represented in the gvcf file.  This is useful for verifying the hvcf file created by HvcfFromPhgGvcf
        // is correct.  The hvcf file should have a ref range for every range in the bed file, and the gvcf file

        val gvcfLines = gvcfFile.readLines()
        val bedLines = bedFile.readLines()
        val gvcfPositions = mutableListOf<Position>()

        // get set of interval ranges
        val intervalRanges: RangeSet<Position> = TreeRangeSet.create()
        for (line in bedLines) {
            if (line.startsWith("#")) continue
            val fields = line.split("\t")
            val interval =
                Range.closed(
                    Position(
                        fields[0],
                        fields[1].toInt() + 1
                    ), Position(fields[0], fields[2].toInt())
                )
            intervalRanges.add(interval)
        }
        println("\nintervalRanges.size: ${intervalRanges.asRanges().size}")
        for (line in gvcfLines) {
            if (line.startsWith("#")) continue
            val fields = line.split("\t")
            val position = Position(fields[0], fields[1].toInt())
            gvcfPositions.add(position)
            if (fields[7].contains(";END=")) {
                val end = fields[7].split(";END=")[1].split(";")[0].toInt()
                val start = position.position+1
                for (idx in start until end) {
                    val interval = Range.closed(Position(fields[0], idx), Position(fields[0], end))
                    intervalRanges.add(interval)
                }
            }
        }
        // Create a list of ranges from  intervalRanges that have no positions from gvcfPositions
        val missingRanges = mutableListOf<Range<Position>>()
        for (range in intervalRanges.asRanges()) {
            var found = false
            if (range.lowerEndpoint().contig == "1" && range.lowerEndpoint().position == 44001) {
                //println("range: $range")
            }
            for (position in gvcfPositions) {
                if (range.contains(position)) {
                    found = true
                    break
                }
            }
            if (!found) missingRanges.add(range)
        }
        println("missingRanges.size: ${missingRanges.size}")
        // print ranges that are missing from the gvcf file
        println("Ranges missing from gvcf file:")
        for (range in missingRanges) {
            println(range)
        }
        return missingRanges.size
    }
}