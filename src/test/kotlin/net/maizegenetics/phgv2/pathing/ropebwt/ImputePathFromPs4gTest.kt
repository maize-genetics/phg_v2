package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File

class ImputePathFromPs4gTest {
    companion object {
        val tempTestDir = "${TestExtension.tempDir}ropebwtTest/"

        @JvmStatic
        @BeforeAll
        fun setup() {
            resetDirs()
            setupDebugLogging()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            resetDirs()
        }

        private fun resetDirs() {
            File(TestExtension.tempDir).deleteRecursively()
            File(TestExtension.testOutputFastaDir).deleteRecursively()
            File(TestExtension.testOutputDir).deleteRecursively()
            File(tempTestDir).deleteRecursively()

            File(TestExtension.tempDir).mkdirs()
            File(TestExtension.testOutputFastaDir).mkdirs()
            File(TestExtension.testOutputDir).mkdirs()
            File(tempTestDir).mkdirs()
        }

        private fun createPs4gFile(outputPath: String, gametes: Map<SampleGamete, Int>, data: List<PS4GData>) {
            val counts = gametes.keys.associateWith { sg ->
                val idx = gametes[sg]!!
                data.filter { it.gameteList.contains(idx) }.sumOf { it.count }
            }
            PS4GUtils.writeOutPS4GFile(data, counts, gametes, outputPath, emptyList(), "test command")
        }

        private fun createKeyFile(keyFilePath: String, sampleName: String, ps4gPath: String) {
            File(keyFilePath).writeText("sampleName\tfilename\n$sampleName\t$ps4gPath\n")
        }
    }

    @Test
    fun testCliktParams() {
        // Missing both --path-keyfile and --read-file
        val noInput = ImputePathFromPs4g().test("--out-path-dir $tempTestDir")
        assertEquals(1, noInput.statusCode)
        assertTrue(noInput.stderr.contains("must provide one of --path-keyfile, --read-file"),
            "Expected missing-input error but got: ${noInput.stderr}")

        // Missing --out-path-dir
        val noOutputDir = ImputePathFromPs4g().test("--read-file someFile.ps4g")
        assertEquals(1, noOutputDir.statusCode)
        assertTrue(noOutputDir.stderr.contains("missing option --out-path-dir"),
            "Expected missing-output-dir error but got: ${noOutputDir.stderr}")

        // Invalid --path-type value
        val badPathType = ImputePathFromPs4g().test(
            "--read-file someFile.ps4g --out-path-dir $tempTestDir --path-type invalid"
        )
        assertEquals(1, badPathType.statusCode)
        assertTrue(badPathType.stderr.contains("path-type"),
            "Expected path-type error but got: ${badPathType.stderr}")
    }

    @Test
    fun testImputeHaploidPathOutput() {
        // All reads support gamete 0 (lineA:0) at every position.
        // The HMM should assign the entire chromosome to lineA:0, and the 5 bins should
        // merge into a single record spanning the contig.
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = (1..5).map { bin ->
            PS4GData(listOf(0), Position("chr1", bin), 10)
        }
        val ps4gFile = "$tempTestDir/singleGamete.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)

        val keyFile = "$tempTestDir/singleGameteKey.txt"
        createKeyFile(keyFile, "sampleX", ps4gFile)

        val outputDir = "$tempTestDir/haploidOut/"
        val result = ImputePathFromPs4g().test(
            "--path-keyfile $keyFile --out-path-dir $outputDir --prob-correct 0.99 --prob-same 0.9999 --bin-size 1"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val bedFile = File("$outputDir/sampleX_imputed_path.bed")
        assertTrue(bedFile.exists(), "Output BED file was not created")

        val lines = bedFile.readLines()
        assertEquals("chrom\tstart\tend\tparent1", lines[0])

        val dataLines = lines.drop(1)
        assertEquals(1, dataLines.size, "Expected the 5 identical bins to merge into 1 record")
        assertEquals(listOf("chr1", "0", "5", "lineA:0"), dataLines[0].split("\t"))
    }

    @Test
    fun testImputeHaploidPathOutputWithHighReadCount() {
        // Since the EmissionProbability caches values for readCount <= 50, need to test read counts > 50 as well
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = (1..5).map { bin ->
            PS4GData(listOf(0), Position("chr1", bin), 60)
        }
        val ps4gFile = "$tempTestDir/singleGamete.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)

        val keyFile = "$tempTestDir/singleGameteKey.txt"
        createKeyFile(keyFile, "sampleX", ps4gFile)

        val outputDir = "$tempTestDir/haploidOut/"
        val result = ImputePathFromPs4g().test(
            "--path-keyfile $keyFile --out-path-dir $outputDir --prob-correct 0.99 --prob-same 0.9999 --bin-size 1"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val bedFile = File("$outputDir/sampleX_imputed_path.bed")
        assertTrue(bedFile.exists(), "Output BED file was not created")

        val lines = bedFile.readLines()
        assertEquals("chrom\tstart\tend\tparent1", lines[0])

        val dataLines = lines.drop(1)
        assertEquals(1, dataLines.size, "Expected the 5 identical bins to merge into 1 record")
        assertEquals(listOf("chr1", "0", "5", "lineA:0"), dataLines[0].split("\t"))
    }

    @Test
    fun testImputeHaploidPathExpandBins() {
        // Same data as testImputeHaploidPathOutput, but --expand-bins should suppress merging
        // and emit one record per bin.
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = (1..5).map { bin ->
            PS4GData(listOf(0), Position("chr1", bin), 10)
        }
        val ps4gFile = "$tempTestDir/expandBins.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)

        val keyFile = "$tempTestDir/expandBinsKey.txt"
        createKeyFile(keyFile, "expandSample", ps4gFile)

        val outputDir = "$tempTestDir/expandOut/"
        val result = ImputePathFromPs4g().test(
            "--path-keyfile $keyFile --out-path-dir $outputDir --prob-correct 0.99 --prob-same 0.9999 " +
                    "--bin-size 1 --expand-bins"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val dataLines = File("$outputDir/expandSample_imputed_path.bed").readLines().drop(1)
        assertEquals(5, dataLines.size, "Expected one record per bin with --expand-bins")
        dataLines.forEach { line ->
            val cols = line.split("\t")
            assertEquals("chr1", cols[0])
            assertEquals("lineA:0", cols[3], "Expected lineA:0 but got ${cols[3]} in: $line")
        }
    }

    @Test
    fun testImputeHaploidPathBedCoordinates() {
        // Use bins 10, 20, 30, 40, 50 with binSize=100 so we can verify the midpoint arithmetic
        // exactly. Coordinates are 0-based half-open, and --expand-bins keeps every cut so the
        // arithmetic is visible:
        //   Row 0: [0,    (10+20)/2*100=1500)
        //   Row 1: [1500, (20+30)/2*100=2500)
        //   Row 2: [2500, (30+40)/2*100=3500)
        //   Row 3: [3500, (40+50)/2*100=4500)
        //   Row 4: [4500, 50*100=5000)
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val bins = listOf(10, 20, 30, 40, 50)
        val ps4gData = bins.map { bin -> PS4GData(listOf(0), Position("chr1", bin), 5) }

        val ps4gFile = "$tempTestDir/coordTest.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)

        val keyFile = "$tempTestDir/coordTestKey.txt"
        createKeyFile(keyFile, "coordSample", ps4gFile)

        val outputDir = "$tempTestDir/coordOut/"
        val result = ImputePathFromPs4g().test(
            "--path-keyfile $keyFile --out-path-dir $outputDir --prob-correct 0.99 --prob-same 0.9999 " +
                    "--bin-size 100 --expand-bins"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val lines = File("$outputDir/coordSample_imputed_path.bed").readLines().drop(1) // skip header
        assertEquals(5, lines.size, "Expected 5 data rows")

        fun cols(row: Int) = lines[row].split("\t")

        assertEquals(listOf("0", "1500"), listOf(cols(0)[1], cols(0)[2]))
        assertEquals(listOf("1500", "2500"), listOf(cols(1)[1], cols(1)[2]))
        assertEquals(listOf("2500", "3500"), listOf(cols(2)[1], cols(2)[2]))
        assertEquals(listOf("3500", "4500"), listOf(cols(3)[1], cols(3)[2]))
        assertEquals(listOf("4500", "5000"), listOf(cols(4)[1], cols(4)[2]))

        // Merging is the default, so the same input collapses to one record spanning the contig.
        val mergedDir = "$tempTestDir/coordMergedOut/"
        val mergedResult = ImputePathFromPs4g().test(
            "--path-keyfile $keyFile --out-path-dir $mergedDir --prob-correct 0.99 --prob-same 0.9999 --bin-size 100"
        )
        assertEquals(0, mergedResult.statusCode, "Command failed:\n${mergedResult.stderr}")

        val mergedLines = File("$mergedDir/coordSample_imputed_path.bed").readLines().drop(1)
        assertEquals(1, mergedLines.size, "Expected the 5 identical bins to merge into 1 record")
        assertEquals(listOf("chr1", "0", "5000", "lineA:0"), mergedLines[0].split("\t"))
    }

    @Test
    fun testImputeHaploidPathWithRecombination() {
        // First 4 bins: all reads support gamete 0 (lineA:0)
        // Last 4 bins: all reads support gamete 1 (lineB:0)
        // With very low probSame the Viterbi path should switch to lineB:0 after bin 4.
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = (1..4).map { bin ->
            PS4GData(listOf(0), Position("chr1", bin), 20)
        } + (5..8).map { bin ->
            PS4GData(listOf(1), Position("chr1", bin), 20)
        }

        val ps4gFile = "$tempTestDir/recombTest.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)

        val keyFile = "$tempTestDir/recombKey.txt"
        createKeyFile(keyFile, "recombSample", ps4gFile)

        val outputDir = "$tempTestDir/recombOut/"
        val result = ImputePathFromPs4g().test(
            "--path-keyfile $keyFile --out-path-dir $outputDir --prob-correct 0.99 --prob-same 0.5 --bin-size 1"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        // The 8 bins collapse to two records that meet at the recombination point,
        // which falls at the midpoint of bins 4 and 5: (4+5)/2*1 = 4.
        val lines = File("$outputDir/recombSample_imputed_path.bed").readLines().drop(1)
        assertEquals(2, lines.size, "Expected one record per side of the recombination")

        assertEquals(listOf("chr1", "0", "4", "lineA:0"), lines[0].split("\t"))
        assertEquals(listOf("chr1", "4", "8", "lineB:0"), lines[1].split("\t"))
    }

    @Test
    fun testImputeDiploidPathOutput() {
        // With inbreedCoef=1.0 and all reads supporting gamete 0, the diploid
        // path is forced homozygous and should be (lineA:0, lineA:0) throughout.
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = (1..5).map { bin ->
            PS4GData(listOf(0), Position("chr1", bin), 10)
        }
        val ps4gFile = "$tempTestDir/diploidTest.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)

        val keyFile = "$tempTestDir/diploidKey.txt"
        createKeyFile(keyFile, "diploidSample", ps4gFile)

        val outputDir = "$tempTestDir/diploidOut/"
        val result = ImputePathFromPs4g().test(
            "--path-keyfile $keyFile --out-path-dir $outputDir --path-type diploid " +
                    "--prob-correct 0.99 --prob-same 0.9999 --inbreed-coef 1.0 --bin-size 1"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val bedFile = File("$outputDir/diploidSample_imputed_path.bed")
        assertTrue(bedFile.exists(), "Output BED file was not created")

        val lines = bedFile.readLines()
        assertEquals("chrom\tstart\tend\tparent1\tparent2", lines[0])

        val dataLines = lines.drop(1)
        assertEquals(1, dataLines.size, "Expected the 5 identical bins to merge into 1 record")
        assertEquals(listOf("chr1", "0", "5", "lineA:0", "lineA:0"), dataLines[0].split("\t"))
    }

    @Test
    fun testImputeDiploidPathOutputWithHighReadCount() {
        // With inbreedCoef=1.0 and all reads supporting gamete 0, the diploid
        // path is forced homozygous and should be (lineA:0, lineA:0) throughout.
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = (1..5).map { bin ->
            PS4GData(listOf(0), Position("chr1", bin), 60)
        }
        val ps4gFile = "$tempTestDir/diploidTest.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)

        val keyFile = "$tempTestDir/diploidKey.txt"
        createKeyFile(keyFile, "diploidSample", ps4gFile)

        val outputDir = "$tempTestDir/diploidOut/"
        val result = ImputePathFromPs4g().test(
            "--path-keyfile $keyFile --out-path-dir $outputDir --path-type diploid " +
                    "--prob-correct 0.99 --prob-same 0.9999 --inbreed-coef 1.0 --bin-size 1"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val bedFile = File("$outputDir/diploidSample_imputed_path.bed")
        assertTrue(bedFile.exists(), "Output BED file was not created")

        val lines = bedFile.readLines()
        assertEquals("chrom\tstart\tend\tparent1\tparent2", lines[0])

        val dataLines = lines.drop(1)
        assertEquals(1, dataLines.size, "Expected the 5 identical bins to merge into 1 record")
        assertEquals(listOf("chr1", "0", "5", "lineA:0", "lineA:0"), dataLines[0].split("\t"))
    }

    @Test
    fun testImputeHaploidPathMultiContig() {
        // Verify that imputation runs independently per contig and produces
        // a single BED file with rows from both contigs.
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = listOf(
            PS4GData(listOf(0), Position("chr1", 1), 10),
            PS4GData(listOf(0), Position("chr1", 2), 10),
            PS4GData(listOf(1), Position("chr2", 1), 10),
            PS4GData(listOf(1), Position("chr2", 2), 10)
        )

        val ps4gFile = "$tempTestDir/multiContig.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)

        val keyFile = "$tempTestDir/multiContigKey.txt"
        createKeyFile(keyFile, "multiSample", ps4gFile)

        val outputDir = "$tempTestDir/multiContigOut/"
        val result = ImputePathFromPs4g().test(
            "--path-keyfile $keyFile --out-path-dir $outputDir --prob-correct 0.99 --prob-same 0.9999 --bin-size 1"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        // 2 bins per contig, each contig uniform, so each contig merges to a single record.
        // Merging must not run across the contig boundary even though both contigs start at bin 1.
        val lines = File("$outputDir/multiSample_imputed_path.bed").readLines().drop(1)
        assertEquals(2, lines.size, "Expected 1 merged record per contig")

        assertEquals(listOf("chr1", "0", "2", "lineA:0"), lines[0].split("\t"))
        assertEquals(listOf("chr2", "0", "2", "lineB:0"), lines[1].split("\t"))
    }

}
