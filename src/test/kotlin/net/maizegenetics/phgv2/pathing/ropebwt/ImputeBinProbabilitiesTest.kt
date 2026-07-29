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

class ImputeBinProbabilitiesTest {
    companion object {
        val tempTestDir = "${TestExtension.tempDir}imputeBinProbabilitiesTest/"

        // The column header written to every output file. NB: the command writes this string
        // without a trailing newline, so the first data row is physically concatenated onto it.
        // The parsing helpers below account for that.
        const val HEADER = "contig\tposition\tstate\tprobability\n"

        @JvmStatic
        @BeforeAll
        fun setup() {
            resetDirs()
            setupDebugLogging()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(tempTestDir).deleteRecursively()
        }

        private fun resetDirs() {
            File(tempTestDir).deleteRecursively()
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

        /**
         * Reads an output probability file into rows of [contig, position, state, probability],
         * stripping the (newline-less) header that the first data row is concatenated onto.
         */
        private fun readProbabilityRows(file: File): List<List<String>> {
            val body = file.readText().removePrefix(HEADER)
            return body.split("\n").filter { it.isNotBlank() }.map { it.split("\t") }
        }

        /** For every distinct position, returns the state whose probability is the highest. */
        private fun topStateByPosition(rows: List<List<String>>): Map<String, String> {
            return rows.groupBy { it[1] }
                .mapValues { (_, group) -> group.maxByOrNull { it[3].toDouble() }!![2] }
        }
    }

    @Test
    fun testCliktParams() {
        // Missing both --path-keyfile and --read-file
        val noInput = ImputeBinProbabilities().test("--output-dir $tempTestDir")
        assertEquals(1, noInput.statusCode)
        assertTrue(
            noInput.stderr.contains("--path-keyfile") && noInput.stderr.contains("--read-file"),
            "Expected missing-input error but got: ${noInput.stderr}"
        )

        // Missing --output-dir
        val noOutputDir = ImputeBinProbabilities().test("--read-file someFile.ps4g")
        assertEquals(1, noOutputDir.statusCode)
        assertTrue(
            noOutputDir.stderr.contains("missing option --output-dir"),
            "Expected missing-output-dir error but got: ${noOutputDir.stderr}"
        )

        // Invalid --impute-type value
        val badImputeType = ImputeBinProbabilities().test(
            "--read-file someFile.ps4g --output-dir $tempTestDir --impute-type invalid"
        )
        assertEquals(1, badImputeType.statusCode)
        assertTrue(
            badImputeType.stderr.contains("impute-type"),
            "Expected impute-type error but got: ${badImputeType.stderr}"
        )
    }

    @Test
    fun testHaploidProbabilities() {
        // Almost all reads support gamete 0 (lineA:0). A single low-count read hits gamete 1 so
        // that MostLikelyPs4gParents can still select two parents (required when --n-parents > 0).
        // The posterior should favor lineA:0 at every position.
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = (1..5).map { bin ->
            PS4GData(listOf(0), Position("chr1", bin), 50)
        } + PS4GData(listOf(1), Position("chr1", 1), 1)

        val ps4gFile = "$tempTestDir/haploid.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)
        val keyFile = "$tempTestDir/haploidKey.txt"
        createKeyFile(keyFile, "haploidSample", ps4gFile)

        val outputDir = "$tempTestDir/haploidOut/"
        val result = ImputeBinProbabilities().test(
            "--path-keyfile $keyFile --output-dir $outputDir --n-parents 2 " +
                    "--prob-correct 0.99 --prob-same 0.9999 --bin-size 1"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val outFile = File("$outputDir/haploidSample_imputed_probabilities.txt")
        assertTrue(outFile.exists(), "Output probabilities file was not created")
        assertTrue(outFile.readText().startsWith(HEADER), "Output file should start with the column header")

        val rows = readProbabilityRows(outFile)
        assertTrue(rows.isNotEmpty(), "Expected at least one probability row")

        // Every reported state name is one of the two haploid parents.
        rows.forEach { row ->
            assertEquals("chr1", row[0])
            assertTrue(row[2] == "lineA:0" || row[2] == "lineB:0", "Unexpected state ${row[2]}")
        }

        // All five bins are present (start = binSize * bin = bin, with bin-size 1).
        val topState = topStateByPosition(rows)
        assertEquals(setOf("1", "2", "3", "4", "5"), topState.keys)
        topState.forEach { (position, state) ->
            assertEquals("lineA:0", state, "Expected lineA:0 to dominate at position $position")
        }

        // Probabilities are valid values in [0, 1].
        rows.forEach { row ->
            val prob = row[3].toDouble()
            assertTrue(prob in 0.0..1.0, "Probability ${row[3]} out of range")
        }
    }

    @Test
    fun testHaploidProbabilitiesDefaultAllParents() {
        // With --n-parents omitted (default 0) the full gamete set is used and the initial-state
        // distribution is sized from the actual parent count, not the nParents option. Every read
        // supports gamete 0, so lineA:0 should dominate at every position.
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = (1..4).map { bin ->
            PS4GData(listOf(0), Position("chr1", bin), 50)
        }

        val ps4gFile = "$tempTestDir/allParents.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)
        val keyFile = "$tempTestDir/allParentsKey.txt"
        createKeyFile(keyFile, "allParentsSample", ps4gFile)

        val outputDir = "$tempTestDir/allParentsOut/"
        val result = ImputeBinProbabilities().test(
            "--path-keyfile $keyFile --output-dir $outputDir " +
                    "--prob-correct 0.99 --prob-same 0.9999 --bin-size 1"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val rows = readProbabilityRows(File("$outputDir/allParentsSample_imputed_probabilities.txt"))
        assertTrue(rows.isNotEmpty(), "Expected at least one probability row")
        val topState = topStateByPosition(rows)
        assertEquals(setOf("1", "2", "3", "4"), topState.keys)
        topState.forEach { (position, state) ->
            assertEquals("lineA:0", state, "Expected lineA:0 to dominate at position $position")
        }
    }

    @Test
    fun testDiploidProbabilities() {
        // inbreedCoef = 1.0 forces homozygous states. Reads mostly hit gamete 0, so the dominant
        // diploid state should be the homozygote (lineA:0, lineA:0).
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = (1..5).map { bin ->
            PS4GData(listOf(0), Position("chr1", bin), 50)
        } + PS4GData(listOf(1), Position("chr1", 1), 1)

        val ps4gFile = "$tempTestDir/diploid.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)
        val keyFile = "$tempTestDir/diploidKey.txt"
        createKeyFile(keyFile, "diploidSample", ps4gFile)

        val outputDir = "$tempTestDir/diploidOut/"
        val result = ImputeBinProbabilities().test(
            "--path-keyfile $keyFile --output-dir $outputDir --impute-type diploid --n-parents 2 " +
                    "--prob-correct 0.99 --prob-same 0.9999 --inbreed-coef 1.0 --bin-size 1"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val outFile = File("$outputDir/diploidSample_imputed_probabilities.txt")
        assertTrue(outFile.exists(), "Output probabilities file was not created")

        val rows = readProbabilityRows(outFile)
        assertTrue(rows.isNotEmpty(), "Expected at least one probability row")

        // Diploid state names are ordered pairs of parent names.
        rows.forEach { row ->
            assertTrue(row[2].startsWith("(") && row[2].endsWith(")"), "Unexpected diploid state ${row[2]}")
        }

        val topState = topStateByPosition(rows)
        assertEquals(setOf("1", "2", "3", "4", "5"), topState.keys)
        topState.forEach { (position, state) ->
            assertEquals("(lineA:0,lineA:0)", state, "Expected homozygous lineA:0 to dominate at position $position")
        }
    }

    @Test
    fun testMinProbRatioSuppressesAllRows() {
        // A min-prob-ratio greater than 1 means no state can exceed maxProb * ratio, so no data
        // rows are written and the file contains only the header.
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = (1..3).map { bin ->
            PS4GData(listOf(0), Position("chr1", bin), 50)
        } + PS4GData(listOf(1), Position("chr1", 1), 1)

        val ps4gFile = "$tempTestDir/minprob.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)
        val keyFile = "$tempTestDir/minprobKey.txt"
        createKeyFile(keyFile, "minprobSample", ps4gFile)

        val outputDir = "$tempTestDir/minprobOut/"
        val result = ImputeBinProbabilities().test(
            "--path-keyfile $keyFile --output-dir $outputDir --n-parents 2 " +
                    "--prob-correct 0.99 --min-prob-ratio 1.1 --bin-size 1"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val outFile = File("$outputDir/minprobSample_imputed_probabilities.txt")
        assertTrue(outFile.exists(), "Output probabilities file was not created")
        assertEquals(HEADER, outFile.readText(), "Expected only the header when min-prob-ratio > 1")
        assertTrue(readProbabilityRows(outFile).isEmpty(), "Expected no data rows")
    }

    @Test
    fun testMultiContig() {
        // chr1 reads support gamete 0, chr2 reads support gamete 1. Each contig should be imputed
        // independently, and the pooled read counts let MostLikelyPs4gParents pick both parents.
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = (1..3).map { bin -> PS4GData(listOf(0), Position("chr1", bin), 50) } +
                (1..3).map { bin -> PS4GData(listOf(1), Position("chr2", bin), 50) }

        val ps4gFile = "$tempTestDir/multi.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)
        val keyFile = "$tempTestDir/multiKey.txt"
        createKeyFile(keyFile, "multiSample", ps4gFile)

        val outputDir = "$tempTestDir/multiOut/"
        val result = ImputeBinProbabilities().test(
            "--path-keyfile $keyFile --output-dir $outputDir --n-parents 2 " +
                    "--prob-correct 0.99 --prob-same 0.9999 --bin-size 1"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val rows = readProbabilityRows(File("$outputDir/multiSample_imputed_probabilities.txt"))
        val contigs = rows.map { it[0] }.toSet()
        assertEquals(setOf("chr1", "chr2"), contigs, "Both contigs should appear in the output")

        // chr1 should be dominated by lineA:0 and chr2 by lineB:0.
        val chr1Top = topStateByPosition(rows.filter { it[0] == "chr1" })
        val chr2Top = topStateByPosition(rows.filter { it[0] == "chr2" })
        chr1Top.values.forEach { assertEquals("lineA:0", it) }
        chr2Top.values.forEach { assertEquals("lineB:0", it) }
    }

    @Test
    fun testScafContigsExcluded() {
        // Contigs whose name starts with "scaf" must be excluded from imputation and output.
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = (1..3).map { bin -> PS4GData(listOf(0), Position("chr1", bin), 50) } +
                PS4GData(listOf(1), Position("chr1", 1), 1) +
                (1..3).map { bin -> PS4GData(listOf(0), Position("scaf1", bin), 50) }

        val ps4gFile = "$tempTestDir/scaf.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)
        val keyFile = "$tempTestDir/scafKey.txt"
        createKeyFile(keyFile, "scafSample", ps4gFile)

        val outputDir = "$tempTestDir/scafOut/"
        val result = ImputeBinProbabilities().test(
            "--path-keyfile $keyFile --output-dir $outputDir --n-parents 2 " +
                    "--prob-correct 0.99 --prob-same 0.9999 --bin-size 1"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val rows = readProbabilityRows(File("$outputDir/scafSample_imputed_probabilities.txt"))
        assertTrue(rows.isNotEmpty(), "Expected data rows for chr1")
        assertTrue(rows.none { it[0].startsWith("scaf") }, "scaf contigs must be excluded from output")
        assertEquals(setOf("chr1"), rows.map { it[0] }.toSet())
    }

    @Test
    fun testBinSizeScalesPosition() {
        // The reported position column equals binSize * bin index. With bin-size 100 and bins
        // 1..3, the positions should be 100, 200, 300.
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = (1..3).map { bin -> PS4GData(listOf(0), Position("chr1", bin), 50) } +
                PS4GData(listOf(1), Position("chr1", 1), 1)

        val ps4gFile = "$tempTestDir/binsize.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)
        val keyFile = "$tempTestDir/binsizeKey.txt"
        createKeyFile(keyFile, "binsizeSample", ps4gFile)

        val outputDir = "$tempTestDir/binsizeOut/"
        val result = ImputeBinProbabilities().test(
            "--path-keyfile $keyFile --output-dir $outputDir --n-parents 2 " +
                    "--prob-correct 0.99 --prob-same 0.9999 --bin-size 100"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val rows = readProbabilityRows(File("$outputDir/binsizeSample_imputed_probabilities.txt"))
        val positions = rows.map { it[1] }.toSet()
        assertEquals(setOf("100", "200", "300"), positions, "Positions should be binSize * bin index")
    }

    @Test
    fun testReadFileOptionDerivesSampleName() {
        // When --read-file is used instead of a key file, the sample name is derived from the
        // file name, so the output file is named <fileName>_imputed_probabilities.txt.
        val gametes = mapOf(
            SampleGamete("lineA", 0) to 0,
            SampleGamete("lineB", 0) to 1
        )
        val ps4gData = (1..3).map { bin -> PS4GData(listOf(0), Position("chr1", bin), 50) } +
                PS4GData(listOf(1), Position("chr1", 1), 1)

        val ps4gFile = "$tempTestDir/readFileSample.ps4g"
        createPs4gFile(ps4gFile, gametes, ps4gData)

        val outputDir = "$tempTestDir/readFileOut/"
        val result = ImputeBinProbabilities().test(
            "--read-file $ps4gFile --output-dir $outputDir --n-parents 2 " +
                    "--prob-correct 0.99 --prob-same 0.9999 --bin-size 1"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        // ReadFiles derives the sample name from the file's base name (only fastq-style suffixes
        // are stripped, so the ".ps4g" extension remains part of the name).
        val outFile = File("$outputDir/readFileSample.ps4g_imputed_probabilities.txt")
        assertTrue(outFile.exists(), "Output file for --read-file sample was not created")
        assertTrue(readProbabilityRows(outFile).isNotEmpty(), "Expected data rows")
    }
}
