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
         * Writes a PS4G file with reads on chr1 (supporting lineA:0), chr2 (supporting lineB:0),
         * and scaf1 (supporting lineA:0), and returns its path. Every contig carries reads for
         * both gametes so that MostLikelyPs4gParents can pick two parents from any single contig.
         */
        private fun createContigFilterPs4gFile(fileName: String): String {
            val gametes = mapOf(
                SampleGamete("lineA", 0) to 0,
                SampleGamete("lineB", 0) to 1
            )
            val ps4gData = (1..3).map { bin -> PS4GData(listOf(0), Position("chr1", bin), 50) } +
                    PS4GData(listOf(1), Position("chr1", 1), 1) +
                    (1..3).map { bin -> PS4GData(listOf(1), Position("chr2", bin), 50) } +
                    PS4GData(listOf(0), Position("chr2", 1), 1) +
                    (1..3).map { bin -> PS4GData(listOf(0), Position("scaf1", bin), 50) } +
                    PS4GData(listOf(1), Position("scaf1", 1), 1)

            val ps4gFile = "$tempTestDir/$fileName"
            createPs4gFile(ps4gFile, gametes, ps4gData)
            return ps4gFile
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

        // Invalid --prob-correct value (not between 0.5 and 1.0)
        val badProbCorrectValue = ImputeBinProbabilities().test(
            "--read-file someFile.ps4g --output-dir $tempTestDir --prob-correct 0.0"
        )
        assertEquals(1, badProbCorrectValue.statusCode)
        assertTrue(
            badProbCorrectValue.stderr.contains("prob-correct"),
            "Expected impute-type error but got: ${badProbCorrectValue.stderr}"
        )

        // Invalid --prob-same value (not between 0.5 and 1.0)
        val badProbSameValue = ImputeBinProbabilities().test(
            "--read-file someFile.ps4g --output-dir $tempTestDir --prob-same 0.0"
        )
        assertEquals(1, badProbSameValue.statusCode)
        assertTrue(
            badProbSameValue.stderr.contains("prob-same"),
            "Expected impute-type error but got: ${badProbSameValue.stderr}"
        )

        // Invalid --inbreed-coef value (not between 0.0 and 1.0)
        val badInbreedCoefValue = ImputeBinProbabilities().test(
            "--read-file someFile.ps4g --output-dir $tempTestDir --inbreed-coef 1.5"
        )
        assertEquals(1, badInbreedCoefValue.statusCode)
        assertTrue(
            badInbreedCoefValue.stderr.contains("inbreed-coef"),
            "Expected impute-type error but got: ${badInbreedCoefValue.stderr}"
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
    fun testAllContigsIncludedByDefault() {
        // Without --contigs-to-use every contig in the ps4g file is imputed, including scaffold
        // contigs, which used to be dropped by a "scaf" name check.
        val ps4gFile = createContigFilterPs4gFile("scaf.ps4g")
        val keyFile = "$tempTestDir/scafKey.txt"
        createKeyFile(keyFile, "scafSample", ps4gFile)

        val outputDir = "$tempTestDir/scafOut/"
        val result = ImputeBinProbabilities().test(
            "--path-keyfile $keyFile --output-dir $outputDir --n-parents 2 " +
                    "--prob-correct 0.99 --prob-same 0.9999 --bin-size 1"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val rows = readProbabilityRows(File("$outputDir/scafSample_imputed_probabilities.txt"))
        assertTrue(rows.isNotEmpty(), "Expected data rows")
        assertEquals(setOf("chr1", "chr2", "scaf1"), rows.map { it[0] }.toSet())
    }

    @Test
    fun testContigsToUseCommaSeparatedList() {
        // --contigs-to-use restricts imputation to the listed contigs; the others produce no output.
        val ps4gFile = createContigFilterPs4gFile("contigList.ps4g")
        val keyFile = "$tempTestDir/contigListKey.txt"
        createKeyFile(keyFile, "contigListSample", ps4gFile)

        val outputDir = "$tempTestDir/contigListOut/"
        val result = ImputeBinProbabilities().test(
            "--path-keyfile $keyFile --output-dir $outputDir --n-parents 2 " +
                    "--prob-correct 0.99 --prob-same 0.9999 --bin-size 1 --contigs-to-use chr1,chr2"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val rows = readProbabilityRows(File("$outputDir/contigListSample_imputed_probabilities.txt"))
        assertEquals(setOf("chr1", "chr2"), rows.map { it[0] }.toSet(), "scaf1 should be filtered out")

        // The retained contigs are still imputed correctly: chr1 reads support lineA:0 and chr2
        // reads support lineB:0.
        topStateByPosition(rows.filter { it[0] == "chr1" }).values.forEach { assertEquals("lineA:0", it) }
        topStateByPosition(rows.filter { it[0] == "chr2" }).values.forEach { assertEquals("lineB:0", it) }
    }

    @Test
    fun testContigsToUseFromFile() {
        // The option value may also name a file holding one contig per line.
        val ps4gFile = createContigFilterPs4gFile("contigFile.ps4g")
        val keyFile = "$tempTestDir/contigFileKey.txt"
        createKeyFile(keyFile, "contigFileSample", ps4gFile)

        val contigFile = "$tempTestDir/contigsToUse.txt"
        File(contigFile).writeText("chr1\nscaf1\n")

        val outputDir = "$tempTestDir/contigFileOut/"
        val result = ImputeBinProbabilities().test(
            "--path-keyfile $keyFile --output-dir $outputDir --n-parents 2 " +
                    "--prob-correct 0.99 --prob-same 0.9999 --bin-size 1 --contigs-to-use $contigFile"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val rows = readProbabilityRows(File("$outputDir/contigFileSample_imputed_probabilities.txt"))
        assertEquals(setOf("chr1", "scaf1"), rows.map { it[0] }.toSet())
    }

    @Test
    fun testContigsToUseRestrictsParentSelection() {
        // The parent set is chosen from the filtered reader, so MostLikelyPs4gParents only counts
        // reads on the requested contigs. Only chr2 is requested, and all of its reads support
        // lineB:0, so lineB:0 must dominate everywhere in the output.
        val ps4gFile = createContigFilterPs4gFile("parentFilter.ps4g")
        val keyFile = "$tempTestDir/parentFilterKey.txt"
        createKeyFile(keyFile, "parentFilterSample", ps4gFile)

        val outputDir = "$tempTestDir/parentFilterOut/"
        val result = ImputeBinProbabilities().test(
            "--path-keyfile $keyFile --output-dir $outputDir " +
                    "--prob-correct 0.99 --prob-same 0.9999 --bin-size 1 --contigs-to-use chr2"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val rows = readProbabilityRows(File("$outputDir/parentFilterSample_imputed_probabilities.txt"))
        assertTrue(rows.isNotEmpty(), "Expected data rows for chr2")
        assertEquals(setOf("chr2"), rows.map { it[0] }.toSet())
        topStateByPosition(rows).values.forEach { assertEquals("lineB:0", it) }
    }

    @Test
    fun testContigsToUseMatchingNoContigs() {
        // A contig list that matches nothing in the ps4g file leaves the reader empty, so only the
        // header is written.
        val ps4gFile = createContigFilterPs4gFile("noMatch.ps4g")
        val keyFile = "$tempTestDir/noMatchKey.txt"
        createKeyFile(keyFile, "noMatchSample", ps4gFile)

        val outputDir = "$tempTestDir/noMatchOut/"
        val result = ImputeBinProbabilities().test(
            "--path-keyfile $keyFile --output-dir $outputDir " +
                    "--prob-correct 0.99 --prob-same 0.9999 --bin-size 1 --contigs-to-use chr99"
        )
        assertEquals(0, result.statusCode, "Command failed:\n${result.stderr}")

        val outFile = File("$outputDir/noMatchSample_imputed_probabilities.txt")
        assertEquals(HEADER, outFile.readText(), "Expected only the header when no contig matches")
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
