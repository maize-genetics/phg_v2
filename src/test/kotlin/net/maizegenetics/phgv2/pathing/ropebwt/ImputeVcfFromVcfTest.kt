package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import htsjdk.variant.vcf.VCFFileReader
import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.io.TempDir
import java.io.File

/**
 * Uses the same fixtures as [ImputePathFromVcfTest], plus `panelHighDensity.vcf`: the same three
 * founders with `founderD` added, 98 sites against the path panel's 15. `founderA` still carries REF
 * everywhere and `founderB` ALT, so a sample assigned to a founder has a predictable output genotype at
 * every one of the dense panel's sites -- which is the whole point being tested, that the output carries
 * the dense panel's sites rather than the sample's.
 */
class ImputeVcfFromVcfTest {

    private val testDir = "data/test/imputeVcfTests"
    private val sampleVcf = "$testDir/sample.vcf"
    private val panelVcf = "$testDir/panel.vcf"
    private val densePanelVcf = "$testDir/panelHighDensity.vcf"

    private fun run(outputFile: File, vararg extra: String) =
        ImputeVcfFromVcf().test(
            "--to-impute-vcf $sampleVcf --panel-vcf $panelVcf " +
                    "--high-density-panel-vcf $densePanelVcf --output-file ${outputFile.absolutePath} " +
                    extra.joinToString(" ")
        )

    private fun genotypesOf(vcf: File, sample: String): List<String> =
        VCFFileReader(vcf, false).use { reader -> reader.map { it.getGenotype(sample).genotypeString } }

    private fun sitesOf(vcf: File): List<String> =
        VCFFileReader(vcf, false).use { reader -> reader.map { "${it.contig}:${it.start}" } }

    @Test
    fun theOutputCarriesTheDensePanelsSitesNotTheSamplesOwn(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        assertEquals(0, run(out).statusCode)
        assertEquals(sitesOf(File(densePanelVcf)), sitesOf(out),
            "raising density is the point: the dense panel's sites are the output's sites")
        assertEquals(98, sitesOf(out).size, "against the sample VCF's own 16")
    }

    @Test
    fun aSampleMatchingOneFounderGetsThatFoundersAllelesThroughout(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        assertEquals(0, run(out, "--extend-to-contig-ends").statusCode)
        // pureA is 0/0 at its 15 shared sites, which only founderA explains; founderA is REF at every
        // dense site, so every one of the 98 output genotypes must be homozygous REF. The flag is what
        // makes that "every": without it the dense sites outside the shared span are not called.
        val pureA = genotypesOf(out, "pureA")
        assertEquals(98, pureA.size)
        assertTrue(pureA.all { it == "A/A" || it == "G/G" },
            "pureA should be homozygous REF at every dense site: ${pureA.distinct()}")
        // pureB is the mirror image
        val pureB = genotypesOf(out, "pureB")
        assertTrue(pureB.all { it == "C/C" || it == "T/T" },
            "pureB should be homozygous ALT: ${pureB.distinct()}")
    }

    @Test
    fun aHeterozygousSampleComesOutHeterozygousAtEveryDenseSite(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        assertEquals(0, run(out, "--path-type diploid", "--extend-to-contig-ends").statusCode)
        val hetAB = genotypesOf(out, "hetAB")
        assertTrue(hetAB.all { it == "A/C" || it == "G/T" },
            "one REF and one ALT at every site: ${hetAB.distinct()}")
    }

    @Test
    fun aRecombinantSwitchesAllelesWhereItsPathSwitches(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        assertEquals(0, run(out).statusCode)
        // recomb follows founderA over chr1's first cluster and founderB over the second. The path's
        // boundary is the midpoint of 600 kb and 1.6 Mb, so dense sites below 1.1 Mb take founderA's
        // REF and those above take founderB's ALT.
        val byPosition = VCFFileReader(out, false).use { reader ->
            reader.filter { it.contig == "chr1" }
                .associate { it.start to it.getGenotype("recomb").genotypeString }
        }
        assertEquals("A/A", byPosition[500_000], "well before the boundary")
        assertEquals("A/A", byPosition[1_100_000 - 50_000], "the last dense site before it")
        assertEquals("C/C", byPosition[1_150_000], "the first dense site after it")
        assertEquals("C/C", byPosition[2_000_000], "well after")
    }

    @Test
    fun theSamplesAreTheOnesImputedNotThePanelFounders(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        assertEquals(0, run(out).statusCode)
        VCFFileReader(out, false).use { reader ->
            assertEquals(listOf("gappy", "hetAB", "pureA", "pureB", "recomb"),
                reader.fileHeader.sampleNamesInOrder.sorted())
        }
    }

    // ------------------------------------------------------------------ intermediates

    @Test
    fun nothingIntermediateIsWrittenUnlessAsked(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        assertEquals(0, run(out).statusCode)
        assertEquals(listOf("imputed.vcf"), dir.listFiles()!!.map { it.name },
            "the output and nothing else")
    }

    @Test
    fun bedDirWritesThePathsAsWell(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        val beds = File(dir, "beds")
        assertEquals(0, run(out, "--bed-dir ${beds.absolutePath}").statusCode)
        assertEquals(
            listOf("gappy", "hetAB", "pureA", "pureB", "recomb").map { "${it}_imputed_path.bed" },
            beds.listFiles()!!.map { it.name }.sorted()
        )
        assertTrue(out.exists(), "and the VCF is still written")
    }

    @Test
    fun theChainedCommandMatchesRunningTheTwoStepsSeparately(@TempDir dir: File) {
        // The chained form must be a shortcut, not a second implementation. Run the two commands by
        // hand through a BED on disk and compare against the in-memory chain.
        val bedDir = File(dir, "bed").also { it.mkdirs() }
        assertEquals(0, ImputePathFromVcf().test(
            "--to-impute-vcf $sampleVcf --panel-vcf $panelVcf --path-type diploid " +
                    "--out-path-dir ${bedDir.absolutePath}").statusCode)
        val viaFiles = File(dir, "viaFiles.vcf")
        assertEquals(0, BedToVcf().test(
            "--bed-dir ${bedDir.absolutePath} --reference-panel-vcf $densePanelVcf " +
                    "--output-file ${viaFiles.absolutePath}").statusCode)

        val chained = File(dir, "chained.vcf")
        assertEquals(0, run(chained, "--path-type diploid").statusCode)

        assertEquals(viaFiles.readLines(), chained.readLines(),
            "the in-memory chain and the two-file route must agree line for line")
    }

    @Test
    fun theInMemoryPathAdapterAgreesWithReadingTheSameBed(@TempDir dir: File) {
        // pathsFromIntervals and readPath have to apply the same coordinate conversion, or a path would
        // land differently depending on whether it went through a file.
        val bedDir = File(dir, "bed").also { it.mkdirs() }
        assertEquals(0, ImputePathFromVcf().test(
            "--to-impute-vcf $sampleVcf --panel-vcf $panelVcf --path-type diploid " +
                    "--out-path-dir ${bedDir.absolutePath}").statusCode)
        val fromFiles = BedToVcf().readPaths(bedDir)

        val result = imputeFounderPaths(File(sampleVcf), File(panelVcf),
            VcfPathParameters(pathType = "diploid"))
        val fromMemory = BedToVcf().pathsFromIntervals(result.paths)

        assertEquals(fromFiles.keys.sorted(), fromMemory.keys.sorted())
        for (sample in fromFiles.keys) {
            for (contig in fromFiles[sample]!!.keys) {
                val viaFile = fromFiles[sample]!![contig]!!
                val viaMemory = fromMemory[sample]!![contig]!!
                for (position in listOf(1, 100_000, 600_000, 1_100_000, 2_100_000)) {
                    assertEquals(viaFile.get(position), viaMemory.get(position),
                        "$sample $contig:$position")
                }
            }
        }
    }

    // ------------------------------------------------------------------ validation

    @Test
    fun aFounderMissingFromTheDensePanelIsReportedBeforeAnyWork(@TempDir dir: File) {
        // Composition looks founders up by name and emits a no-call where one is absent, so a naming
        // mismatch would otherwise produce a whole VCF of no-calls with nothing to say why.
        val thin = File(dir, "thin.vcf")
        thin.writeText(
            """
            ##fileformat=VCFv4.2
            ##contig=<ID=chr1,length=3000000>
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	founderA
            chr1	100000	.	A	C	.	.	.	GT	0
            """.trimIndent() + "\n"
        )
        val out = File(dir, "imputed.vcf")
        val result = ImputeVcfFromVcf().test(
            "--to-impute-vcf $sampleVcf --panel-vcf $panelVcf " +
                    "--high-density-panel-vcf ${thin.absolutePath} --output-file ${out.absolutePath}")
        assertEquals(1, result.statusCode)
        assertTrue(result.stderr.contains("founderB"), result.stderr)
        assertTrue(result.stderr.contains("must be a subset"), result.stderr)
        assertFalse(out.exists(), "and nothing is written")
    }

    @Test
    fun aDensePanelWithExtraFoundersIsFine(@TempDir dir: File) {
        // panelHighDensity.vcf carries founderD, which the path panel does not. A superset is allowed:
        // the extra founder is simply never asked for.
        val out = File(dir, "imputed.vcf")
        assertEquals(0, run(out).statusCode)
        assertTrue(out.exists())
    }

    @Test
    fun usingOnePanelForBothStepsWorks(@TempDir dir: File) {
        // Allowed, and nothing is gained -- the output has the path panel's own sites.
        val out = File(dir, "imputed.vcf")
        val result = ImputeVcfFromVcf().test(
            "--to-impute-vcf $sampleVcf --panel-vcf $panelVcf " +
                    "--high-density-panel-vcf $panelVcf --output-file ${out.absolutePath}")
        assertEquals(0, result.statusCode, result.stderr)
        assertEquals(15, sitesOf(out).size, "the path panel's own 15 sites")
    }

    @Test
    fun noSharedSitesIsReportedRatherThanWritingNoCalls(@TempDir dir: File) {
        // A sample and panel on different coordinates share nothing, which would otherwise compose into
        // a complete VCF of no-calls.
        val elsewhere = File(dir, "elsewhere.vcf")
        elsewhere.writeText(
            """
            ##fileformat=VCFv4.2
            ##contig=<ID=chr1,length=3000000>
            ##contig=<ID=chr2,length=1000000>
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	founderA	founderB	founderC
            chr1	999999	.	A	C	.	.	.	GT	0	1	0
            """.trimIndent() + "\n"
        )
        val out = File(dir, "imputed.vcf")
        val result = ImputeVcfFromVcf().test(
            "--to-impute-vcf $sampleVcf --panel-vcf ${elsewhere.absolutePath} " +
                    "--high-density-panel-vcf $densePanelVcf --output-file ${out.absolutePath}")
        assertEquals(1, result.statusCode)
        assertTrue(result.stderr.contains("shared no sites"), result.stderr)
    }

    @Test
    fun testCliktParams(@TempDir dir: File) {
        val out = "${dir.path}/out.vcf"
        for ((missing, flag) in listOf(
            "--to-impute-vcf" to "--panel-vcf $panelVcf --high-density-panel-vcf $densePanelVcf --output-file $out",
            "--panel-vcf" to "--to-impute-vcf $sampleVcf --high-density-panel-vcf $densePanelVcf --output-file $out",
            "--high-density-panel-vcf" to "--to-impute-vcf $sampleVcf --panel-vcf $panelVcf --output-file $out",
            "--output-file" to "--to-impute-vcf $sampleVcf --panel-vcf $panelVcf --high-density-panel-vcf $densePanelVcf"
        )) {
            val result = ImputeVcfFromVcf().test(flag)
            assertEquals(1, result.statusCode, "$missing should be required")
            assertTrue(result.stderr.contains("missing option $missing"), result.stderr)
        }
    }

    @Test
    fun anIntermediateInbreedingCoefficientIsRefused(@TempDir dir: File) {
        val result = run(File(dir, "out.vcf"), "--path-type diploid --inbreed-coef 0.5")
        assertEquals(1, result.statusCode)
        assertTrue(result.stderr.contains("0.0 or 1.0"), result.stderr)
    }

    @Test
    fun aPathClaimsNothingBeyondItsTerminalSites(@TempDir dir: File) {
        // The default. pathToIntervals is asymmetric -- its first interval starts at 0 while its last
        // ends at the last observation -- so both ends are corrected here: the path spans the first to
        // the last site shared with the panel, and dense sites outside that span get no call. There is
        // no evidence of ancestry past the terminal markers, and a telomere-proximal recombination is
        // exactly where it would be missed.
        val out = File(dir, "imputed.vcf")
        val beds = File(dir, "beds")
        assertEquals(0, run(out, "--bed-dir ${beds.absolutePath}").statusCode)

        val bed = File(beds, "pureA_imputed_path.bed").readLines().drop(1).map { it.split('\t') }
        fun firstStart(contig: String) = bed.first { it[0] == contig }[1].toInt()
        fun lastEnd(contig: String) = bed.last { it[0] == contig }[2].toInt()
        // shared sites run chr1 100 kb to 2.1 Mb and chr2 50 kb to 250 kb; BED is 0-based half-open,
        // so a start one below the first site makes that site the first base covered.
        assertEquals(99_999, firstStart("chr1"))
        assertEquals(2_100_000, lastEnd("chr1"))
        assertEquals(49_999, firstStart("chr2"))
        assertEquals(250_000, lastEnd("chr2"))

        val called = VCFFileReader(out, false).use { reader ->
            reader.associate { "${it.contig}:${it.start}" to (it.getGenotype("pureA").genotypeString != "./.") }
        }
        assertEquals(false, called["chr1:50000"], "before the first shared site")
        assertEquals(true, called["chr1:100000"], "the first shared site itself is covered")
        assertEquals(true, called["chr1:2100000"], "and so is the last")
        assertEquals(false, called["chr1:2150000"], "past the last shared site")
        assertTrue(called.values.count { !it } > 0, "some dense sites fall outside the span")
    }

    @Test
    fun extendToContigEndsFillsBothEnds(@TempDir dir: File) {
        // The opposite choice, for a caller who would rather assume the terminal ancestry continues.
        // Both ends, never one: filling the leading edge while leaving the trailing edge uncalled was
        // the asymmetry that left 46 of these 98 dense sites no-call.
        val out = File(dir, "imputed.vcf")
        val beds = File(dir, "beds")
        assertEquals(0, run(out, "--extend-to-contig-ends", "--bed-dir ${beds.absolutePath}").statusCode)

        val bed = File(beds, "pureA_imputed_path.bed").readLines().drop(1).map { it.split('\t') }
        assertEquals("0", bed.first { it[0] == "chr1" }[1], "chr1's path starts at the contig start")
        assertEquals("3000000", bed.last { it[0] == "chr1" }[2], "and reaches the declared length")
        assertEquals("0", bed.first { it[0] == "chr2" }[1])
        assertEquals("1000000", bed.last { it[0] == "chr2" }[2])

        assertEquals(0, genotypesOf(out, "pureA").count { it == "./." }, "every dense site is called")
    }

    @Test
    fun withoutDeclaredContigLengthsThePathStillCoversEveryDenseSite(@TempDir dir: File) {
        // ##contig lines are optional and real panels omit them -- the 25-founder maize panel has none,
        // so --extend-to-contig-ends cannot rely on a declared length. Without one the path is carried
        // to Int.MAX_VALUE, which covers anything a denser panel could hold.
        fun stripContigLines(source: String, target: File) {
            target.writeText(File(source).readLines()
                .filterNot { it.startsWith("##contig") }.joinToString("\n") + "\n")
        }
        val panelNoContigs = File(dir, "panelNoContigs.vcf").also { stripContigLines(panelVcf, it) }
        val out = File(dir, "imputed.vcf")
        val result = ImputeVcfFromVcf().test(
            "--to-impute-vcf $sampleVcf --panel-vcf ${panelNoContigs.absolutePath} " +
                    "--high-density-panel-vcf $densePanelVcf --output-file ${out.absolutePath} " +
                    "--extend-to-contig-ends")
        assertEquals(0, result.statusCode, result.stderr)
        assertEquals(0, genotypesOf(out, "pureA").count { it == "./." },
            "every dense site is still called")
    }
}
