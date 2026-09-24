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
        assertEquals(0, run(out).statusCode)
        // pureA is 0/0 at its 15 shared sites, which only founderA explains; founderA is REF at every
        // dense site, so every one of the 98 output genotypes must be homozygous REF.
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
        assertEquals(0, run(out, "--path-type diploid").statusCode)
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
    fun thePathIsCarriedToTheEndOfEachContig(@TempDir dir: File) {
        // pathToIntervals extrapolates the leading edge back to position 0 but leaves the trailing edge
        // at the last observation. On a uniform bin grid that is invisible; on sparse VCF sites it left
        // 46 of these 98 dense sites uncalled, because the shared sites stop at chr1:2.1 Mb while the
        // dense panel runs to 2.95 Mb. Extrapolating the terminal ancestry is the same assumption the
        // leading edge already makes.
        val out = File(dir, "imputed.vcf")
        val beds = File(dir, "beds")
        assertEquals(0, run(out, "--bed-dir ${beds.absolutePath}").statusCode)

        val bed = File(beds, "pureA_imputed_path.bed").readLines().drop(1).map { it.split('\t') }
        fun endOf(contig: String) = bed.last { it[0] == contig }[2]
        assertEquals("3000000", endOf("chr1"), "chr1's path reaches the length the panel declares")
        assertEquals("1000000", endOf("chr2"), "and so does chr2's")

        val noCalls = genotypesOf(out, "pureA").count { it == "./." }
        assertEquals(0, noCalls, "every dense site is called")
    }

    @Test
    fun withoutDeclaredContigLengthsThePathStillCoversEveryDenseSite(@TempDir dir: File) {
        // ##contig lines are optional and real panels omit them -- the 25-founder maize panel has none.
        // Without a length the path is carried to Int.MAX_VALUE, which covers anything a denser panel
        // could hold.
        fun stripContigLines(source: String, target: File) {
            target.writeText(File(source).readLines()
                .filterNot { it.startsWith("##contig") }.joinToString("\n") + "\n")
        }
        val panelNoContigs = File(dir, "panelNoContigs.vcf").also { stripContigLines(panelVcf, it) }
        val out = File(dir, "imputed.vcf")
        val result = ImputeVcfFromVcf().test(
            "--to-impute-vcf $sampleVcf --panel-vcf ${panelNoContigs.absolutePath} " +
                    "--high-density-panel-vcf $densePanelVcf --output-file ${out.absolutePath}")
        assertEquals(0, result.statusCode, result.stderr)
        assertEquals(0, genotypesOf(out, "pureA").count { it == "./." },
            "every dense site is still called")
    }
}
