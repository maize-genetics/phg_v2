package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.io.TempDir
import java.io.File

/**
 * The fixtures in `data/test/imputeVcfTests` are built so the right answer is known by construction.
 *
 * Panel: three haploid founders. `founderA` carries REF at every site, `founderB` carries ALT at every
 * site, and `founderC` alternates, so any sample matching a founder matches exactly one.
 *
 * chr1 has twelve sites in two clusters -- 100 kb to 600 kb, then 1.6 Mb to 2.1 Mb -- with a 1 Mb gap
 * between them. That gap is where a distance-scaled transition makes a switch cheapest.
 *
 * Samples: `pureA` is `0/0` throughout, `pureB` is `1/1`, `hetAB` is `0/1`, `recomb` follows founderA
 * across the first cluster and founderB across the second, and `gappy` is `pureA` with two no-calls.
 * The sample VCF also carries one site the panel lacks, for the diagnostic.
 */
class ImputePathFromVcfTest {

    private val testDir = "data/test/imputeVcfTests"
    private val sampleVcf = "$testDir/sample.vcf"
    private val panelVcf = "$testDir/panel.vcf"

    private fun command(vararg extra: String, outDir: File): ImputePathFromVcf {
        val command = ImputePathFromVcf()
        val result = command.test(
            "--to-impute-vcf $sampleVcf --panel-vcf $panelVcf --out-path-dir ${outDir.absolutePath} " +
                    extra.joinToString(" ")
        )
        assertEquals(0, result.statusCode, "command failed:\n${result.stderr}")
        return command
    }

    /**
     * One sample's BED as (contig, range, parents), the range expressed in the **1-based inclusive**
     * coordinates a VCF position is compared against. BedToVcf makes the same conversion --
     * `(start + 1, end)` -- so probing with a raw VCF position against a 0-based half-open range would
     * be off by one at both edges and would wrongly report the last site as uncovered.
     */
    private fun readBed(file: File): List<Triple<String, IntRange, String>> =
        file.readLines().drop(1).filter { it.isNotBlank() }.map { line ->
            val f = line.split('\t')
            Triple(f[0], (f[1].toInt() + 1)..f[2].toInt(), f.drop(3).joinToString("/"))
        }

    private fun callAt(bed: List<Triple<String, IntRange, String>>, contig: String, position: Int) =
        bed.firstOrNull { it.first == contig && position in it.second }?.third

    // ------------------------------------------------------------------ pairing

    @Test
    fun onlySitesInBothFilesAreUsed() {
        var sites = 0
        val counts = pairVcfSites(File(sampleVcf), File(panelVcf)) { sites += it.nSites }
        assertEquals(15, counts.sitesUsed, "twelve on chr1 and three on chr2")
        assertEquals(15, sites)
        assertEquals(1, counts.sampleSitesNotInPanel, "chr1:2500000 is not in the panel")
        assertEquals(0, counts.panelSitesNotInSample)
    }

    @Test
    fun sitesArriveGroupedByContigWithTheirPositions() {
        val seen = mutableListOf<Pair<String, Int>>()
        pairVcfSites(File(sampleVcf), File(panelVcf)) { seen.add(it.contig to it.nSites) }
        assertEquals(listOf("chr1" to 12, "chr2" to 3), seen)
    }

    @Test
    fun theEncodingPutsTheFoundersAndSamplesOnAConsistentAlleleIndex() {
        lateinit var chr1: ContigSites
        pairVcfSites(File(sampleVcf), File(panelVcf)) { if (it.contig == "chr1") chr1 = it }
        val founderA = chr1.founderNames.indexOf("founderA")
        val founderB = chr1.founderNames.indexOf("founderB")
        val pureA = chr1.sampleNames.indexOf("pureA")
        val pureB = chr1.sampleNames.indexOf("pureB")
        // REF is shared by both records so it indexes to 0; founderA is REF everywhere, founderB ALT.
        assertEquals(chr1.founderAllele1(0, founderA), chr1.sampleAllele1(0, pureA),
            "pureA's allele must encode the same as founderA's")
        assertEquals(chr1.founderAllele1(0, founderB), chr1.sampleAllele1(0, pureB))
        assertNotEquals(chr1.founderAllele1(0, founderA), chr1.founderAllele1(0, founderB))
        // a haploid panel call is stored as homozygous, which is what makes it one gamete
        assertEquals(chr1.founderAllele1(0, founderA), chr1.founderAllele2(0, founderA))
    }

    @Test
    fun aNoCallIsEncodedAsMissingInBothHalves() {
        lateinit var chr1: ContigSites
        pairVcfSites(File(sampleVcf), File(panelVcf)) { if (it.contig == "chr1") chr1 = it }
        val gappy = chr1.sampleNames.indexOf("gappy")
        // sites 2 and 8 are "./." for gappy
        assertEquals(ContigSites.MISSING, chr1.sampleAllele1(2, gappy))
        assertEquals(ContigSites.MISSING, chr1.sampleAllele2(2, gappy))
        assertNotEquals(ContigSites.MISSING, chr1.sampleAllele1(1, gappy), "its other sites are called")
    }

    @Test
    fun contigsToUseRestrictsWhichContigsArePaired() {
        val seen = mutableListOf<String>()
        val counts = pairVcfSites(File(sampleVcf), File(panelVcf), setOf("chr2")) { seen.add(it.contig) }
        assertEquals(listOf("chr2"), seen)
        assertEquals(3, counts.sitesUsed)
    }

    // ------------------------------------------------------------------ the paths

    @Test
    fun aSampleIdenticalToOneFounderIsAssignedThatFounder(@TempDir dir: File) {
        command(outDir = dir)
        assertTrue(readBed(File(dir, "pureA_imputed_path.bed")).all { it.third == "founderA" },
            "pureA carries REF at every site, as only founderA does")
        assertTrue(readBed(File(dir, "pureB_imputed_path.bed")).all { it.third == "founderB" })
    }

    @Test
    fun aHeterozygousSampleIsAssignedThePairThatExplainsIt(@TempDir dir: File) {
        command("--path-type diploid", outDir = dir)
        val bed = readBed(File(dir, "hetAB_imputed_path.bed"))
        assertTrue(bed.all { it.third == "founderA/founderB" || it.third == "founderB/founderA" },
            "0/1 everywhere is explained by the A,B pair and nothing else: ${bed.map { it.third }.distinct()}")
    }

    @Test
    fun aRecombinantIsFollowedAcrossTheSwitch(@TempDir dir: File) {
        command(outDir = dir)
        val bed = readBed(File(dir, "recomb_imputed_path.bed"))
        assertEquals("founderA", callAt(bed, "chr1", 100_000), "first cluster follows founderA")
        assertEquals("founderA", callAt(bed, "chr1", 600_000))
        assertEquals("founderB", callAt(bed, "chr1", 1_600_000), "second cluster follows founderB")
        assertEquals("founderB", callAt(bed, "chr1", 2_100_000))
    }

    @Test
    fun theSwitchFallsInTheGapBetweenTheClusters(@TempDir dir: File) {
        // The midpoint cut rule places the boundary halfway between the last site of one founder and
        // the first of the next, so it should land inside the 1 Mb gap rather than beside a site.
        command(outDir = dir)
        val bed = readBed(File(dir, "recomb_imputed_path.bed")).filter { it.first == "chr1" }
        val boundary = bed.zipWithNext().firstOrNull { it.first.third != it.second.third }
            ?: fail("recomb should have a switch on chr1")
        // first of the 1-based inclusive range, so one past the BED start
        val cut = boundary.second.second.first - 1
        assertTrue(cut in 600_000..1_600_000, "the cut should fall in the gap, was $cut")
        assertEquals((600_000 + 1_600_000) / 2, cut, "and at the midpoint of the two flanking sites")
    }

    @Test
    fun aMissingGenotypeDoesNotDerailThePath(@TempDir dir: File) {
        // gappy is pureA with two no-calls. Those sites contribute nothing, so the answer should be
        // unchanged rather than switching around the gaps.
        command(outDir = dir)
        assertTrue(readBed(File(dir, "gappy_imputed_path.bed")).all { it.third == "founderA" },
            "two uninformative sites must not move the path")
    }

    @Test
    fun everySampleGetsItsOwnFile(@TempDir dir: File) {
        command(outDir = dir)
        assertEquals(
            listOf("gappy", "hetAB", "pureA", "pureB", "recomb").map { "${it}_imputed_path.bed" }.sorted(),
            dir.listFiles()!!.map { it.name }.sorted()
        )
    }

    @Test
    fun bothContigsAppearInEachPath(@TempDir dir: File) {
        command(outDir = dir)
        val contigs = readBed(File(dir, "pureA_imputed_path.bed")).map { it.first }.distinct()
        assertEquals(listOf("chr1", "chr2"), contigs)
    }

    // ------------------------------------------------------------------ output shape

    @Test
    fun aHaploidPathWritesFourColumnsAndADiploidPathFive(@TempDir dir: File) {
        command(outDir = dir)
        val haploid = File(dir, "pureA_imputed_path.bed").readLines()
        assertEquals("chrom\tstart\tend\tparent1", haploid[0])
        assertEquals(4, haploid[1].split('\t').size)

        val diploidDir = File(dir, "diploid").also { it.mkdirs() }
        command("--path-type diploid", outDir = diploidDir)
        val diploid = File(diploidDir, "pureA_imputed_path.bed").readLines()
        assertEquals("chrom\tstart\tend\tparent1\tparent2", diploid[0])
        assertEquals(5, diploid[1].split('\t').size)
    }

    @Test
    fun theBedIsReadableByBedToVcf(@TempDir dir: File) {
        // The two commands have to interoperate -- that is the point of the pipeline -- so the path
        // this writes must come back out of BedToVcf under the right sample name.
        command("--path-type diploid", outDir = dir)
        val paths = BedToVcf().readPaths(dir)
        assertEquals(listOf("gappy", "hetAB", "pureA", "pureB", "recomb"), paths.keys.sorted())
        // hetAB rather than pureA: a fully inbred sample meets the F = 0 initial-state artifact
        // pinned in theFirstSiteOfADiploidPathCannotBeHomozygousAtFZero.
        assertEquals(Pair("founderA", "founderB"), paths["hetAB"]!!["chr1"]!!.get(100_000))
    }

    @Test
    fun theImputedPathComposesIntoAVcfThroughBedToVcf(@TempDir dir: File) {
        // End to end: infer a path from the panel, then compose it back against the same panel. hetAB
        // is heterozygous A/B everywhere, so every composed genotype must carry one of each.
        val bedDir = File(dir, "bed").also { it.mkdirs() }
        command("--path-type diploid", outDir = bedDir)
        val out = File(dir, "imputed.vcf")
        val result = BedToVcf().test(
            "--bed-dir ${bedDir.absolutePath} --reference-panel-vcf $panelVcf " +
                    "--output-file ${out.absolutePath}")
        assertEquals(0, result.statusCode, result.stderr)
        val genotypes = htsjdk.variant.vcf.VCFFileReader(out, false).use { reader ->
            reader.map { it.getGenotype("hetAB").genotypeString }
        }
        assertTrue(genotypes.all { it == "A/C" || it == "G/T" },
            "hetAB should carry one REF and one ALT at every panel site: ${genotypes.distinct()}")
    }

    // ------------------------------------------------------------------ options and failures

    @Test
    fun probSwitchGovernsHowReadilyThePathMoves(@TempDir dir: File) {
        // A severe penalty should refuse the recombinant's switch; a permissive one should take it.
        val strict = File(dir, "strict").also { it.mkdirs() }
        command("--prob-switch 1e-12", outDir = strict)
        val strictCalls = readBed(File(strict, "recomb_imputed_path.bed"))
            .filter { it.first == "chr1" }.map { it.third }.distinct()
        assertEquals(1, strictCalls.size, "at 1e-12 per Mb the path should not move: $strictCalls")

        val loose = File(dir, "loose").also { it.mkdirs() }
        command("--prob-switch 1e-3", outDir = loose)
        val looseCalls = readBed(File(loose, "recomb_imputed_path.bed"))
            .filter { it.first == "chr1" }.map { it.third }.distinct()
        assertEquals(2, looseCalls.size, "at 1e-3 per Mb the switch should be taken: $looseCalls")
    }

    @Test
    fun testCliktParams(@TempDir dir: File) {
        val noSample = ImputePathFromVcf().test("--panel-vcf $panelVcf --out-path-dir ${dir.path}")
        assertEquals(1, noSample.statusCode)
        assertTrue(noSample.stderr.contains("missing option --to-impute-vcf"), noSample.stderr)

        val noPanel = ImputePathFromVcf().test("--to-impute-vcf $sampleVcf --out-path-dir ${dir.path}")
        assertEquals(1, noPanel.statusCode)
        assertTrue(noPanel.stderr.contains("missing option --panel-vcf"), noPanel.stderr)

        val noOut = ImputePathFromVcf().test("--to-impute-vcf $sampleVcf --panel-vcf $panelVcf")
        assertEquals(1, noOut.statusCode)
        assertTrue(noOut.stderr.contains("missing option --out-path-dir"), noOut.stderr)
    }

    @Test
    fun anIntermediateInbreedingCoefficientIsRefusedWithAReason(@TempDir dir: File) {
        // Distance-scaled transitions cannot go through the general scan, so 0.5 is rejected up front
        // rather than silently applying one step's transition everywhere.
        val result = ImputePathFromVcf().test(
            "--to-impute-vcf $sampleVcf --panel-vcf $panelVcf --out-path-dir ${dir.path} " +
                    "--path-type diploid --inbreed-coef 0.5")
        assertEquals(1, result.statusCode)
        assertTrue(result.stderr.contains("0.0 or 1.0"), result.stderr)
    }

    @Test
    fun outOfRangeParametersAreRefused(@TempDir dir: File) {
        for (bad in listOf("--prob-switch 0.0", "--prob-switch 1.0", "--prob-correct 0.4",
            "--prob-switch-distance 0")) {
            val result = ImputePathFromVcf().test(
                "--to-impute-vcf $sampleVcf --panel-vcf $panelVcf --out-path-dir ${dir.path} $bad")
            assertEquals(1, result.statusCode, "$bad should be refused")
        }
    }

    @Test
    fun anEmptyContigsToUseIsReportedRatherThanMeaningEverything(@TempDir dir: File) {
        val contigFile = File(dir, "empty.txt").also { it.writeText("\n\n") }
        val result = ImputePathFromVcf().test(
            "--to-impute-vcf $sampleVcf --panel-vcf $panelVcf --out-path-dir ${dir.path} " +
                    "--contigs-to-use ${contigFile.absolutePath}")
        assertEquals(1, result.statusCode)
        assertTrue(result.stderr.contains("yielded no contigs"), result.stderr)
    }

    @Test
    fun aContigTheSampleHasButThePanelHeaderDoesNotIsReported(@TempDir dir: File) {
        val odd = File(dir, "odd.vcf")
        odd.writeText(
            """
            ##fileformat=VCFv4.2
            ##contig=<ID=chrZ,length=1000>
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1
            chrZ	100	.	A	C	.	.	.	GT	0/0
            """.trimIndent() + "\n"
        )
        val error = assertThrows(IllegalArgumentException::class.java) {
            pairVcfSites(odd, File(panelVcf)) { }
        }
        assertTrue(error.message!!.contains("not in the panel's"), error.message)
    }

    @Test
    fun anUnsortedInputIsReportedRatherThanPairingNothing(@TempDir dir: File) {
        val unsorted = File(dir, "unsorted.vcf")
        unsorted.writeText(
            """
            ##fileformat=VCFv4.2
            ##contig=<ID=chr1,length=3000000>
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1
            chr1	200000	.	A	C	.	.	.	GT	0/0
            chr1	100000	.	A	C	.	.	.	GT	0/0
            """.trimIndent() + "\n"
        )
        val error = assertThrows(IllegalArgumentException::class.java) {
            pairVcfSites(unsorted, File(panelVcf)) { }
        }
        assertTrue(error.message!!.contains("not coordinate sorted"), error.message)
    }

    @Test
    fun theFirstSiteOfADiploidPathCannotBeHomozygousAtFZero(@TempDir dir: File) {
        // Not a defect in this command, but worth pinning because it surprises: at an inbreeding
        // coefficient of zero the initial distribution gives a homozygous state probability
        // F / nFounders = 0, which lnOrFloor turns into -1e6. A homozygous state is therefore a priori
        // impossible at the *first* position of a contig, whatever the evidence, so a fully inbred
        // sample is called heterozygous there and then switches to the right answer.
        //
        // pureA is 0/0 at every site and founderA is REF at every site, so the correct diploid call is
        // (founderA, founderA) throughout. What comes out is one short wrong interval and then the
        // right one. `--path-type haploid`, which is the right mode for an inbred sample, has no such
        // artifact -- aSampleIdenticalToOneFounderIsAssignedThatFounder covers that.
        //
        // The same holds for impute-path-from-ps4g, since the initial distribution is shared.
        command("--path-type diploid", outDir = dir)
        val bed = readBed(File(dir, "pureA_imputed_path.bed")).filter { it.first == "chr1" }
        assertTrue(bed.size >= 2, "expected a short wrong interval then the right one: $bed")
        assertNotEquals("founderA/founderA", bed.first().third,
            "the first interval cannot be homozygous at F = 0")
        assertEquals("founderA/founderA", bed.last().third,
            "and the path corrects itself once transitions are available")
        assertTrue(bed.first().second.last < 200_000,
            "the artifact is confined to the start: ${bed.first().second}")
    }
}
