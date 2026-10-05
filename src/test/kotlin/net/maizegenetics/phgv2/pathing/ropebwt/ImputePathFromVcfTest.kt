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
    fun aHalfCallIsKeptHalfMissingRatherThanReadAsHomozygous(@TempDir dir: File) {
        // A half call has one allele unknown, which is not the same as a haploid call's one allele.
        // Reading 0/. as 0/0 would invent an observed homozygote.
        val header = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n" +
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        val panel = File(dir, "panel.vcf").apply {
            writeText("$header\tfRef\tfAlt\tfHalf\tfFlip\nchr1\t100\t.\tA\tC\t.\t.\t.\tGT\t0\t1\t1/.\t./1\n")
        }
        val sample = File(dir, "sample.vcf").apply {
            writeText("$header\thalf\tflip\thet\thap\nchr1\t100\t.\tA\tC\t.\t.\t.\tGT\t0/.\t./1\t0/1\t1\n")
        }
        lateinit var chr1: ContigSites
        pairVcfSites(sample, panel) { chr1 = it }
        val missing = ContigSites.MISSING
        val ref: Byte = 0
        val alt = chr1.founderAllele1(0, chr1.founderNames.indexOf("fAlt"))

        // haploid calls fill both halves; diploid calls are stored as written
        assertEquals(ref, chr1.founderAllele2(0, chr1.founderNames.indexOf("fRef")), "haploid 0 is 0/0")
        assertEquals(alt, chr1.sampleAllele2(0, chr1.sampleNames.indexOf("hap")), "haploid 1 is 1/1")
        val het = chr1.sampleNames.indexOf("het")
        assertEquals(ref, chr1.sampleAllele1(0, het))
        assertEquals(alt, chr1.sampleAllele2(0, het))

        // half calls: the called allele first, MISSING second, whichever side was missing
        for ((name, allele) in listOf("half" to ref, "flip" to alt)) {
            val s = chr1.sampleNames.indexOf(name)
            assertEquals(allele, chr1.sampleAllele1(0, s), "$name keeps its called allele")
            assertEquals(missing, chr1.sampleAllele2(0, s), "$name is not read as homozygous")
        }
        for (name in listOf("fHalf", "fFlip")) {
            val f = chr1.founderNames.indexOf(name)
            assertEquals(alt, chr1.founderAllele1(0, f), "$name keeps its called allele first")
            assertEquals(missing, chr1.founderAllele2(0, f))
        }

        // A half-called sample says nothing about any state.
        val halfSample = VcfGenotypeEmissionProbability(chr1, chr1.sampleNames.indexOf("half"), 0.98)
        assertTrue(halfSample.getDiploidEmissionProbabilityArray(0).all { it == 0.0 },
            "a half call is neutral")

        // A half-called founder still contributes its one gamete: paired with itself it emits exactly
        // what the homozygous fAlt does, for a sample that is called.
        val n = chr1.nFounders
        val hap = VcfGenotypeEmissionProbability(chr1, chr1.sampleNames.indexOf("hap"), 0.98)
            .getDiploidEmissionProbabilityArray(0)
        val fAlt = chr1.founderNames.indexOf("fAlt")
        for (name in listOf("fHalf", "fFlip")) {
            val f = chr1.founderNames.indexOf(name)
            assertEquals(hap[fAlt * n + fAlt], hap[f * n + f], "$name/$name emits as fAlt/fAlt")
        }
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
        command("--path-type haploid", outDir = dir)
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
        command("--path-type haploid", outDir = dir)
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
        command("--path-type haploid", outDir = dir)
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
        command("--path-type haploid", outDir = dir)
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
    fun thePathTypeSetsTheInbreedingCoefficientAndDefaultsToDiploid(@TempDir dir: File) {
        // --path-type is the only control: haploid is the diploid path at an inbreeding coefficient of 1,
        // diploid at 0. Diploid is the default because it imputes inbred lines as well as haploid does,
        // while haploid cannot represent a heterozygote.
        assertEquals(1.0, VcfPathParameters(pathType = "haploid").coefficient)
        assertEquals(0.0, VcfPathParameters(pathType = "diploid").coefficient)
        assertEquals("diploid", VcfPathParameters().pathType)

        command(outDir = dir)
        assertEquals("chrom\tstart\tend\tparent1\tparent2",
            File(dir, "pureA_imputed_path.bed").readLines()[0], "no --path-type gives a diploid path")
    }

    @Test
    fun inbreedCoefIsNoLongerAnOption(@TempDir dir: File) {
        // It duplicated --path-type: only 0 and 1 were ever accepted, and those are diploid and haploid.
        val result = ImputePathFromVcf().test(
            "--to-impute-vcf $sampleVcf --panel-vcf $panelVcf --out-path-dir ${dir.path} --inbreed-coef 1.0")
        assertEquals(1, result.statusCode)
        assertTrue(result.stderr.contains("no such option"), result.stderr)
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
    fun aFullyInbredSampleIsHomozygousFromTheFirstSite(@TempDir dir: File) {
        // This used to fail. The initial distribution weighted a homozygous state by F / nParents,
        // which is zero at F = 0 and floors to -1e6, so a homozygous state was impossible at the first
        // position of a contig whatever the evidence: pureA came out heterozygous there and only
        // corrected itself once a transition was available. The distribution is now uniform, since the
        // coefficient's job is the transitions.
        command("--path-type diploid", outDir = dir)
        val bed = readBed(File(dir, "pureA_imputed_path.bed"))
        assertTrue(bed.all { it.third == "founderA/founderA" },
            "pureA is 0/0 at every site and founderA is REF at every site: $bed")
        assertEquals(2, bed.size, "one interval per contig, with no spurious segment at either start")
    }

    // ------------------------------------------------------------------ contig bounds

    @Test
    fun aPathSpansTheFirstToTheLastSharedSiteByDefault(@TempDir dir: File) {
        command("--path-type haploid", outDir = dir)
        val bed = readBed(File(dir, "pureA_imputed_path.bed"))
        // shared sites: chr1 100 kb to 2.1 Mb, chr2 50 kb to 250 kb. Read as 1-based inclusive, the
        // path's first base is the first shared site and its last is the last shared site.
        val chr1 = bed.filter { it.first == "chr1" }
        assertEquals(100_000, chr1.first().second.first, "starts at the first shared site, not at 1")
        assertEquals(2_100_000, chr1.last().second.last, "and stops at the last")
        assertNull(callAt(bed, "chr1", 99_999), "nothing is claimed before the first shared site")
        assertNull(callAt(bed, "chr1", 2_100_001), "nor after the last")
        assertEquals("founderA", callAt(bed, "chr1", 100_000))
        assertEquals("founderA", callAt(bed, "chr1", 2_100_000))
    }

    @Test
    fun extendToContigEndsClaimsWholeContigs(@TempDir dir: File) {
        command("--path-type haploid", "--extend-to-contig-ends", outDir = dir)
        val bed = readBed(File(dir, "pureA_imputed_path.bed"))
        val chr1 = bed.filter { it.first == "chr1" }
        assertEquals(1, chr1.first().second.first, "the leading edge reaches the contig start")
        assertEquals(3_000_000, chr1.last().second.last, "and the trailing edge the declared length")
        assertEquals("founderA", callAt(bed, "chr1", 1))
        assertEquals("founderA", callAt(bed, "chr1", 3_000_000))
    }

    @Test
    fun applyContigBoundsTrimsOrExtendsTheTerminalIntervals() {
        // Unit level, because the two ends are easy to get independently wrong and the command-level
        // tests above cannot show the degenerate case.
        val ab = Pair("founderA", "founderB")
        val cd = Pair("founderC", "founderD")
        val intervals = listOf(
            PathInterval("chr1", 0, 500, ab),
            PathInterval("chr1", 500, 900, cd)
        )
        // first site at 101, so the default start is 100: 0-based half-open, first base covered 101
        val trimmed = applyContigBounds(intervals, 101, 4000, false)
        assertEquals(listOf(100, 500), trimmed.map { it.start })
        assertEquals(listOf(500, 900), trimmed.map { it.end }, "the trailing edge is left alone")

        val extended = applyContigBounds(intervals, 101, 4000, true)
        assertEquals(listOf(0, 500), extended.map { it.start })
        assertEquals(listOf(500, 4000), extended.map { it.end })

        // no declared length: carried as far as any panel could reach
        assertEquals(Int.MAX_VALUE, applyContigBounds(intervals, 101, null, true).last().end)

        // An interval that trimming empties is dropped rather than emitted zero-length, which has no
        // 1-based inclusive form and BedToVcf would reject. The midpoint rule does not produce one --
        // the first cut never falls below the first site -- so this guards intervals from elsewhere.
        val degenerate = listOf(PathInterval("chr1", 0, 100, ab), PathInterval("chr1", 100, 400, cd))
        assertEquals(listOf(cd), applyContigBounds(degenerate, 101, 4000, false).map { it.call })

        assertEquals(emptyList<PathInterval<Pair<String, String>>>(),
            applyContigBounds(emptyList(), 101, 4000, true))
    }
}
