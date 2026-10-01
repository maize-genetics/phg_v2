package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import htsjdk.variant.vcf.VCFFileReader
import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.io.TempDir
import java.io.File

/**
 * The test data in `data/test/bedToVcfTests` is built so each case below is checkable by hand.
 *
 * Panel, four founders over two contigs:
 *
 * | site | REF | ALT | founderA | founderB | founderC | founderD |
 * |---|---|---|---|---|---|---|
 * | chr1:10 | A | C | A | C | A | C |
 * | chr1:60 | G | T | T | G | T | G |
 * | chr1:150 | C | A,T | C | A | T | A |
 * | chr1:400 | T | G | T | T | G | G |
 * | chr1:600 | C | T | C | T | T | het `0/1` |
 * | chr1:700 | A | C | A | C | het `0/1` | het `1/0` |
 * | chr1:800 | G | A | G | A | half call `1/.` | G |
 * | chr1:900 | A | G | G | G | A | no-call |
 * | chr2:50 | G | C | G | C | C | G |
 * | chr2:500 | T | A | A | T | T | A |
 *
 * Paths:
 *
 * - `sample1_imputed_path.bed` -- PHG's own naming, one file, diploid.
 *   chr1 [0,200) = (A,B), chr1 [200,1000) = (C,D), chr2 all = (A,A)
 * - `sample2_chr1_imputed.bed` + `sample2_chr2_imputed.bed` -- per-contig naming, two files that
 *   must merge into one sample. chr1 = (B,B), chr2 = (C,D)
 * - `sample3_imputed_path.bed` -- haploid, four columns, chr1 only. chr1 = founderC
 */
class BedToVcfTest {

    private val testDir = "data/test/bedToVcfTests"
    private val panel = "$testDir/panel.vcf"

    /** sample -> position -> the genotype string the VCF should carry. */
    private fun genotypesByPosition(vcf: File): Map<String, Map<String, String>> {
        val result = mutableMapOf<String, MutableMap<String, String>>()
        VCFFileReader(vcf, false).use { reader ->
            for (variant in reader) {
                for (genotype in variant.genotypes) {
                    result.getOrPut(genotype.sampleName) { mutableMapOf() }
                        .put("${variant.contig}:${variant.start}", genotype.genotypeString)
                }
            }
        }
        return result
    }

    private fun runCommand(outputFile: File) {
        val result = BedToVcf().test(
            "--bed-dir $testDir --reference-panel-vcf $panel --output-file ${outputFile.absolutePath}"
        )
        assertEquals(0, result.statusCode, "command failed:\n${result.stderr}")
    }

    @Test
    fun theThreeBedNamingConventionsAllResolveToTheirSampleName() {
        // impute-path-from-ps4g writes <sample>_imputed_path.bed. The grits naming this was ported
        // from only recognised the per-contig form, so a path written by this project's own imputer
        // would have come back as a sample called "sample1_imputed_path".
        val command = BedToVcf()
        assertEquals("sample1", command.sampleNameOf(File("sample1_imputed_path.bed")))
        assertEquals("sample2", command.sampleNameOf(File("sample2_chr1_imputed.bed")))
        assertEquals("sample2", command.sampleNameOf(File("sample2_chr10_imputed.bed")))
        assertEquals("plain", command.sampleNameOf(File("plain.bed")))
        // a sample name containing an underscore must survive both forms
        assertEquals("B73xCML103", command.sampleNameOf(File("B73xCML103_imputed_path.bed")))
        assertEquals("B73xCML103", command.sampleNameOf(File("B73xCML103_chr1_imputed.bed")))
    }

    @Test
    fun severalFilesForOneSampleMergeIntoOnePath() {
        val paths = BedToVcf().readPaths(File(testDir))
        assertEquals(listOf("sample1", "sample2", "sample3"), paths.keys.sorted())
        // sample2's two per-contig files land in one path covering both contigs
        assertEquals(setOf("chr1", "chr2"), paths["sample2"]!!.keys)
        assertEquals(setOf("chr1", "chr2"), paths["sample1"]!!.keys)
        assertEquals(setOf("chr1"), paths["sample3"]!!.keys, "sample3's path names only chr1")
    }

    @Test
    fun bedIsReadAsZeroBasedHalfOpenAndConvertedToOneBasedInclusive() {
        val paths = BedToVcf().readPaths(File(testDir))
        val chr1 = paths["sample1"]!!["chr1"]!!
        // BED [0,200) becomes 1..200 inclusive, so 200 is the last base of the first interval
        assertEquals(Pair("founderA", "founderB"), chr1.get(1))
        assertEquals(Pair("founderA", "founderB"), chr1.get(200))
        assertEquals(Pair("founderC", "founderD"), chr1.get(201))
        assertNull(chr1.get(0), "position 0 does not exist in 1-based coordinates")
        assertNull(chr1.get(1001), "past the end of the path")
    }

    @Test
    fun aFourColumnHaploidBedIsReadAsHomozygous() {
        val paths = BedToVcf().readPaths(File(testDir))
        assertEquals(Pair("founderC", "founderC"), paths["sample3"]!!["chr1"]!!.get(10),
            "a haploid path names one founder; both haplotypes take it")
    }

    @Test
    fun zeroLengthIntervalsAreSkipped(@TempDir dir: File) {
        // pathToIntervals no longer emits these, but a hand-made or older BED can contain one and
        // it has no 1-based inclusive form.
        val bed = File(dir, "degenerate_imputed_path.bed")
        bed.writeText("chrom\tstart\tend\tparent1\tparent2\nchr1\t0\t0\tfounderA\tfounderB\n" +
                "chr1\t0\t100\tfounderA\tfounderB\n")
        val paths = BedToVcf().readPaths(dir)
        assertEquals(Pair("founderA", "founderB"), paths["degenerate"]!!["chr1"]!!.get(1))
    }

    @Test
    fun genotypesComeFromTheFoundersTheyAreAssignedTo(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        runCommand(out)
        val byPosition = genotypesByPosition(out)

        // chr1:10 REF A ALT C; founders A=A B=C C=A D=C
        assertEquals("A/C", byPosition["sample1"]!!["chr1:10"], "sample1 is (founderA, founderB)")
        assertEquals("C/C", byPosition["sample2"]!!["chr1:10"], "sample2 is (founderB, founderB)")
        assertEquals("A/A", byPosition["sample3"]!!["chr1:10"], "sample3 is haploid founderC")

        // chr1:60 REF G ALT T; A=T B=G C=T D=G
        assertEquals("T/G", byPosition["sample1"]!!["chr1:60"])
        // chr1:150 is multi-allelic, REF C ALT A,T; A=C B=A C=T D=A
        assertEquals("C/A", byPosition["sample1"]!!["chr1:150"], "multi-allelic site")
        assertEquals("T/T", byPosition["sample3"]!!["chr1:150"])
    }

    @Test
    fun thePathSwitchMovesWhichFoundersAreRead(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        runCommand(out)
        val sample1 = genotypesByPosition(out)["sample1"]!!
        // sample1 switches from (A,B) to (C,D) at BED 200, i.e. after 1-based position 200
        assertEquals("C/A", sample1["chr1:150"], "before the switch, founders A and B")
        assertEquals("G/G", sample1["chr1:400"], "after the switch, founders C and D")
    }

    @Test
    fun aHeterozygousPanelFounderContributesANoCall(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        runCommand(out)
        val sample1 = genotypesByPosition(out)["sample1"]!!
        // sample1 is (founderC, founderD) over chr1 [200,1000).
        // chr1:600 REF C ALT T. founderC is 1 (T); founderD is het "0/1". Which of founderD's two
        // alleles a descendant inherited is unknown, so that haplotype is not called. Taking the
        // first allele -- what the grits original did -- would give T/C here and T/T for the same
        // unphased genotype written "1/0".
        assertEquals("T/.", sample1["chr1:600"], "the heterozygous founder's haplotype is unknown")
        // chr1:700 REF A ALT C, both founders het: neither haplotype is known.
        assertEquals("./.", sample1["chr1:700"], "both founders heterozygous")
    }

    @Test
    fun aHalfCalledPanelFounderIsReadAsItsOneCalledAllele(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        runCommand(out)
        // chr1:800 REF G ALT A. founderC is "1/." -- one allele called, so there is nothing
        // ambiguous about it and it reads as A. founderD is 0, so sample1 is (A, G).
        assertEquals("A/G", genotypesByPosition(out)["sample1"]!!["chr1:800"])
        // sample3 is haploid founderC, so both its haplotypes take that same called allele.
        assertEquals("A/A", genotypesByPosition(out)["sample3"]!!["chr1:800"])
    }

    @Test
    fun aFounderWithNoCallInThePanelGivesANoCallHaplotype(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        runCommand(out)
        // chr1:900 REF A ALT G. founderC is 0 (A), founderD has no call. sample1 is (C,D), so one
        // haplotype is known and the other is not -- the sample is still reported.
        val genotype = genotypesByPosition(out)["sample1"]!!["chr1:900"]!!
        assertTrue(genotype.contains("A"), "founderC's allele is still called: $genotype")
        assertTrue(genotype.contains("."), "founderD contributes a no-call: $genotype")
    }

    @Test
    fun aSampleWhosePathDoesNotCoverAContigIsMissingThere(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        runCommand(out)
        val sample3 = genotypesByPosition(out)["sample3"]!!
        assertTrue(sample3.containsKey("chr1:10"), "sample3's path covers chr1")
        val chr2 = sample3["chr2:50"]
        assertTrue(chr2 == null || chr2.all { it == '.' || it == '/' || it == '|' },
            "sample3 has no chr2 path, so chr2 must not be called: $chr2")
    }

    @Test
    fun everyPanelSiteAppearsInTheOutput(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        runCommand(out)
        // The panel's sites are the output's sites -- this is what raises density when a dense
        // panel is paired with a path inferred from sparse data.
        fun sitesOf(file: String) = VCFFileReader(File(file), false).use { reader ->
            reader.map { "${it.contig}:${it.start}" }
        }
        assertEquals(sitesOf(panel), sitesOf(out.path))
    }

    @Test
    fun theSampleColumnsAreThePathSamplesNotThePanelFounders(@TempDir dir: File) {
        val out = File(dir, "imputed.vcf")
        runCommand(out)
        VCFFileReader(out, false).use { reader ->
            assertEquals(listOf("sample1", "sample2", "sample3"),
                reader.fileHeader.sampleNamesInOrder.sorted())
        }
    }

    @Test
    fun testCliktParams() {
        val missingBedDir = BedToVcf().test("--reference-panel-vcf $panel --output-file out.vcf")
        assertEquals(1, missingBedDir.statusCode)
        assertTrue(missingBedDir.stderr.contains("missing option --bed-dir"), missingBedDir.stderr)

        val missingPanel = BedToVcf().test("--bed-dir $testDir --output-file out.vcf")
        assertEquals(1, missingPanel.statusCode)
        assertTrue(missingPanel.stderr.contains("missing option --reference-panel-vcf"), missingPanel.stderr)

        val missingOut = BedToVcf().test("--bed-dir $testDir --reference-panel-vcf $panel")
        assertEquals(1, missingOut.statusCode)
        assertTrue(missingOut.stderr.contains("missing option --output-file"), missingOut.stderr)

        val badBedDir = BedToVcf().test(
            "--bed-dir $testDir/nope --reference-panel-vcf $panel --output-file out.vcf")
        assertEquals(1, badBedDir.statusCode)
        assertTrue(badBedDir.stderr.contains("not a directory"), badBedDir.stderr)
    }

    @Test
    fun anEmptyBedDirectoryIsReportedRatherThanWritingAnEmptyVcf(@TempDir dir: File) {
        val result = BedToVcf().test(
            "--bed-dir ${dir.absolutePath} --reference-panel-vcf $panel --output-file ${dir.path}/out.vcf")
        assertEquals(1, result.statusCode)
        assertTrue(result.stderr.contains("No usable BED files"), result.stderr)
    }
}
