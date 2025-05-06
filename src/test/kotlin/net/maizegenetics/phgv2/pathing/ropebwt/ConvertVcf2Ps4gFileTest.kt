package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.utils.Position
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.fail
import org.junit.jupiter.api.Test
import kotlin.test.junit5.JUnit5Asserter.fail

class ConvertVcf2Ps4gFileTest {

    @Test
    fun testCliktParams() {
        val convertVcf2Ps4gFile = ConvertVcf2Ps4gFile()

        // val sampleVcf by option(help = "Sample VCF file")
        //        .required()
        //
        //    val gameteVcf by option(help = "Gamete VCF file")
        //        .required()
        //
        //    val outputDir by option(help = "Output file")
        //        .required()

        //val noVcfDir = buildSplineKnots.test("--output-file dummyFile.json.gz")
        //        assertEquals(1, noVcfDir.statusCode)
        //        assertEquals("Usage: build-spline-knots [<options>]\n\n" +
        //                "Error: missing option --vcf-dir\n", noVcfDir.stderr)

        val noSampleVCF = convertVcf2Ps4gFile.test("--gamete-vcf testDir.vcf --output-dir testDir")
        assertEquals(1, noSampleVCF.statusCode)
        assertEquals("Usage: convert-vcf2ps4g-file [<options>]\n\n" +
                "Error: missing option --sample-vcf\n", noSampleVCF.stderr)

        val noGameteVCF = convertVcf2Ps4gFile.test("--sample-vcf testDir.vcf --output-dir testDir")
        assertEquals(1, noGameteVCF.statusCode)
        assertEquals("Usage: convert-vcf2ps4g-file [<options>]\n\n" +
                "Error: missing option --gamete-vcf\n", noGameteVCF.stderr)

        val noOutputDir = convertVcf2Ps4gFile.test("--sample-vcf testDir.vcf --gamete-vcf testDir.vcf")
        assertEquals(1, noOutputDir.statusCode)
        assertEquals("Usage: convert-vcf2ps4g-file [<options>]\n\n" +
                "Error: missing option --output-dir\n", noOutputDir.stderr)
    }


    @Test
    fun testCreatePositionSampleGameteLookup() {
        fail("Not yet implemented")
    }

    @Test
    fun testCreatePS4GData() {
        fail("Not yet implemented")
    }

    //fun createGameteToIdxMap(positionSampleGameteLookup: Map<Position, Map<String, List<SampleGamete>>>): Map<SampleGamete, Int>
    @Test
    fun testCreateGameteToIdxMap() {
        val convertVcf2Ps4gFile = ConvertVcf2Ps4gFile()

        val positionSampleGameteLookup = mapOf<Position, Map<String, List<SampleGamete>>>(
            Position("chr1", 1) to mapOf(
                "imputeSample1" to listOf(SampleGamete("sample1", 0)),
                "imputeSample2" to listOf(SampleGamete("sample2", 0))),
            Position("chr1", 2) to mapOf(
                "imputeSample1" to listOf(SampleGamete("sample1", 0), SampleGamete("sample1", 1)),
                "imputeSample2" to listOf(SampleGamete("sample1", 1))),
            Position("chr1", 3) to mapOf(
                "imputeSample1" to listOf(SampleGamete("sample1", 0), SampleGamete("sample1", 1)),
                "imputeSample2" to listOf(SampleGamete("sample4", 0), SampleGamete("sample3", 0))),
        )

        val truthMap = mapOf<SampleGamete, Int>(
            SampleGamete("sample1", 0) to 0,
            SampleGamete("sample1", 1) to 1,
            SampleGamete("sample2", 0) to 2,
            SampleGamete("sample4", 0) to 4,
            SampleGamete("sample3", 0) to 3
        )



        val gameteToIdxMap = convertVcf2Ps4gFile.createGameteToIdxMap(positionSampleGameteLookup)
        assertEquals(truthMap.size, gameteToIdxMap.size)
        truthMap.forEach { (sampleGamete, idx) ->
            assertEquals(idx, gameteToIdxMap[sampleGamete])
        }

    }

    @Test
    fun testConvertCountMapsToData() {
        fail("Not yet implemented")
    }



}