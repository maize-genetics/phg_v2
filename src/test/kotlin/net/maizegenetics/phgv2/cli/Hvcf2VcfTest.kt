package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.Position
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.fail

class Hvcf2VcfTest {

    @Test
    fun testCliktParams() {
        val hvcf2Vcf = Hvcf2Vcf()

        //check parameters
        //val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        //        .default("")
        //
        //
        //    val hvcfDir by option(help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        //        .required()
        //
        //    val donorVcfFile by option(help = "Path to the VCF file containing all the PHG SNPs.  This is typically created by running merge-gvcf on the Assembly Gvcf files.")
        //        .required()
        //
        //    val outputFile by option(help = "Output file.")
        //        .required()
        //
        //    //This is needed to create the refSeqDictionary
        //    val referenceFile by option(help = "Path to local Reference FASTA file needed for sequence dictionary")
        //        .required()

        //"--hvcf-dir /path/to/hvcf --donor-vcf-file /path/to/donor.vcf --output-file /path/to/output.vcf --reference-file /path/to/reference.fasta"

        val noHvcfDir = hvcf2Vcf.test("--donor-vcf-file /path/to/donor.vcf --output-file /path/to/output.vcf --reference-file /path/to/reference.fasta")
        assertEquals(noHvcfDir.statusCode, 1)
        assertEquals(
            "Usage: hvcf2vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --hvcf-dir\n", noHvcfDir.output
        )

        val noDonorVcf = hvcf2Vcf.test("--hvcf-dir /path/to/hvcf --output-file /path/to/output.vcf --reference-file /path/to/reference.fasta")
        assertEquals(noDonorVcf.statusCode, 1)
        assertEquals(
            "Usage: hvcf2vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --donor-vcf-file\n", noDonorVcf.output
        )

        val noOutputFile = hvcf2Vcf.test("--hvcf-dir /path/to/hvcf --donor-vcf-file /path/to/donor.vcf --reference-file /path/to/reference.fasta")
        assertEquals(noOutputFile.statusCode, 1)
        assertEquals(
            "Usage: hvcf2vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --output-file\n", noOutputFile.output
        )

        val noRefFile = hvcf2Vcf.test("--hvcf-dir /path/to/hvcf --donor-vcf-file /path/to/donor.vcf --output-file /path/to/output.vcf")
        assertEquals(noRefFile.statusCode, 1)
        assertEquals(
            "Usage: hvcf2vcf [<options>]\n" +
                    "\n" +
                    "Error: missing option --reference-file\n", noRefFile.output
        )
    }

    @Test
    fun processHVCFAndBuildVCFTest() {
        fail("Not yet implemented")
    }

    @Test
    fun createRangeHapMapToSampleGameteTest() {
        fail("Not yet implemented")
    }

    @Test
    fun processSingleHvcfFileTest() {
        fail("Not yet implemented")
    }

    @Test
    fun processSingleHvcfVariantTest() {
        fail("Not yet implemented")
    }

    @Test
    fun buildPositionToRefRangeMapTest() {
        //buildPositionToRefRangeMap(ranges: Set<ReferenceRange>): TreeMap<Position, ReferenceRange> {

        val hvcf2Vcf = Hvcf2Vcf()

        //Build a simple reference range set
        val ranges = setOf(ReferenceRange("chr1", 1, 100), ReferenceRange("chr1", 101, 200), ReferenceRange("chr2", 1, 100))

        val treeMap = hvcf2Vcf.buildPositionToRefRangeMap(ranges)
        assertEquals(6, treeMap.size)
        assertEquals(ReferenceRange("chr1", 1, 100), treeMap[Position("chr1", 1)])
        assertEquals(null, treeMap[Position("chr1", 50)])
        assertEquals(ReferenceRange("chr1", 1, 100), treeMap.floorEntry(Position("chr1", 50)).value)
        assertEquals(ReferenceRange("chr1", 101, 200), treeMap[Position("chr1", 101)])
        assertEquals(null, treeMap[Position("chr1", 150)])
        assertEquals(ReferenceRange("chr1", 101, 200), treeMap.floorEntry(Position("chr1", 150)).value)
        assertEquals(ReferenceRange("chr2", 1, 100), treeMap[Position("chr2", 1)])
        assertEquals(null, treeMap[Position("chr2", 50)])
        assertEquals(ReferenceRange("chr2", 1, 100), treeMap.floorEntry(Position("chr2", 50)).value)

    }

    @Test
    fun extractVcfAndExportTest() {
        fail("Not yet implemented")
    }

    @Test
    fun extractVariantsAndWriteTest() {
        fail("Not yet implemented")
    }

}