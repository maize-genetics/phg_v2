package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.brapi.createSmallSeqTiledb
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.getChecksum
import net.maizegenetics.phgv2.utils.seqFromAGC
import org.apache.logging.log4j.LogManager
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class CreateFastaFromHvcfRangeFastaTest {

    private val myLogger = LogManager.getLogger(CreateFastaFromHvcfRangeFastaTest::class.java)

    companion object {

        private val exportHvcfDir = "${TestExtension.tempDir}/ref-range-fasta/"
        private val multiInputDir = "${TestExtension.tempDir}/multi-input/"
        private val outputHvcfDir = "${exportHvcfDir}/output/"

        private val bedfile = "data/test/refRangeFasta/refRangeFasta.bed"
        private val expectedFasta1 = "data/test/refRangeFasta/testRefRangeFasta.vcf-1_1-1000.fasta"
        private val expectedFasta2 = "data/test/refRangeFasta/testRefRangeFasta.vcf-2_23001-27500.fasta"

        private val dbPath = "${exportHvcfDir}/tiledb_ref_range_fasta"

        val HVCF_PATTERN = Regex("""(\.hvcf|\.h\.vcf|\.hvcf\.gz|\.h\.vcf\.gz)$""")

        @BeforeAll
        @JvmStatic
        fun setup() {
            File(exportHvcfDir).mkdirs()
            File(outputHvcfDir).mkdirs()
            File(multiInputDir).mkdirs()

            File(TestExtension.smallseqRefHvcfFile).copyTo(File(multiInputDir + File(TestExtension.smallseqRefHvcfFile).name))
            File(TestExtension.smallseqLineAHvcfFile).copyTo(File(multiInputDir + File(TestExtension.smallseqLineAHvcfFile).name))
            File(TestExtension.smallseqLineBHvcfFile).copyTo(File(multiInputDir + File(TestExtension.smallseqLineBHvcfFile).name))

            createSmallSeqTiledb(dbPath)
        }

    }

    @Test
    fun testRefRangeFasta() {

        val outputFile = "${outputHvcfDir}testRefRangeFasta.vcf"

        // phg create-fasta-from-hvcf --db-path /workdir/wl748/db_huehue
        // --hvcf-dir /local/workdir/wl748/huehue_hvcf_v2 --output-dir ref-range-fasta/Terry
        // --range-bedfile hpc1_RefRange.bed --fasta-type rangeFasta
        val result = CreateFastaFromHvcf().test(
            "--db-path $dbPath --hvcf-dir $multiInputDir --output-dir $outputFile --range-bedfile $bedfile --fasta-type rangeFasta"
        )

        myLogger.info("testRefRangeFasta: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        var checksum1 = getChecksum(expectedFasta1)
        var checksum2 = getChecksum("${outputHvcfDir}testRefRangeFasta.vcf-1_1-1000.fasta")

        myLogger.info("testRefRangeFasta1 expected checksum1: $checksum1")
        myLogger.info("testRefRangeFasta1 actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testRefRangeFasta1 checksums do not match")

        checksum1 = getChecksum(expectedFasta2)
        checksum2 = getChecksum("${outputHvcfDir}testRefRangeFasta.vcf-2_23001-27500.fasta")

        myLogger.info("testRefRangeFasta2 expected checksum1: $checksum1")
        myLogger.info("testRefRangeFasta2 actual checksum2: $checksum2")

        assertEquals(checksum1, checksum2, "testRefRangeFasta2 checksums do not match")

    }

    @Test
    fun testHaplotypeGraphHapidToSeqLength() {

        val graph = HaplotypeGraph(multiInputDir)

        val hapidToSeqLength = graph.hapidToSeqLength()

        graph.ranges()
            .forEach { range ->
                graph.hapIdToSampleGametes(range).keys
                    .forEach { hapid ->
                        seqFromAGC(
                            dbPath,
                            graph,
                            hapid,
                            Pair(Position(range.contig, range.start), Position(range.contig, range.end))
                        )
                            .let { (seq, _) ->
                                assertEquals(
                                    hapidToSeqLength.getValue(hapid),
                                    seq.length,
                                    "hapid: $hapid seq length: ${seq.length} != ${hapidToSeqLength.getValue(hapid)}"
                                )
                            }
                    }
            }

    }

}