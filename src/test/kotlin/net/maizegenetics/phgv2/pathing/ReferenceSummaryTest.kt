package net.maizegenetics.phgv2.pathing

import biokotlin.seqIO.NucSeqIO
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.sequenceDictionary
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import java.io.File
import java.util.zip.GZIPOutputStream
import kotlin.test.assertEquals

class ReferenceSummaryTest {

    private val ranges = HaplotypeGraph(listOf(TestExtension.smallseqRefHvcfFile)).ranges()

    /** Streaming the fasta gives the same lengths, contig order and range start bases as NucSeqIO. */
    @Test
    fun matchesNucSeqIO() {
        val sequences = NucSeqIO(TestExtension.smallseqRefFile).readAll()
        val summary = ReferenceSummary.read(TestExtension.smallseqRefFile, ranges)

        assertEquals(sequences.keys.toList(), summary.contigLengths.keys.toList())
        assertEquals(sequenceDictionary(sequences).sequences.map { it.sequenceName to it.sequenceLength },
            summary.sequenceDictionary().sequences.map { it.sequenceName to it.sequenceLength })
        for (range in ranges) {
            assertEquals(sequences[range.contig]!!.sequence[range.start - 1].name, summary.startBase(range), "$range")
        }
    }

    /**
     * Lower case, header words after the name, uneven line lengths, a range starting on the first base of a
     * line, and a gzipped file.
     */
    @Test
    fun handlesFastaDetails() {
        val fasta = ">c1 description words\nacgT\nGG\nttAc\n>c2\nNNA\n"
        val plain = File.createTempFile("ref", ".fa").also { it.deleteOnExit(); it.writeText(fasta) }
        val gzipped = File.createTempFile("ref", ".fa.gz").also { file ->
            file.deleteOnExit()
            GZIPOutputStream(file.outputStream()).bufferedWriter().use { it.write(fasta) }
        }
        val testRanges = listOf(ReferenceRange("c1", 1, 3), ReferenceRange("c1", 5, 6),
            ReferenceRange("c1", 7, 10), ReferenceRange("c2", 3, 3))
        for (file in listOf(plain, gzipped)) {
            val summary = ReferenceSummary.read(file.path, testRanges)
            assertEquals(mapOf("c1" to 10, "c2" to 3), summary.contigLengths)
            assertEquals(listOf("A", "G", "T", "A"), testRanges.map { summary.startBase(it) })
            assertThrows<IllegalArgumentException> { summary.startBase(ReferenceRange("c3", 1, 2)) }
        }
    }
}
