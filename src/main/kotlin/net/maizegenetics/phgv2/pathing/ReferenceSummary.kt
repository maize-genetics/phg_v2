package net.maizegenetics.phgv2.pathing

import biokotlin.seqIO.NucSeqIO
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.SAMSequenceRecord
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.utils.getBufferedReader

/**
 * The parts of a reference genome that [FindPaths] needs to write an h.vcf: each contig's length, for the
 * header's sequence dictionary, and the reference base at the start of each reference range, for the REF
 * allele. [read] gets them in one pass over a fasta file, so the sequence is never held in memory and the
 * reference is read once per run rather than once per sample.
 *
 * @property contigLengths contig name to length, in the order the contigs appear in the reference.
 */
class ReferenceSummary(val contigLengths: Map<String, Int>, private val startBases: Map<ReferenceRange, String>) {

    /** The reference base at [range]'s start position, upper case. */
    fun startBase(range: ReferenceRange): String = startBases[range]
        ?: throw IllegalArgumentException("The reference genome has no base at ${range.contig}:${range.start}")

    fun sequenceDictionary() = SAMSequenceDictionary(contigLengths.map { (name, length) -> SAMSequenceRecord(name, length) })

    companion object {

        /**
         * Reads [referenceFile] (fasta, optionally gzipped) for [ranges]. A contig's name is the first
         * word of its header line. A file that is not fasta (one that does not start with '>') is read
         * whole with NucSeqIO instead.
         */
        fun read(referenceFile: String, ranges: Collection<ReferenceRange>): ReferenceSummary {
            val firstChar = getBufferedReader(referenceFile).use { reader ->
                generateSequence { reader.readLine() }.firstOrNull { it.isNotBlank() }?.trimStart()?.firstOrNull()
            }
            return if (firstChar == '>') readFasta(referenceFile, ranges) else readWithNucSeqIO(referenceFile, ranges)
        }

        private fun readFasta(referenceFile: String, ranges: Collection<ReferenceRange>): ReferenceSummary {
            val rangesByContig = ranges.groupBy { it.contig }.mapValues { (_, list) -> list.sortedBy { it.start } }
            val contigLengths = LinkedHashMap<String, Int>()
            val startBases = HashMap<ReferenceRange, String>()

            var contig: String? = null
            var length = 0L
            var contigRanges = emptyList<ReferenceRange>()
            var next = 0

            fun endContig() {
                val name = contig ?: return
                check(length <= Int.MAX_VALUE) { "Contig $name is longer than ${Int.MAX_VALUE} bases" }
                contigLengths[name] = length.toInt()
            }

            getBufferedReader(referenceFile).use { reader ->
                var line = reader.readLine()
                while (line != null) {
                    if (line.startsWith(">")) {
                        endContig()
                        contig = line.substring(1).trim().split(Regex("\\s+"))[0]
                        length = 0L
                        contigRanges = rangesByContig[contig] ?: emptyList()
                        next = 0
                    } else {
                        val sequence = line.trim()
                        val lineEnd = length + sequence.length
                        // a range's 1-based start lies on this line when length < start <= lineEnd
                        while (next < contigRanges.size && contigRanges[next].start <= lineEnd) {
                            val range = contigRanges[next++]
                            if (range.start > length) {
                                startBases[range] = sequence[(range.start - 1 - length).toInt()].uppercaseChar().toString()
                            }
                        }
                        length = lineEnd
                    }
                    line = reader.readLine()
                }
            }
            endContig()
            return ReferenceSummary(contigLengths, startBases)
        }

        private fun readWithNucSeqIO(referenceFile: String, ranges: Collection<ReferenceRange>): ReferenceSummary {
            val sequences = NucSeqIO(referenceFile).readAll()
            val contigLengths = LinkedHashMap<String, Int>()
            for ((name, seq) in sequences) contigLengths[name] = seq.size()
            val startBases = ranges.filter { sequences.containsKey(it.contig) && it.start <= sequences[it.contig]!!.size() }
                .associateWith { sequences[it.contig]!!.sequence[it.start - 1].name }
            return ReferenceSummary(contigLengths, startBases)
        }
    }
}
