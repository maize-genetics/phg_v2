package net.maizegenetics.phgv2.rphg

import net.maizegenetics.phgv2.api.HaplotypeGraph

/**
 * Main class for interfacing PHGv2 data return methods with R
 */
class RMethods {
    /**
     * Return a 2d haplotype ID matrix (sample gametes by ref range)
     * from a HaplotypeGraph object
     *
     * @param hapGraph A `HaplotypeGraph` object
     *
     * @return A `StringMatrix` object that will be converted to a
     *         `character` matrix in R
     */
    fun getHapIdMatrixFromGraph(hapGraph: HaplotypeGraph): RCharacterMatrix {
        val allRanges = hapGraph.ranges()
        val sampleGametes = hapGraph.sampleGametesInGraph()
        val array2D = Array(hapGraph.samples().size) { Array(allRanges.size) { "" } }

        allRanges.forEach { range ->
            val hapIdRange = hapGraph.sampleGameteToHaplotypeId(range)
            hapIdRange.forEach { sgEntry ->
                array2D[sampleGametes.indexOf(sgEntry.key)][allRanges.indexOf(range)] = sgEntry.value
            }
        }

        return RCharacterMatrix (
            colNames = allRanges.map { "R$it" }.toTypedArray(),
            rowNames = sampleGametes.map { "${it.name}_G${it.gameteId + 1}" }.toTypedArray(),
            matrixData = array2D
        )
    }

    /**
     * Return reference range data (contig, start, end) from a
     * HaplotypeGraph object
     *
     * @param hapGraph A `HaplotypeGraph` object
     *
     * @return A `RPHGList` object that will be converted to a
     *         `GenomicRanges` object in R
     */
    fun getRefRangesFromGraph(hapGraph: HaplotypeGraph): RList {
        val refRanges = hapGraph.ranges()
        return RList(
            colNames = arrayOf("seqname", "start", "end"),
            rowNames = null,
            matrixData = arrayOf(
                refRanges.map { it.contig }.toTypedArray(),
                refRanges.map { it.start }.toTypedArray(),
                refRanges.map { it.end }.toTypedArray()
            )
        )
    }

    /**
     * Return haplotype ID metadata (alt headers) from a
     * HaplotypeGraph object
     *
     * @param hapGraph A `HaplotypeGraph` object
     *
     * @return A `RPHGList` object that will be converted to a
     *         `data.frame` object in R
     */
    fun getAltHeadersFromGraph(hapGraph: HaplotypeGraph): RList {
        val altHeaders = hapGraph.altHeaders()

        // Return all possible positions for each hash - returns as
        // a list of `RList` objects (e.g. a list of dataframe-like
        // objects when evaluated)
        val positions = altHeaders
            .map { it.value.regions }
            .map { posPair ->
                RList(
                    colNames = arrayOf(
                        "contig_start",
                        "contig_end",
                        "start",
                        "end"
                    ),
                    matrixData = arrayOf(
                        posPair.map { it.first.contig }.toTypedArray(),
                        posPair.map { it.second.contig }.toTypedArray(),
                        posPair.map { it.first.position }.toTypedArray(),
                        posPair.map { it.second.position }.toTypedArray()
                    )
                )
            }

        return RList(
            colNames = arrayOf(
                "hap_id",
                "sample_name",
                "description",
                "source",
                "checksum",
                "positions",
                "ref_range_hash"
            ),
            matrixData = arrayOf(
                altHeaders.map { it.value.id }.toTypedArray(),
                altHeaders.map { it.value.sampleName }.toTypedArray(),
                altHeaders.map { it.value.description }.toTypedArray(),
                altHeaders.map { it.value.source }.toTypedArray(),
                altHeaders.map { it.value.checksum }.toTypedArray(),
                positions.toTypedArray(),
                altHeaders.map { it.value.refRange }.toTypedArray()
            )
        )
    }
}

