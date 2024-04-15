package net.maizegenetics.phgv2.rphg

import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
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
     * @param nThreads Number of threads to use
     *
     * @return A `StringMatrix` object that will be converted to a
     *         `character` matrix in R
     */
    fun getHapIdMatrixFromGraph(hapGraph: HaplotypeGraph, nThreads: Int = 100): RCharacterMatrix {
        val allRanges = hapGraph.ranges()
        val sampleGametes = hapGraph.sampleGametesInGraph()
        val sampleChannel = Channel<Deferred<Array<String>>>(nThreads)

        CoroutineScope(Dispatchers.IO).launch {
            sampleGametes
                .forEach {
                    sampleChannel.send(async { hapGraph.sampleGameteToHaplotypeId(it) })
                }
            sampleChannel.close()
        }

        return runBlocking {
            val result = mutableListOf<Array<String>>()
            for (deferred in sampleChannel) {
                result.add(deferred.await())
            }

            RCharacterMatrix (
                colNames = allRanges.map { it.toString() }.toTypedArray(),
                rowNames = sampleGametes.map { "${it.name}_G${it.gameteId + 1}" }.toTypedArray(),
                matrixData = result.toTypedArray()
            )
        }
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
            colNames = arrayOf("seqname", "start", "end", "rr_id"),
            rowNames = null,
            matrixData = arrayOf(
                refRanges.map { it.contig }.toTypedArray(),
                refRanges.map { it.start }.toTypedArray(),
                refRanges.map { it.end }.toTypedArray(),
                refRanges.map { it.toString() }.toTypedArray()
            )
        )
    }

    /**
     * Return haplotype ID positional data from a
     * HaplotypeGraph object's AltHeader metadata
     *
     * @param hapGraph A `HaplotypeGraph` object
     *
     * @return A `RPHGList` object that will be converted to a
     *         `data.frame` object in R
     */
    fun getAltHeaderPositionsFromGraph(hapGraph: HaplotypeGraph): RList {
        val hapIdList = mutableListOf<String>()
        val cOneId = mutableListOf<String>()
        val cTwoId = mutableListOf<String>()
        val cOnePos = mutableListOf<Int>()
        val cTwoPos = mutableListOf<Int>()

        val altHeaders = hapGraph.altHeaders()

        altHeaders.forEach { (key, valueRegion) ->
            hapIdList.addAll(List(valueRegion.regions.size) { key })

            valueRegion.regions.forEach {
                cOneId.add(it.first.contig)
                cTwoId.add(it.second.contig)
                cOnePos.add(it.first.position)
                cTwoPos.add(it.second.position)
            }
        }

        return RList(
            colNames = arrayOf(
                "hap_id",
                "contig_start",
                "contig_end",
                "start",
                "end"
            ),
            matrixData = arrayOf(
                hapIdList.toTypedArray(),
                cOneId.toTypedArray(),
                cTwoId.toTypedArray(),
                cOnePos.toTypedArray(),
                cTwoPos.toTypedArray()
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

        return RList(
            colNames = arrayOf(
                "hap_id",
                "sample_name",
                "description",
                "source",
                "checksum",
                "ref_range_hash"
            ),
            matrixData = arrayOf(
                altHeaders.map { it.value.id }.toTypedArray(),
                altHeaders.map { it.value.sampleName() }.toTypedArray(),
                altHeaders.map { it.value.description }.toTypedArray(),
                altHeaders.map { it.value.source }.toTypedArray(),
                altHeaders.map { it.value.checksum }.toTypedArray(),
                altHeaders.map { it.value.refRange }.toTypedArray()
            )
        )
    }
}

