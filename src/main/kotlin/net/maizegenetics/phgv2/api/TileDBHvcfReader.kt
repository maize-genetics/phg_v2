package net.maizegenetics.phgv2.api

import io.tiledb.libvcfnative.VCFReader
import java.util.*

class TileDBHvcfReader(uri: String, samples: List<String>?, ranges: List<HvcfReader.PositionRange>?): HvcfReader {

    val dbReader: VCFReader

    init {
        if (samples == null) dbReader = VCFReader(uri, null, Optional.empty(), Optional.empty())
        else {
            dbReader = VCFReader(uri, samples.toTypedArray(), Optional.empty(), Optional.empty())
        }
    }

    override fun range(range: List<HvcfReader.PositionRange>): HvcfReader {
        val rangeString = range.map { it.toString() }.toTypedArray()
        dbReader.setRanges(rangeString)
        return this
    }

    override fun header(sampleName: String): String {
        TODO("Not yet implemented")
    }

    override fun headerMap(): Map<String, String> {
        TODO("Not yet implemented")
    }

    override fun sampleNames(): List<String> {
        TODO("Not yet implemented")
    }

    override fun data(): Map<String, List<HvcfReader.SampleData>> {
        TODO("Not yet implemented")
    }

}