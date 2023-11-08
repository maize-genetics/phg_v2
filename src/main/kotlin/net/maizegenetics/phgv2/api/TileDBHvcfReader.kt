package net.maizegenetics.phgv2.api

class TileDBHvcfReader(uri: String): HvcfReader {


    init {

    }

    override fun range(range: List<HvcfReader.PositionRange>): HvcfReader {
        TODO("Not yet implemented")
    }

    override fun samples(sampleNames: List<String>): HvcfReader {
        TODO("Not yet implemented")
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