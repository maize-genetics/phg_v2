package net.maizegenetics.phgv2.api

import io.tiledb.java.api.Context
import io.tiledb.java.api.QueryType
import io.tiledb.libvcfnative.VCFReader
import java.nio.ByteBuffer
import java.nio.ByteOrder
import java.util.*

class TileDBHvcfReader(val uri: String, samples: List<String>? = null, ranges: List<HvcfReader.PositionRange>? = null): HvcfReader {

    val dbReader: VCFReader
    private val missingInt = -1

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
        Context().use {context ->
            io.tiledb.java.api.Array(context, uri, QueryType.TILEDB_READ).use {

            }
        }

        return listOf()
    }

    override fun data(): List<HvcfReader.SampleData> {
        val capacity = 1024
        val attributeNames = arrayOf(
            "sample_name",
            "contig",
            "pos_start",
            "pos_end",
            "alleles",
            "fmt_GT",
            "fmt_AD",
            "fmt_DP"
        )

        //delete any attributes not in this tiledb
        val attributeInfo = dbReader.attributes
        val filteredAttributeNames = attributeNames.filter { attributeInfo.keys.contains(it) }
        val hasAD = filteredAttributeNames.contains("fmt_AD")
        val hasDP = filteredAttributeNames.contains("fmt_DP")

        //set buffers
        //buffers for all candidate attributes that are in the database
        for (attribute in filteredAttributeNames) dbReader.setBuffer(attribute,
            ByteBuffer.allocateDirect(capacity).order(ByteOrder.nativeOrder()))

        //needed for variable length attributes
        attributeNames.filter { attributeInfo[it]?.isVarLen ?: false }.forEach {
            dbReader.setBufferOffsets(it, ByteBuffer.allocateDirect(capacity).order(ByteOrder.nativeOrder()))
        }

        //needed for nullable attributes
        attributeNames.filter { attributeInfo[it]?.isNullable ?: false }.forEach {
            dbReader.setBufferValidityBitmap(it, ByteBuffer.allocateDirect(capacity).order(ByteOrder.nativeOrder()))
        }

        //needed for list attributes
        attributeNames.filter { attributeInfo[it]?.isList ?: false }.forEach {
            dbReader.setBufferListOffsets(it, ByteBuffer.allocateDirect(capacity).order(ByteOrder.nativeOrder()))
        }

        //process data
        val sampleNames = mutableListOf<String>()
        val contigs = mutableListOf<String>()
        val startPositions = mutableListOf<Int>()
        val endPositions = mutableListOf<Int>()
        val genotypes = mutableListOf<List<String>>()
        val ADs = mutableListOf<List<Int>>()
        val DPs = mutableListOf<Int>()

        while (!dbReader.status.equals(VCFReader.Status.COMPLETED)) {
            dbReader.submit()
            val numberOfRecords = dbReader.numRecords.toInt()

            //get sample names
            decodeVarlenChar(
                dbReader.getBuffer("sample_name"),
                dbReader.getOffsets("sample_name"),
                numberOfRecords,
                sampleNames
            )

            //get contig (contig, varlen, CHAR)
            decodeVarlenChar(
                dbReader.getBuffer("contig"),
                dbReader.getOffsets("contig"),
                numberOfRecords,
                contigs,
            )

            //get startPos (pos_start, INT32)
            decodeInt(
                dbReader.getBuffer("pos_start"),
                numberOfRecords,
                startPositions
            )

            //get endPos (pos_end, INT32)
            decodeInt(
                dbReader.getBuffer("pos_end"),
                numberOfRecords,
                endPositions
            )

            //get genotype (alleles, varlen list of CHAR; fmt_GT, varlen, nullable, INT32)
            //get alleles and fmt_GT, then derive genotypes
            val gtList = mutableListOf<List<Int>>()
            val alleles = mutableListOf<List<String>>()

            decodeVarlenNullableInt(dbReader.getBuffer("fmt_GT"),
                dbReader.getOffsets("fmt_GT"),
                BitSet.valueOf(dbReader.getBitMap("fmt_GT")),
                numberOfRecords,
                gtList
            )

            decodeVarlenListChar(dbReader.getBuffer("alleles"),
                dbReader.getOffsets("alleles"),
                dbReader.getListOffsets("alleles"),
                numberOfRecords,
                alleles
            )

            for (ndx in 0 until numberOfRecords) {
                val currAlleles = alleles[ndx]
                val geno = gtList[ndx].map { if( it < 0) "." else currAlleles[it] }
                genotypes.add(geno)
            }

            //get AD (fmt_AD, varlen, nullable, INT32)
            if (hasAD) {
                decodeVarlenNullableInt(dbReader.getBuffer("fmt_AD"),
                    dbReader.getOffsets("fmt_AD"),
                    BitSet.valueOf(dbReader.getBitMap("fmt_AD")),
                    numberOfRecords, ADs
                )
            }

            //get DP (fmt_DP, nullable, INT32)
            if (hasDP) {
                decodeNullableInt(dbReader.getBuffer("fmt_DP"),
                    BitSet.valueOf(dbReader.getBitMap("fmt_DP")),
                    numberOfRecords,
                    DPs
                )
            }

        }

        //convert values to a list of sample data
        return sampleNames.indices.map { HvcfReader.SampleData(sampleNames[it], contigs[it],
            startPositions[it], endPositions[it], genotypes[it],
            if (hasAD) ADs[it] else null, if (hasDP) DPs[it] else null) }

    }

    fun decodeVarlenChar(dataBuffer: ByteBuffer,
                         offsetBuffer: ByteBuffer,
                         numberOfRecords: Int,
                         outputList: MutableList<String>)  {
        val offsets = IntArray(numberOfRecords + 1)
        offsetBuffer.asIntBuffer().get(offsets)
        for (ndx in 0..<numberOfRecords) {
            val sizeOfName = offsets[ndx + 1] - offsets[ndx]
            val nameByteArray = ByteArray(sizeOfName)
            dataBuffer.get(nameByteArray)
            outputList.add(nameByteArray.decodeToString())
        }
    }

    fun decodeInt(dataBuffer: ByteBuffer,
                  numberOfRecords: Int,
                  outputList: MutableList<Int>) {
        val data = IntArray(numberOfRecords)
        dataBuffer.asIntBuffer().get(data)
        outputList.addAll(data.toList())
    }

    fun decodeVarlenListChar(dataBuffer: ByteBuffer,
                             offsetBuffer: ByteBuffer,
                             listOffsetBuffer: ByteBuffer,
                             numberOfRecords: Int,
                             outputList: MutableList<List<String>>) {
        val nameList = mutableListOf<String>()
        decodeVarlenChar(dataBuffer, offsetBuffer, numberOfRecords, nameList)
        val listOffsets = IntArray(numberOfRecords + 1)
        listOffsetBuffer.asIntBuffer().get(listOffsets)
        val numberOfValues = listOffsets[numberOfRecords]
        decodeVarlenChar(dataBuffer, offsetBuffer, numberOfValues, nameList)
        for (ndx in 0 until numberOfRecords) {
            outputList.add(nameList.subList(listOffsets[ndx], listOffsets[ndx + 1]))
        }
    }

    fun decodeVarlenNullableInt(dataBuffer: ByteBuffer,
                                offsetBuffer: ByteBuffer,
                                validityBitset: BitSet,
                                numberOfRecords: Int,
                                outputList: MutableList<List<Int>>) {
        val offsets = IntArray(numberOfRecords + 1)
        offsetBuffer.asIntBuffer().get(offsets)
        val data = IntArray(offsets[numberOfRecords])
        dataBuffer.asIntBuffer().get(data)
        for (offsetIndex in 0..<numberOfRecords) {
            val obsList = mutableListOf<Int>()
            outputList.add(obsList)
            for (ndx in offsets[offsetIndex] until offsets[offsetIndex + 1]) {
                obsList.add(
                    if (validityBitset.get(offsetIndex)) data[ndx] else missingInt
                )
            }

        }

    }

    fun decodeNullableInt(dataBuffer: ByteBuffer,
                          validityBitset: BitSet,
                          numberOfRecords: Int,
                          outputList: MutableList<Int>) {
        val data = IntArray(numberOfRecords)
        dataBuffer.asIntBuffer().get(data)
        for (ndx in 0 until numberOfRecords) {
            outputList.add(
                if ( validityBitset.get(ndx) ) data[ndx] else missingInt
            )
        }

    }
}