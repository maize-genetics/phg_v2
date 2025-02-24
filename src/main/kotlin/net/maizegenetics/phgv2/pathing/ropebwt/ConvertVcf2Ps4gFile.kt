package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.cli.headerCommand
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.utils.Position
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * NOTE This class is currently a WIP.  This will need to be finished eventually but was started as a hackathon project that got pivoted away from.
 */
data class SingleSamplePS4GData(val ps4gData: List<PS4GData>, val gameteCountMap: Map<SampleGamete,Int>)
data class VCFPS4GData(val knownGameteToIdxMap: Map<SampleGamete,Int>, val sampleToPS4GMap : Map<String, SingleSamplePS4GData>)

class ConvertVcf2Ps4gFile: CliktCommand(help = "Convert VCF to PS4G") {

    private val myLogger = LogManager.getLogger(ConvertVcf2Ps4gFile::class.java)

    val sampleVcf by option(help = "Sample VCF file")
        .required()

    val gameteVcf by option(help = "Gamete VCF file")
        .required()

    val outputDir by option(help = "Output file")
        .required()


    override fun run() {
        logCommand(this)
        myLogger.info("Convert VCF to PS4G")

        val command = headerCommand(this)


        val header = listOf<String>() //TODO update this eventually to take the needed parts out of the header

        val positionSampleGameteLookup = createPositionSampleGameteLookup(gameteVcf)

        val toImputeMap = createPS4GData(sampleVcf, positionSampleGameteLookup)

        toImputeMap.forEach { (sampleGamete, triple) ->
            val (ps4GData, sampleGameteCount, gameteToIdxMap) = triple
            val outputFile = PS4GUtils.buildOutputFileName(sampleVcf,outputDir)

            PS4GUtils.writeOutPS4GFile(ps4GData, sampleGameteCount, gameteToIdxMap, outputFile, header, command)
        }
//        val (ps4GData, sampleGameteCount, gameteToIdxMap) = createPS4GData(sampleVcf, positionSampleGameteLookup)


//        PS4GUtils.writeOutPS4GFile(pS4GData, sampleGameteCount, gameteToIdxMap, outputFile, header, command)
    }

    fun createPositionSampleGameteLookup(gameteVcf: String): Map<Position, Map<SampleGamete, String>> {
        val positionMap = mutableMapOf<Position, Map<SampleGamete, String>>()

        //TODO use custom reader
        VCFFileReader(File(gameteVcf)).use { reader ->
            reader.forEach { record ->
                val position = Position(record.contig, record.start)
                val gameteMap = mutableMapOf<SampleGamete,String>()
                val genotypes = record.genotypes
                genotypes.forEach { genotype ->

                    genotype.alleles.forEachIndexed { index, allele ->
                        val sampleGamete = SampleGamete(genotype.sampleName, index)
                        gameteMap[sampleGamete] = allele.baseString
                    }

                }
                positionMap[position] = gameteMap
            }
        }

        return positionMap
    }

    fun createPS4GData(sampleVcf: String, positionSampleGameteLookup: Map<Position, Map<SampleGamete, String>>): Map<SampleGamete,Triple<List<PS4GData>, Map<SampleGamete,Int>,Map<SampleGamete,Int>>> {

        //(ps4GData, sampleGameteCount, gameteToIdxMap)
        val toImputeMap = mutableMapOf<SampleGamete,Triple<List<PS4GData>, Map<SampleGamete,Int>,Map<SampleGamete,Int>>>()

        //Make a global map of gametes to index
        val gameteToIdxMap = positionSampleGameteLookup.values.first().keys.mapIndexed { index, sampleGamete -> sampleGamete to index }.toMap()



        val sampleGameteCount = positionSampleGameteLookup.values.first().size

        val ps4GData = mutableListOf<List<Int>>()

        VCFFileReader(File(sampleVcf)).use { reader ->
            reader.forEach { record ->
                val position = Position(record.contig, record.start)
                val gameteMap = positionSampleGameteLookup[position] ?: throw IllegalStateException("Position $position not found in gamete VCF")
                val data = mutableListOf<Int>()
                gameteMap.forEach { (sampleGamete, allele) ->
                    val idx = gameteToIdxMap[sampleGamete] ?: throw IllegalStateException("SampleGamete $sampleGamete not found in gameteToIdxMap")
                    data.add(idx)
                }
                ps4GData.add(data)
            }
        }
        TODO()
//        return Triple(ps4GData, sampleGameteCount, gameteToIdxMap)
    }
}