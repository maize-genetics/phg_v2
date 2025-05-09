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
import kotlin.collections.mutableMapOf

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
    }


    fun createPositionSampleGameteLookup(gameteVcf: String): Map<Position, Map<String, List<SampleGamete>>> {
       val positionMap = mutableMapOf<Position, Map<String, List<SampleGamete>>>()

        //TODO use custom reader
        VCFFileReader(File(gameteVcf),false).use { reader ->
            reader.map { record ->
                val position = Position(record.contig, record.start)
                val genotypes = record.genotypes
                val alleleMap = genotypes.flatMap { genotype ->
                    genotype.alleles.mapIndexed { index, allele ->
                        val sampleGamete = SampleGamete(genotype.sampleName, index)
                        Pair(allele.baseString, sampleGamete)
                    }
                }.groupBy({ it.first }, { it.second })

                positionMap[position] = alleleMap
            }
        }
        return positionMap
    }



    fun createPS4GData(sampleVcf: String, positionSampleGameteLookup: Map<Position, Map<String, List<SampleGamete>>>): Map<SampleGamete,Triple<List<PS4GData>, Map<SampleGamete,Int>,Map<SampleGamete,Int>>> {

        val gameteToCountMap = mutableMapOf<SampleGamete,Pair<MutableMap<Pair<Int, List<Int>>, Int>, MutableMap<SampleGamete,Int>>>()

        //Make a global map of gametes to index
        val gameteToIdxMap = createGameteToIdxMap(positionSampleGameteLookup)

        val sampleNameToIdxMap = gameteToIdxMap.keys
            .mapIndexed { index, sampleGamete ->
                sampleGamete.name to index
            }.toMap()


        //This map is Map<toImputeSampleGamete, Map<refPanelSampleGamete, totalCountOfPositionsThatMatch>>
        val sampleGameteCount = mutableMapOf<SampleGamete, MutableMap<SampleGamete,Int>>() //Need to build this on the fly
        var missedPositionCount = 0

        VCFFileReader(File(sampleVcf)).use { reader ->
            reader.forEach { record ->
                val position = Position(record.contig, record.start)
                if(!positionSampleGameteLookup.containsKey(position)) {
                    missedPositionCount++
                    return@forEach
                }

                val encodedPosition = PS4GUtils.encodePosition(position, sampleNameToIdxMap)

                //Get out the alleles and lists of SampleGametes for this position
                val alleleMap = positionSampleGameteLookup[position]!!

                //Need to loop through each genotype
                val genotypes = record.genotypes

                genotypes.forEach { genotype ->
                    genotype.alleles.forEachIndexed { index, allele ->
                        val sampleGamete = SampleGamete(genotype.sampleName, index)
                        val alleleString = allele.baseString

                        //Get the list of SampleGametes for this allele
                        val sampleGametes = alleleMap[alleleString] ?: emptyList()

                        //Add to the map of counts
                        sampleGameteCount.getOrPut(sampleGamete) { mutableMapOf() }.let { countMap ->
                            sampleGametes.forEach { refSampleGamete ->
                                countMap[refSampleGamete] = countMap.getOrDefault(refSampleGamete, 0) + 1
                            }
                        }
                        //Need to add this set of SampleGametes to the specific PS4GData
                        val (countMap, sampleGameteCountMap) = gameteToCountMap[sampleGamete] ?: Pair(mutableMapOf(), mutableMapOf())

                        //increment the count for these hits from SampleGametes
                        //First convert SampleGametes to indices
                        val sampleGameteIndices = sampleGametes.map { gameteToIdxMap[it] ?: throw IllegalStateException("SampleGamete $it not found in gameteToIdxMap") }
                        //Get the count for this position and sampleGamete
                        countMap[Pair(encodedPosition, sampleGameteIndices)] = countMap.getOrDefault(Pair(encodedPosition, sampleGameteIndices), 0) + 1

                        //update the count for these SampleGametes
                        for(sampleGamete in sampleGametes) {
                            sampleGameteCountMap[sampleGamete] = sampleGameteCountMap.getOrDefault(sampleGamete, 0) + 1
                        }

                        gameteToCountMap[sampleGamete] = Pair(countMap, sampleGameteCountMap)
                    }
                }
            }
        }
        //Go through the toImputeMap and create the PS4GData
        //First thing we need to do is convert the count map to a list of PS4G data

        return convertCountMapsToData(gameteToCountMap, gameteToIdxMap)
    }

    fun createGameteToIdxMap(positionSampleGameteLookup: Map<Position, Map<String, List<SampleGamete>>>): Map<SampleGamete, Int> {
        return positionSampleGameteLookup.values.flatMap { it.values }
            .flatten()
            .distinct()
            .sorted()
            .mapIndexed { index, sampleGamete ->
                sampleGamete to index
            }.toMap()
    }

    fun convertCountMapsToData(countMaps: Map<SampleGamete,Pair<Map<Pair<Int, List<Int>>, Int>, Map<SampleGamete,Int>>>, gameteToIdxMap: Map<SampleGamete, Int>) : Map<SampleGamete,Triple<List<PS4GData>, Map<SampleGamete,Int>,Map<SampleGamete,Int>>> {
        return countMaps.map { (sampleGamete, pair) ->
            val (countMap, sampleGameteCountMap) = pair
            val ps4GDataList = PS4GUtils.convertCountMapToPS4GData(countMap)
            Pair(sampleGamete, Triple(ps4GDataList, sampleGameteCountMap, gameteToIdxMap))
        }.toMap()
    }
}