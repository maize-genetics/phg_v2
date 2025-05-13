package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.VariantContext
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

data class PS4GDataWithMaps(val ps4gData: List<PS4GData>, val sampleGameteCount: Map<SampleGamete,Int>, val gameteToIdxMap: Map<SampleGamete,Int>)

/**
 * This data class holds the two SampleGameteCountMaps that are used to create the PS4GData
 * The first map is the overall positional count.  The key is Pair<encodedPosition, listOfReferenceSampleGameteIndices>
 * The second is the reference sample gamete count map.  The key is the SampleGamete and the value is the count of how many times this SampleGamete was seen in the reference panel.
 */
data class SampleGameteCountMaps(val countMap: MutableMap<Pair<Int, List<Int>>, Int>, val refSampleGameteCountMap: MutableMap<SampleGamete, Int>)


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



    fun createPS4GData(sampleVcf: String,
                       positionSampleGameteLookup: Map<Position, Map<String, List<SampleGamete>>>): Map<SampleGamete, PS4GDataWithMaps> {

        val gameteToCountMap = mutableMapOf<SampleGamete, SampleGameteCountMaps>()

        //Make a global map of gametes to index
        val gameteToIdxMap = createGameteToIdxMap(positionSampleGameteLookup)

        val contigNameToIdxMap = positionSampleGameteLookup.keys.map { position -> position.contig }.distinct()
            .mapIndexed { index, contig ->
                contig to index
            }.toMap()


        //This map is Map<toImputeSampleGamete, Map<refPanelSampleGamete, totalCountOfPositionsThatMatch>>
        val sampleGameteCount = mutableMapOf<SampleGamete, MutableMap<SampleGamete,Int>>() //Need to build this on the fly
        var missedPositionCount = 0

        VCFFileReader(File(sampleVcf), false).use { reader ->
            reader.forEach { record ->
                val position = Position(record.contig, record.start)
                if(!positionSampleGameteLookup.containsKey(position)) {
                    missedPositionCount++
                    return@forEach
                }

                processVariantPosition(
                    position,
                    contigNameToIdxMap,
                    positionSampleGameteLookup,
                    record,
                    sampleGameteCount,
                    gameteToCountMap,
                    gameteToIdxMap
                )
            }
        }

        myLogger.info("Number of positions not found in the reference panel: $missedPositionCount")

        //Go through the toImputeMap and create the PS4GData
        //First thing we need to do is convert the count map to a list of PS4G data

        return convertCountMapsToData(gameteToCountMap, gameteToIdxMap)
    }

    fun processVariantPosition(
        position: Position,
        contigNameToIdxMap: Map<String, Int>,
        positionSampleGameteLookup: Map<Position, Map<String, List<SampleGamete>>>,
        record: VariantContext,
        sampleGameteCount: MutableMap<SampleGamete, MutableMap<SampleGamete, Int>>,
        gameteToCountMap: MutableMap<SampleGamete, SampleGameteCountMaps>,
        gameteToIdxMap: Map<SampleGamete, Int>
    ) {
        val encodedPosition = PS4GUtils.encodePosition(position, contigNameToIdxMap)

        //Get out the alleles and lists of SampleGametes for this position
        val alleleMap = positionSampleGameteLookup[position]!!

        //Need to loop through each genotype
        val genotypes = record.genotypes

        genotypes.forEach { genotype ->
            genotype.alleles.forEachIndexed { alleleIndex, allele ->
                processSingleGenotype(
                    genotype,
                    alleleIndex,
                    allele,
                    alleleMap,
                    sampleGameteCount,
                    gameteToCountMap,
                    gameteToIdxMap,
                    encodedPosition
                )
            }
        }
    }

    /**
     * Function to process a single genotype and add it to the count map
     * We update teh sampleGameteCount to have a running count of each of the reference gametes that were hit by this SNP
     * We also update the gameteToIdxMap as we may have seen a new gamete this time.
     */
    fun processSingleGenotype(
        genotype: Genotype,
        alleleIndex: Int,
        allele: Allele,
        alleleMap: Map<String, List<SampleGamete>>, //This is coming from the refPanel
        sampleGameteCount: MutableMap<SampleGamete, MutableMap<SampleGamete, Int>>,
        gameteToCountMap: MutableMap<SampleGamete, SampleGameteCountMaps>,
        gameteToIdxMap: Map<SampleGamete, Int>,
        encodedPosition: Int
    ) {
        val sampleGamete = SampleGamete(genotype.sampleName, alleleIndex)
        val alleleString = allele.baseString

        //Get the list of RefPanel SampleGametes for this allele
        val refPanelSampleGametes = alleleMap[alleleString] ?: emptyList()

        //Add to the map of counts
        sampleGameteCount.getOrPut(sampleGamete) { mutableMapOf() }.let { countMap ->
            refPanelSampleGametes.forEach { refSampleGamete ->
                countMap[refSampleGamete] = countMap.getOrDefault(refSampleGamete, 0) + 1
            }
        }
        //Need to add this set of SampleGametes to the specific PS4GData
        val (countMap, refSampleGameteCountMap) = gameteToCountMap[sampleGamete] ?: SampleGameteCountMaps(
            mutableMapOf(),
            mutableMapOf()
        )

        //increment the count for these hits from RefPanel SampleGametes
        //First convert SampleGametes to indices
        val sampleGameteIndices = refPanelSampleGametes.map {
            gameteToIdxMap[it] ?: throw IllegalStateException("SampleGamete $it not found in gameteToIdxMap")
        }
        //Get the count for this position and sampleGamete
        countMap[Pair(encodedPosition, sampleGameteIndices)] =
            countMap.getOrDefault(Pair(encodedPosition, sampleGameteIndices), 0) + 1

        //update the count for these SampleGametes
        for (refSampleGamete in refPanelSampleGametes) {
            refSampleGameteCountMap[refSampleGamete] = refSampleGameteCountMap.getOrDefault(refSampleGamete, 0) + 1
        }

        gameteToCountMap[sampleGamete] = SampleGameteCountMaps(countMap, refSampleGameteCountMap)
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

    fun convertCountMapsToData(countMaps: Map<SampleGamete, SampleGameteCountMaps>, gameteToIdxMap: Map<SampleGamete, Int>) : Map<SampleGamete, PS4GDataWithMaps> {
        return countMaps.map { (sampleGamete, pair) ->
            val (countMap, sampleGameteCountMap) = pair
            val ps4GDataList = PS4GUtils.convertCountMapToPS4GData(countMap)
            Pair(sampleGamete, PS4GDataWithMaps(ps4GDataList, sampleGameteCountMap, gameteToIdxMap))
        }.toMap()
    }
}