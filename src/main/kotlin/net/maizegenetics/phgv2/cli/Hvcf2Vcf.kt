package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.VCFConversionUtils
import org.apache.logging.log4j.LogManager
import java.util.TreeMap

class Hvcf2Vcf:
    CliktCommand(help = "Create vcf file for a PHG pathing h.vcf using data from existing PHG created vcf files") {

    private val myLogger = LogManager.getLogger(Hvcf2Vcf::class.java)

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        .default("")


    val hvcfDir by option(help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        .required()

    val donorVcfFile by option(help = "Path to the VCF file containing all the PHG SNPs.  This is typically created by running merge-gvcf on the Assembly Gvcf files.")
        .required()

    val outputFile by option(help = "Output file.")
        .required()



    override fun run() {
        processHVCFAndBuildVCF(dbPath, hvcfDir, donorVcfFile, outputFile)
    }

    fun processHVCFAndBuildVCF(dbPath: String, hvcfDir: String, donorVcfFile: String, outputFile: String) {
        //Load in the ASM hapId map so we can make sure we pull out the right SNPs.
        //This is in form (RefRange,ASMName) -> hapID to allow for easy lookup
        val asmHapIdMap = VCFConversionUtils.createASMNameAndRefRangeMap(dbPath)

        //We then need to load in the hvcfFiles and create a map of refRange+SampleName -> HapId1 + HapId2
        //Or maybe it needs to be refRange+HapId -> List<Pair<SampleName,Gamete>>
        val refRangeAndHapIdMap = createRangeHapMapToSampleGamete(hvcfDir)

        val ranges = refRangeAndHapIdMap.map {it.key.first}.toSet()

        //We also need to build a fast Position -> refRange lookup
        val positionToRangeMap = buildPositionToRefRangeMap(ranges)
        //Build an output VCF With the header and sampleNames coming from the outputSampleNames

        //Then we can walk through the mergedVCF file
        //For each position, lookup the RefRange,
        //Then lookup RefRange+ASMName in asmHapIdMap to get out the HapId for each sampleName in the VCF
            //Then for each RefRange+ASMNames hapIds we lookup the SampleName and Gamete and write it to that genotype


    }

    fun createRangeHapMapToSampleGamete(hvcfDir: String): Map<Pair<ReferenceRange, String>, SampleGamete> {
        TODO("CODE THIS UP")
    }

    /**
     * Build a TreeMap for both the start and the end point for each reference range.
     * This should allow us to easily hit the correct reference ranges.
     *
     * TODO Check to see if we can just use the start...We should be able to just makes it a bit trickier logic
     */
    fun buildPositionToRefRangeMap(ranges: Set<ReferenceRange>): TreeMap<Position, ReferenceRange> {
        //Use a TreeMap to allow for fast lookup of the position  A treeRangeMap would also work, but they tend to be slower than simple TreeMap
        //TODO implement primitive map

        val positionMap = TreeMap<Position, ReferenceRange>()
        for (range in ranges) {
            positionMap[Position(range.contig, range.start)] = range
            positionMap[Position(range.contig, range.end)] = range
        }
        return positionMap
    }


}