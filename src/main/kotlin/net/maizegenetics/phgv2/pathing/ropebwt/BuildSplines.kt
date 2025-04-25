package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.choice
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.cli.logCommand
import org.apache.logging.log4j.LogManager

class BuildSplines: CliktCommand(help = "Build Splines from gVCFs or hVCFs")  {

    val myLogger = LogManager.getLogger(BuildSplines::class.java)

    val vcfDir by option(help = "Directory containing the hvcf or gvcf files")
        .required()

    val vcfType by option(help = "Type of vcfs to build the splines")
        .choice("hvcf","gvcf")
        .default("hvcf")

    val outputFile by option(help = "Output file")
        .required()

    val minIndelLength by option(help="Minimum length of an indel to break up the running block for spline creation of gvcfs.  If --vcf-type is hvcf this option is ignored.")
        .int()
        .default(10)

    val maxNumPointsPerChrom by option(help = "Number of points per chrom.  If there are more points for each sample's chromosomes we will downsample randomly..")
        .int()
        .default(250_000)

    override fun run() {
        logCommand(this)

        val (splineLookup, chrIndexMap, gameteIndexMap) = SplineUtils.buildSplineLookup(vcfDir, vcfType, minIndelLength, maxNumPointsPerChrom)

        myLogger.info("Writing out splines to $outputFile")

        SplineUtils.writeSplinesToFile(splineLookup, chrIndexMap, gameteIndexMap, outputFile)
    }

}