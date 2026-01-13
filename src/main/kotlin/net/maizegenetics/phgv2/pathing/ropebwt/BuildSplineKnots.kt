package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.choice
import com.github.ajalt.clikt.parameters.types.int
import com.github.ajalt.clikt.parameters.types.long
import net.maizegenetics.phgv2.cli.logCommand
import org.apache.logging.log4j.LogManager
import java.io.File

class BuildSplineKnots: CliktCommand(help = "Build Spline Knot points from gVCFs or hVCFs")  {

    val myLogger = LogManager.getLogger(BuildSplineKnots::class.java)

    val vcfDir by option(help = "Directory containing the hvcf or gvcf files")
        .required()

    val vcfType by option(help = "Type of vcfs to build the splines")
        .choice("hvcf","gvcf")
        .default("hvcf")

    val outputDir by option(help = "Output Directory to write the spline knots to.")
        .required()

    val minIndelLength by option(help="Minimum length of an indel to break up the running block for spline creation of gvcfs.  If --vcf-type is hvcf this option is ignored.")
        .int()
        .default(10)

    val numBpsPerKnot by option(help = "Number of bps per knot for each contig's spline.  This is the maximum number of bps per knot, so the actual number may be less if there are not enough bps in the contig.")
        .int()
        .default(50_000)

    val disableSplineDownsampling by option(help = "Keep all spline knots regardless of --num-bps-per-knot")
        .flag()

    val contigList by option(help = "List of chromosomes to include in the splines.  If not provided, all chromosomes will be included.")
        .default("")

    val randomSeed by option(help = "Random seed for downsampling the number of points per chromosome.  If not provided, the seed 12345 will be used.")
        .long()
        .default(12345)

    override fun run() {
        logCommand(this)

        //Check to see if outputDir exists, if not create it
        if( !File(outputDir).exists() ) {
            myLogger.info("Output directory $outputDir does not exist, creating it.")
            File(outputDir).mkdirs()
        }


        myLogger.info("Building Spline Knots from $vcfType files in $vcfDir and writing Knot Files to $outputDir")
        val contigSet = contigList.split(",").map { it.trim() }.filter { it.isNotEmpty() }.toSet()


        //This now will write out the spline knots for each assembly to the output directory automatically.
        SplineUtils.buildSplineKnots(
            vcfDir,
            vcfType,
            outputDir,
            minIndelLength,
            numBpsPerKnot,
            contigSet,
            disableSplineDownsampling,
            randomSeed
        )

        myLogger.info("Spline Knot building complete.")
    }

}