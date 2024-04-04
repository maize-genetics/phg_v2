package net.maizegenetics.phgv2.utils

import net.maizegenetics.phgv2.cli.AlignAssemblies
import org.jetbrains.kotlinx.dataframe.DataFrame
import org.jetbrains.kotlinx.dataframe.api.print
import org.jetbrains.kotlinx.dataframe.api.toMap
import org.jetbrains.kotlinx.dataframe.io.readDelim
import org.jetbrains.kotlinx.dataframe.io.writeCSV
import org.jetbrains.letsPlot.export.ggsave
import org.jetbrains.letsPlot.facet.facetGrid
import org.jetbrains.letsPlot.geom.geomDensity
import org.jetbrains.letsPlot.geom.geomPoint
import org.jetbrains.letsPlot.ggsize
import org.jetbrains.letsPlot.label.labs
import org.jetbrains.letsPlot.letsPlot
import org.jetbrains.letsPlot.scale.scaleXContinuous
import org.jetbrains.letsPlot.scale.scaleYContinuous
import org.junit.jupiter.api.Test
import java.io.File

class MetricUtilsTest {

    @Test
    fun basicLetsPlotTest() {
        // This copied from this jupyter notebook: https://github.com/JetBrains/lets-plot-kotlin/blob/master/docs/examples/jupyter-notebooks/export_to_file.ipynb

        // this example works, so the problem is with the data I'm giving it.
        val rand = java.util.Random(123)
        val n = 400
        val data = mapOf (
            "rating" to List(n/2) { rand.nextGaussian() } + List(n/2) { rand.nextGaussian() * 1.5 + 1.5 },
            "cond" to List(n/2) { "A" } + List(n/2) { "B" }
        )

        var plot1 = letsPlot(data) +
                geomDensity { x = "rating"; color = "cond" } + ggsize(500, 250)
        val pathSVG = ggsave(plot1, "/Users/lcj34/notes_files/phg_v2/newFeatures/stats/letsPlotResults/density.svg")


    }

    @Test
    fun test2BasicAnchorwave() {
        //val p = letsPlot(data) {x="x"; y="y"} + ggsize(600,300) +
        //        geomPoint(color="black", alpha=.1)
        // From https://github.com/JetBrains/lets-plot-kotlin/blob/master/docs/examples/jupyter-notebooks/density_2d.ipynb

        val origFile = File("/Users/lcj34/notes_files/phg_v2/newFeatures/stats/dummy_anchors_smallCSV.anchorspro")
        //val origFile = File("/Users/lcj34/notes_files/phg_v2/newFeatures/stats/CG44_assembly.bp.p_ctg_Zm-B73-REFERENCE-NAM-5.0.anchorspro")
        // Filter out lines that start with '#' and join the rest with newline characters
        val cleanContent = origFile.useLines { lines ->
            lines.filterNot { it.startsWith("#") }
                .map { it.replace("\t", ",") }
                .joinToString("\n")
        }

        val dfAnchorWave = DataFrame.readDelim(cleanContent.reader())
        val data = dfAnchorWave.toMap()
        // THis works.  LetsPlot does not like the color option based on strand in geomPoint,
        // but takes it in the initial letsPlot call.
        // or scrore that is used in the R version.
        val toMb = { x: Number -> x.toDouble() / 1e6 }
        //val toMbLabel: (Double) -> String = { value -> "${toMb(value)}" }
        val toMbLabel:  (Double) -> String = { value -> "${toMb(value)}" }
        val plot2 = letsPlot(data) {x="queryStart"; y="referenceStart";color = "strand" } +
                geomPoint(size = 1.5) +
                scaleXContinuous(labels = listOf(toMbLabel.toString())) +
                scaleYContinuous(labels = listOf(toMbLabel.toString())) +
                //facetGrid(x="queryChr", scales = "fixed")
                facetGrid(x="queryChr", y="refChr", scales = "fixed") +
                labs(x = "Query (Mbp)", y = "Reference (Mbp)")

        //val pathSVG = ggsave(plot2, "/Users/lcj34/notes_files/phg_v2/newFeatures/stats/letsPlotResults/test2BasicAnchorwaveCG44.png")
        val pathSVG = ggsave(plot2, "/Users/lcj34/notes_files/phg_v2/newFeatures/stats/letsPlotResults/test2BasicAnchorwaveFacetXplusY.png")


    }
    @Test
    fun basicLetsPlotForAnchorwave() {
        val origFile = File("/Users/lcj34/notes_files/phg_v2/newFeatures/stats/dummy_anchors_smallCSV.anchorspro")
        // Filter out lines that start with '#' and join the rest with newline characters
        val cleanContent = origFile.useLines { lines ->
            lines.filterNot { it.startsWith("#") }
                .joinToString("\n")
        }

        val dfAnchorWave = DataFrame.readDelim(cleanContent.reader())
        val data = dfAnchorWave.toMap()
        var p = letsPlot(data) +
                geomDensity { x = "referenceStart"; color = "refChr" } + ggsize(500, 250)
        val pathSVG = ggsave(p, "/Users/lcj34/notes_files/phg_v2/newFeatures/stats/letsPlotResults/testBasicAnchorwave.svg")

    }

    @Test
    fun testConvertTxtToCsv() {
        val origFile = File("/Users/lcj34/notes_files/phg_v2/newFeatures/stats/dummy_anchors_small.anchorspro")
        // give the file name origFile, filter out lines that start with '#', for the remaining lines,
        // change all tab characters to commas, and then join the lines with newline characters
        val cleanContent = origFile.useLines { lines ->
            lines.filterNot { it.startsWith("#") }
                .map { it.replace("\t", ",") }
                .joinToString("\n")
        }
        val outputFile = "/Users/lcj34/notes_files/phg_v2/newFeatures/stats/dummy_anchors_small.anchorspro.csv"
        // write the output file
        File(outputFile).writeText(cleanContent)
    }

    @Test
    fun testDataFrame() {
       // val origFile = File("/Users/lcj34/notes_files/phg_v2/newFeatures/stats/dummy_anchors_small.anchorspro.txt")
        //val origFile = File("/Users/lcj34/git/rPHG/inst/extdata/dummy_anchors_small.anchorspro")
        val origFile = File("/Users/lcj34/notes_files/phg_v2/newFeatures/stats/dummy_anchors_smallCSV.anchorspro")
        // Filter out lines that start with '#' and join the rest with newline characters
        val cleanContent = origFile.useLines { lines ->
            lines.filterNot { it.startsWith("#") }
                .map { it.replace("\t", ",") }
                .joinToString("\n")
        }

        val dfAnchorWave = DataFrame.readDelim(cleanContent.reader())
        dfAnchorWave.print()
        val outputFile = "/Users/lcj34/notes_files/phg_v2/newFeatures/stats/testFile_kotlinDF_OutCSV.txt"
        dfAnchorWave.writeCSV(outputFile)

        println("DataFrame written to $outputFile")
        println("Try to plot it with plotDot!")
        // Now try letsplot!

        val plot = AlignAssemblies().plotDot(dfAnchorWave)
        //println("plotDot was successful!  try to show it")
        //plot.show()

        println("plot was created, try to save to a file")
        //val pathSVG = ggsave(plot, "alignmentDotPlot.png",path="/Users/lcj34/notes_files/phg_v2/newFeatures/stats")
        val pathSVG = ggsave(plot, "dotPlotFromAlignAssemblies.png",path="/Users/lcj34/notes_files/phg_v2/newFeatures/stats/letsPlotResults")
        //HTML(File(pathSVG).readText()) // from an example notebook, doesn' work here.

    }
}