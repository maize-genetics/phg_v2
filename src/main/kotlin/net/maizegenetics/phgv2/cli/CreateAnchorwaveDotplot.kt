package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import org.apache.logging.log4j.LogManager
import org.jetbrains.kotlinx.dataframe.DataFrame
import org.jetbrains.kotlinx.dataframe.io.readDelim
import org.jetbrains.letsPlot.export.ggsave
import java.io.File

/**
 * This class provides a means of testing  the dotplot creation functionality of the AlignAssemblies class.
 * It provides for creating the dotplot independently of the AlignAssemblies class.  This is useful
 * when users have an anchorspro file created from anchorwave outside of phgv2 pipeline.
 *
 * It is recommended the output file ends with .svg.  While .png will work, for some assemblies
 * png does not render well.  This has been the case when the assembly is not at a chromosomal level.
 * SVG files have been more reliable in these cases.
 */
class CreateAnchorwaveDotplot: CliktCommand(help = "create a dot plot stored in PNG file from anchowave anchorspro file") {

    private val myLogger = LogManager.getLogger(CreateAnchorwaveDotplot::class.java)

    val inputFile by option(help = "Full path to the anchorwave generated anchorspro file that will be analyzed.")
        .required()
    val outputFile by option(help = "Full path to the SVG file where the dotplot data will be stored - should end with .svg or .png")
        .required()
    override fun run() {

        logCommand(this)

        val startTime = System.nanoTime()

        val origFile = File(inputFile)
        // Filter out lines that start with '#', change tabs to comma, and join the rest with newline characters
        // we change tabs to commas as the Kotlin DataFrame reader appears to be expecting CSV format,
        // tabs were not working as a delimiter.
        val cleanContent = origFile.useLines { lines ->
            lines.filterNot { it.startsWith("#") }
                .map { it.replace("\t", ",") }
                .joinToString("\n")
        }

        val dfAnchorWave = DataFrame.readDelim(cleanContent.reader())

        // Three, two, one, plot!
        val plot = AlignAssemblies().plotDot(dfAnchorWave)


        // We specified in the help message that the output file shoudl be a full path.
        // But since we are handling relative paths in AlignAssemblies, we will check for that here.
        val plotFile = if (outputFile.startsWith("/")) {
            outputFile
        } else {
            "${System.getProperty("user.dir")}/${outputFile}"
        }
        myLogger.info("plot was created, call ggsave")
        val pathSVG = ggsave(plot, plotFile)

        myLogger.info("dotplot output file written to $pathSVG ")
        val totalTime = (System.nanoTime() - startTime)/1e9
        myLogger.info("Time to process stats: $totalTime seconds")

    }
}
