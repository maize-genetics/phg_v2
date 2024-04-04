package net.maizegenetics.phgv2.utils

import org.jetbrains.kotlinx.dataframe.DataFrame
import org.jetbrains.kotlinx.dataframe.api.filter
import org.jetbrains.kotlinx.dataframe.api.toMap
import org.jetbrains.letsPlot.GGBunch

import org.jetbrains.letsPlot.facet.facetGrid
import org.jetbrains.letsPlot.geom.geomPoint
import org.jetbrains.letsPlot.ggsize
import org.jetbrains.letsPlot.intern.Plot
import org.jetbrains.letsPlot.label.labs
import org.jetbrains.letsPlot.letsPlot
import org.jetbrains.letsPlot.scale.scaleColorManual
import org.jetbrains.letsPlot.scale.scaleColorViridis
import org.jetbrains.letsPlot.scale.scaleXContinuous
import org.jetbrains.letsPlot.scale.scaleYContinuous

//fun plotDot(
//    data: DataFrame<*>,
//    querySeqId: List<String>? = null,
//    refSeqId: List<String>? = null,
//    queryLab: String = "Query",
//    refLab: String = "Reference",
//    colorId: String = "strand"
//): Plot {
////): GGBunch {
//    val filteredData = data.filter {
//        (querySeqId == null || it["queryChr"] in querySeqId) &&
//                (refSeqId == null || it["refChr"] in refSeqId)
//    } // chatGPT filtered based on Brandon's code. Might be ok ....
//
//    val toMb = { x: Number -> x.toDouble() / 1e6 }
//    //val toMbLabel: (Double) -> String = { value -> "${toMb(value)}" }
//    val toMbLabel:  (Double) -> String = { value -> "${toMb(value)}" }
//    val plotData = filteredData.toMap()
//    val plot = letsPlot(plotData ) { // this was from chatGpt based on Brandon's R code.
//    //val p = letsPlot(data.toMap()){
//    //val p = letsPlot(filteredData.toMap()) { // this was from chatGpt based on Brandon's R code.
//        x = "queryStart"; y = "referenceStart"
//        color = when (colorId) {
//            "score" -> "score"
//            else -> "strand" // Default to "strand" if anything else is specified
//        }
//    } +
//            geomPoint(size = 1.5) +
//            //ggsize(800, 350) +
//            // THese were chatGPT suggestions, but they don't compile
//            //scaleXContinuous(labels = { listOf(toMb(it)) }) +
//            //scaleYContinuous(labels = { listOf(toMb(it)) }) +
//            //scaleXContinuous(labels = { it.map { value -> "${toMb(value)} Mbp" } }) +
//            //scaleYContinuous(labels = { it.map { value -> "${toMb(value)} Mbp" } }) +
//            // Below compiles, but at runtime gives the error "The scale 'labels' parameter should be specified with a list or dictionary."
//            //scaleXContinuous(labels = toMbLabel) +
//            //scaleYContinuous(labels = toMbLabel) +
//            scaleXContinuous(labels = listOf(toMbLabel.toString())) +
//            scaleYContinuous(labels = listOf(toMbLabel.toString())) +
//            facetGrid(x="queryChr", y="refChr", scales="fixed") +
//            //facetGrid(facets = "refChr ~ queryChr", scales = "free", space = "free") +
//            //facetGrid(x = "Reference (Mbp)", y= "Query (Mbp)", scales = "free") +
//            labs(x = "$queryLab (Mbp)", y = "$refLab (Mbp)") +
//            when (colorId) {
//                "score" -> scaleColorViridis()
//                "strand" -> scaleColorManual(
//                    values = mapOf("+" to "#DA897C", "-" to "#0D6A82")
//                )
//
//                else -> null
//            }!!
//
//   // return GGBunch().addPlot(p, 0, 0)
//    return plot
//}
