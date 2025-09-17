package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.pathing.AlignmentUtils
import org.apache.logging.log4j.LogManager
import java.io.File

class MapReadsFromBed : CliktCommand(help = "Map reads to a ropebwt3 index from a bed file.") {
    private val myLogger = LogManager.getLogger(MapReads::class.java)

    val bedFile by option(help = "The full path of the bed file.").required()

    val outputDir by option("-o", "--output-dir", help = "Name for output ReadMapping file Directory (Required)")
        .required()

    val hvcfDir by option(help = "Directory for the haplotype VCF files. These hvcfs will be used to filter the hits to only those that hit a single referenceRange.")
        .required()

    val threads by option(help = "Number of threads to use.")
        .int()
        .default(5)

    val maxNumHits by option(help = "Number of hits to report.  Note ropebwt can hit more than --max-num-hits but any alignment hitting more haplotypes than this will be ignored.")
        .int()
        .default(50)

    val condaEnvPrefix by option (help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
        .default("")

    val maxStart by option(help = "Maximum start position for a read to be considered a match. Any alignments with a start above this will be ignored.")
        .int()
        .default(0)

    val minEnd by option(help = "Minimum end position for a read to be considered a match. Any alignments with an end below this will be ignored.")
        .int()
        .default(70)

    val sampleName by option(help = "Sample name to use in the output read mapping file.")
        .required()

    override fun run() {
        logCommand(this)

        myLogger.info("Building the Graph")

        val graph = HaplotypeGraph(hvcfDir)

        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()

        myLogger.info("Mapping reads to pangenome")

        val mapReads = MapReads()
        val bedFileReader = bufferedReader(bedFile)
        val readMapping = mapReads.createReadMappingsForFileReader(bedFileReader, maxNumHits, hapIdToRefRangeMap, maxStart, minEnd)

        bedFileReader.close()

        val outputFile = "$outputDir/${File(bedFile).nameWithoutExtension}_readMapping.txt"

        myLogger.info("Writing read mapping to $outputFile")
        AlignmentUtils.exportReadMapping(outputFile,readMapping, sampleName,Pair(bedFile,""))


    }
}