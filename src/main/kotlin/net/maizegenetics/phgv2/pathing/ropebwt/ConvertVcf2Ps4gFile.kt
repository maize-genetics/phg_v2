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

        val outputFile = PS4GUtils.buildOutputFileName(sampleVcf,outputDir)

        val header = listOf<String>() //TODO update this eventually to take the needed parts out of the header

        val positionSampleGameteLookup = createPositionSampleGameteLookup(gameteVcf)

//        val (ps4GData, sampleGameteCount, gameteToIdxMap) = convert


//        PS4GUtils.writeOutPS4GFile(pS4GData, sampleGameteCount, gameteToIdxMap, outputFile, header, command)
    }

    fun createPositionSampleGameteLookup(gameteVcf: String): Map<Position, Map<SampleGamete, String>> {
        val positionMap = mutableMapOf<Position, Map<SampleGamete, String>>()

        //TODO use custom reader
        VCFFileReader(File(gameteVcf)).use { reader ->
            reader.forEach { record ->
                val position = Position(record.contig, record.start)

                val genotypes = record.genotypes
                genotypes.forEach { genotype ->

                    genotype.alleles.mapIndexed { index, allele ->
                        val sampleGamete = SampleGamete(genotype.sampleName, index)

                    }

                }

            }
        }

        return positionMap
    }
}