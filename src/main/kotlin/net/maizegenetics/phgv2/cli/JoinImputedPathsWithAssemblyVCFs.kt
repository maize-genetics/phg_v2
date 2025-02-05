package net.maizegenetics.phgv2.cli

import biokotlin.genome.Position
import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.boolean
import com.github.ajalt.clikt.parameters.types.int
import com.google.common.collect.ArrayListMultimap
import com.google.common.collect.Multimap
import net.maizegenetics.analysis.association.FixedEffectLMPlugin
import net.maizegenetics.analysis.data.IntersectionAlignmentPlugin
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin
import net.maizegenetics.dna.snp.*
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder
import net.maizegenetics.phenotype.Phenotype
import net.maizegenetics.phenotype.PhenotypeBuilder
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.Datum
import net.maizegenetics.taxa.TaxaListBuilder
import net.maizegenetics.util.TableReport
import net.maizegenetics.util.TableReportUtils
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.*

class JoinImputedPathsWithAssemblyVCFs :
    CliktCommand(help = "Joins imputed paths with assembly VCFs, and runs GLM for each reference range.") {

    private val myLogger = LogManager.getLogger(JoinImputedPathsWithAssemblyVCFs::class.java)

    val vcfFilesPerRangeForAssembliesDir by option(help = "Directory containing VCF files per range for assemblies")
        .required()

    val traitFile by option(help = "Phenotype file")
        .required()

    // Set to null if no population structure file
    val populationStructureFile by option(help = "Population structure file")
        .default("")

    // CHROM   POS     ID      REF     ALT  Samples...
    // /local/workdir/wl748/merge_pangenome_v2.txt
    val pangenomeTableFile by option(help = "Pangenome table file")
        .required()

    // #CHROM  POS Samples...
    // /local/workdir/wl748/merge_hapID/merge_SeeD.txt
    val imputedTableFile by option(help = "Imputed table file")
        .required()

    val indelsToMissing by option(help = "Convert indels to missing genotypes")
        .boolean()
        .default(false)

    val siteMinCount: Int? by option(help = "Site minimum count filter")
        .int()

    override fun run() {

        // Create a map of key (i.e. chr10_154585427) to VCF file name (i.e. Zh-chr10_154585427-154627028.vcf)
        val vcfFileNamesPerRangeForAssemblies = File(vcfFilesPerRangeForAssembliesDir)
            .walk()
            .map { it.absolutePath }
            .filter { it.endsWith(".vcf") || it.endsWith(".vcf.gz") }
            .map { it.substringBeforeLast("-").substringAfter("Zh-") to it }
            .toMap()

        val phenotype = PhenotypeBuilder().fromFile(traitFile).build()[0]

        val populationStructure = if (populationStructureFile.isNotEmpty()) {
            PhenotypeBuilder().fromFile(populationStructureFile).build()[0]
        } else {
            null
        }

        val pangenomeTable = readTable(pangenomeTableFile, 4)

        val imputedTable = readTable(imputedTableFile, 2)

        imputedTable.posToLine.forEach { (pos, line) ->

            val key = "${pos.contig}_${pos.position}"
            val vcfFilename = vcfFileNamesPerRangeForAssemblies[key]
            if (vcfFilename == null) {
                println("WARNING: No VCF file found for key: $key")
                return@forEach
            }

            val genotypeTable = processRange(pos, line, vcfFilename, indelsToMissing)

            val filteredGenotypeTable = filterGenotypeTable(genotypeTable, siteMinCount, siteMinAlleleFreq, bedfile)

            if (filteredGenotypeTable != null && filteredGenotypeTable.numberOfSites() > 0) {

                if (writeVCFFiles) writeVCF(
                    filteredGenotypeTable,
                    "impute-by-range/${vcfFilename.substringAfterLast('/').replace("Zh", "Impute")}"
                )

                val glmOutput = runGLM(filteredGenotypeTable, phenotype, populationStructure)

                if (writeGLMResults) {
                    glmOutput.forEachIndexed { i, table ->
                        writeTable(
                            table,
                            "glm-by-range/${
                                vcfFilename.substringAfterLast('/').replace("Zh", "GLM").replace(".vcf", "-$i.txt")
                            }"
                        )
                    }
                }

            }

        }

    }

    // include Site Min Count, Site Min Allele Freq and Bed file filters
    private fun filterGenotypeTable(
        genotype: GenotypeTable,
        siteMinCount: Int? = null,
        siteMinAlleleFreq: Double? = null,
        bedfile: String? = null
    ): GenotypeTable? {

        // If no filters are set, return the original genotype table
        if (siteMinCount == null && siteMinAlleleFreq == null && bedfile == null) return genotype

        val builder = FilterSiteBuilderPlugin()

        siteMinCount?.let { builder.siteMinCount(it) }
        siteMinAlleleFreq?.let { builder.siteMinAlleleFreq(it) }
        bedfile?.let { builder.bedFile(it) }

        return builder.runPlugin(genotype)

    }

    fun runGLM(
        genotype: GenotypeTable,
        phenotype: Phenotype,
        populationStructure: Phenotype? = null
    ): List<TableReport> {

        val inputDatums = mutableListOf<Datum>()
        inputDatums.add(Datum("genotype", genotype, null))
        inputDatums.add(Datum("phenotype", phenotype, null))
        populationStructure?.let { inputDatums.add(Datum("populationStructure", it, null)) }

        val input = DataSet(inputDatums, null)

        val intersect = IntersectionAlignmentPlugin(null, false)

        val glm = FixedEffectLMPlugin(null, false)

        val result = glm.performFunction(intersect.performFunction(input))

        return result.getDataOfType(TableReport::class.java).map { it.data as TableReport }

    }

    /**
     * @param samples List of samples in the table
     * @param posToLine Map of positions to lines in the file
     */
    data class SummaryTable(val samples: List<String>, val posToLine: SortedMap<Position, String>)

    private fun readTable(filename: String, nonSampleHeaders: Int): SummaryTable {
        bufferedReader(filename).use { reader ->
            val samples = reader.readLine().split("\t").drop(nonSampleHeaders)
            val posToLine = reader.lineSequence()
                .map {
                    val tokens = it.split("\t", limit = 3)
                    Position(tokens[0], tokens[1].toInt()) to it
                }.toMap().toSortedMap()
            return SummaryTable(samples, posToLine)
        }
    }

    private fun writeVCF(genotype: GenotypeTable, filename: String) {
        ExportUtils.writeToVCF(genotype, filename, false)
    }

    private fun writeTable(table: TableReport, filename: String) {
        TableReportUtils.saveDelimitedTableReport(table, "\t", File(filename))
    }

}