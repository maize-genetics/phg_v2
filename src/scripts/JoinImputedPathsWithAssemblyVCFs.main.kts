@file:DependsOn("net.maizegenetics:tassel:5.2.94")
@file:DependsOn("org.biokotlin:biokotlin:0.22")

import biokotlin.genome.Position
import biokotlin.util.bufferedReader
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
import net.maizegenetics.util.TableReportBuilder
import net.maizegenetics.util.TableReportUtils
import java.io.File
import java.util.*


//
// Variables to be set by the user
//

val vcfFilesPerRangeForAssembliesDir = "vcf-by-range"

// /workdir/wl748/SeeD_env/allClimateSoilGeoVariables.txt
// Start with trait AltM
val traitFile = "AltM.txt"

// Set to null if no population structure file
val populationStructureFile: String? = "SeeD_gPCs2.txt"

// CHROM   POS     ID      REF     ALT  Samples...
// /local/workdir/wl748/merge_pangenome_v2.txt
val pangenomeTableFile = "merge_pangenome_v2.txt"

// #CHROM  POS Samples...
// /local/workdir/wl748/merge_hapID/merge_SeeD.txt
val imputedTableFile = "merge_SeeD.txt"

val indelsToMissing = false

val writeVCFFiles = true

val writeGLMResults = true

// Filtering options
// Set to null if option not used
val siteMinCount: Int? = null
val siteMinAlleleFreq: Double? = null
val bedfile: String? = null


// Create a map of key (i.e. chr10_154585427) to VCF file name (i.e. Zh-chr10_154585427-154627028.vcf)
val vcfFileNamesPerRangeForAssemblies = File(vcfFilesPerRangeForAssembliesDir)
    .walk()
    .map { it.absolutePath }
    .filter { it.endsWith(".vcf") || it.endsWith(".vcf.gz") }
    .map { it.substringBeforeLast("-").substringAfter("Zh-") to it }
    .toMap()

val phenotype = PhenotypeBuilder().fromFile(traitFile).build()[0]

val populationStructure = populationStructureFile?.let { PhenotypeBuilder().fromFile(it).build()[0] }

val pangenomeTable = readTable(pangenomeTableFile, 4)

val imputedTable = readTable(imputedTableFile, 2)

imputedTable.posToLine.forEach { (pos, line) ->

    val key = "${pos.contig}_${pos.position}"
    val vcfFilename = vcfFileNamesPerRangeForAssemblies[key]
    if (vcfFilename == null) {
        println("WARNING: No VCF file found for key: $key")
        return@forEach
    }

    var genotypeTable = processRange(pos, line, vcfFilename, indelsToMissing)

    genotypeTable = filterGenotypeTable(genotypeTable, siteMinCount, siteMinAlleleFreq, bedfile)

    if (writeVCFFiles) writeVCF(
        genotypeTable,
        "impute-by-range/${vcfFilename.substringAfterLast('/').replace("Zh", "Impute")}"
    )

    val glmOutput = runGLM(genotypeTable, phenotype, populationStructure)

    if (writeGLMResults) {
        glmOutput.forEachIndexed { i, table ->
            writeTable(
                table,
                "glm-by-range/${vcfFilename.substringAfterLast('/').replace("Zh", "GLM").replace(".vcf", "-$i.txt")}"
            )
        }
    }

}


fun processRange(pos: Position, line: String, vcfFilename: String, indelToMissing: Boolean = false): GenotypeTable {

    val pangenomeLine = pangenomeTable.posToLine[pos] ?: error("No pangenome entry found for $pos")
    val pangenomeHapids = pangenomeLine.split("\t").drop(4)
    require(pangenomeHapids.size == pangenomeTable.samples.size) {
        "Number of hapids (${pangenomeHapids.size}) does not match number of samples (${pangenomeTable.samples.size}) at $pos"
    }

    val hapids = line.split("\t").drop(2)
    val numSamples = imputedTable.samples.size
    require(hapids.size == numSamples) {
        "Number of hapids (${hapids.size}) does not match number of samples ($numSamples) at $pos"
    }

    val rangeGenotypeTableAssemblies = ImportUtils.readFromVCF(vcfFilename, null, false, true)
    val numSites = rangeGenotypeTableAssemblies.numberOfSites()
    val vcfSamples = rangeGenotypeTableAssemblies.taxa()

    val pangenomeHapidToSample = pangenomeHapids
        .mapIndexed { sampleIndex, hapid ->
            val sampleName = pangenomeTable.samples[sampleIndex]
            val newSampleIndex = vcfSamples.indexOf(sampleName)
            Pair(hapid, newSampleIndex)
        }
        .toMap()

    val genotype = GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(numSamples, numSites)

    val pangenomeGenotypesBySample = vcfSamples.indices
        .map { rangeGenotypeTableAssemblies.genotypeAllSites(it) }
        .map { genotypes ->
            if (indelToMissing) {
                genotypes.map { genotype ->
                    val diploid = GenotypeTableUtils.getDiploidValues(genotype)
                    if (diploid[0] == NucleotideAlignmentConstants.GAP_ALLELE || diploid[0] == NucleotideAlignmentConstants.INSERT_ALLELE)
                        diploid[0] = GenotypeTable.UNKNOWN_ALLELE
                    if (diploid[1] == NucleotideAlignmentConstants.GAP_ALLELE || diploid[1] == NucleotideAlignmentConstants.INSERT_ALLELE)
                        diploid[1] = GenotypeTable.UNKNOWN_ALLELE
                    GenotypeTableUtils.getDiploidValue(diploid[0], diploid[1])
                }.toByteArray()
            } else {
                genotypes
            }
        }
        .toTypedArray()

    hapids.forEachIndexed { i, hapid ->
        if (hapid == ".") return@forEachIndexed
        val sampleIndex =
            pangenomeHapidToSample[hapid] ?: error("No pangenome sample found for hapid: $hapid at position: $pos")
        genotype.setBaseRangeForTaxon(i, 0, pangenomeGenotypesBySample[sampleIndex])
    }

    val taxon = TaxaListBuilder().addAll(imputedTable.samples).build()
    return GenotypeTableBuilder.getInstance(genotype.build(), rangeGenotypeTableAssemblies.positions(), taxon)

}

// include Site Min Count, Site Min Allele Freq and Bed file filters
fun filterGenotypeTable(
    genotype: GenotypeTable,
    siteMinCount: Int? = null,
    siteMinAlleleFreq: Double? = null,
    bedfile: String? = null
): GenotypeTable {

    // If no filters are set, return the original genotype table
    if (siteMinCount == null && siteMinAlleleFreq == null && bedfile == null) return genotype

    val builder = FilterSiteBuilderPlugin()

    siteMinCount?.let { builder.siteMinCount(it) }
    siteMinAlleleFreq?.let { builder.siteMinAlleleFreq(it) }
    bedfile?.let { builder.bedFile(it) }

    return builder.runPlugin(genotype)

}

fun runGLM(genotype: GenotypeTable, phenotype: Phenotype, populationStructre: Phenotype? = null): List<TableReport> {

    val inputDatums = mutableListOf<Datum>()
    inputDatums.add(Datum("genotype", genotype, null))
    inputDatums.add(Datum("phenotype", phenotype, null))
    populationStructre?.let { inputDatums.add(Datum("populationStructure", it, null)) }

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

fun writeVCF(genotype: GenotypeTable, filename: String) {
    ExportUtils.writeToVCF(genotype, filename, false)
}

fun writeTable(table: TableReport, filename: String) {
    TableReportUtils.saveDelimitedTableReport(table, "\t", File(filename))
}

// log10(exactP) =  pf(F-statistics, df1 = Numerator df, df2 = Denominator df, lower.tail = F, log.p = TRUE)/log(10)
fun addPValueColumn(table: TableReport): TableReport {

    val columnNames = table.tableColumnNames.map { it.toString() }

    // Trait   Marker  Chr     Pos     marker_F        p       marker_Rsq      add_F   add_p   dom_F
    // dom_p   marker_df       marker_MS       error_df        error_MS        model_df
    // model_MS        minorObs
    val fStatsIndex = columnNames.indexOf("marker_F")
    val df1Index = columnNames.indexOf("marker_df")
    val df2Index = columnNames.indexOf("error_df")

    val builder = TableReportBuilder.getInstance(table.tableTitle, arrayOf(columnNames + "exactP"))

    (0 until table.rowCount).forEach { row ->

    }

    return builder.build()

}