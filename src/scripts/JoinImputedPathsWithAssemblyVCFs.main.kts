@file:DependsOn("net.maizegenetics:tassel:5.2.94")
@file:DependsOn("org.biokotlin:biokotlin:0.22")

import biokotlin.genome.Position
import biokotlin.util.bufferedReader
import net.maizegenetics.dna.snp.GenotypeTable
import net.maizegenetics.dna.snp.GenotypeTableBuilder
import net.maizegenetics.dna.snp.ImportUtils
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder
import net.maizegenetics.phenotype.Phenotype
import net.maizegenetics.phenotype.PhenotypeBuilder
import net.maizegenetics.taxa.TaxaListBuilder
import java.io.File
import java.util.*


val vcfFilesPerRangeForAssembliesDir = "vcf-by-range"
// Create a map of key (i.e. chr10_154585427) to VCF file name (i.e. Zh-chr10_154585427-154627028.vcf)
val vcfFileNamesPerRangeForAssemblies = File(vcfFilesPerRangeForAssembliesDir)
    .walk()
    .map { it.absolutePath }
    .filter { it.endsWith(".vcf") || it.endsWith(".vcf.gz") }
    .map { it.substringBeforeLast("-").substringAfter("Zh-") to it }
    .toMap()

val traitFile = "allClimateSoilGeoVariables.txt"
val phenotype = PhenotypeBuilder().fromFile(traitFile).build()[0]

// CHROM   POS     ID      REF     ALT  Samples...
// /local/workdir/wl748/merge_pangenome_v2.txt
val pangenomeTableFile = "merge_pangenome_v2.txt"
val pangenomeTable = readTable(pangenomeTableFile, 4)

// #CHROM  POS Samples...
// /local/workdir/wl748/merge_hapID/merge_SeeD.txt
val imputedTableFile = "merge_SeeD2.txt"
val imputedTable = readTable(imputedTableFile, 2)

imputedTable.posToLine.forEach { (pos, line) ->
    val genotypeTable = processRange(pos, line)
    runGLM(genotypeTable, phenotype)
}


fun processRange(pos: Position, line: String): GenotypeTable {

    val pangenomeLine = pangenomeTable.posToLine[pos] ?: error("No pangenome entry found for $pos")
    val pangenomeHapids = pangenomeLine.split("\t").drop(4)
    require(pangenomeHapids.size == pangenomeTable.samples.size) {
        "Number of hapids (${pangenomeHapids.size}) does not match number of samples (${pangenomeTable.samples.size}) at $pos"
    }
    val pangenomeHapidToSample = pangenomeHapids
        .mapIndexed { sampleIndex, hapid -> Pair(hapid, sampleIndex) }
        .toMap()

    val hapids = line.split("\t").drop(2)
    val numSamples = imputedTable.samples.size
    require(hapids.size == numSamples) {
        "Number of hapids (${hapids.size}) does not match number of samples ($numSamples) at $pos"
    }

    val key = "${pos.contig}_${pos.position}"
    val vcf = vcfFileNamesPerRangeForAssemblies[key] ?: error("No VCF file found for $key")
    val rangeGenotypeTableAssemblies = ImportUtils.readFromVCF(vcf, null, false, true)
    val numSites = rangeGenotypeTableAssemblies.numberOfSites()
    val vcfSamples = rangeGenotypeTableAssemblies.taxa()

    val genotype = GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(numSamples, numSites)

    val pangenomeGenotypesBySample = vcfSamples.indices
        .map { rangeGenotypeTableAssemblies.genotypeAllSites(it) }
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

fun runGLM(genotype: GenotypeTable, phenotype: Phenotype) {
    // Run GLM
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
