package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import org.apache.logging.log4j.LogManager
import java.io.File

class JoinImputedPathsWithAssemblyVCFs :
    CliktCommand(help = "Joins imputed paths with assembly VCFs, and runs GLM for each reference range.") {

    private val myLogger = LogManager.getLogger(JoinImputedPathsWithAssemblyVCFs::class.java)

    val vcfFilesPerRangeForAssembliesDir by option(help = "Directory containing VCF files per range for assemblies")
        .required()

    val traitFile by option(help = "Phenotype file")
        .required()

    val populationStructureFile by option(help = "Population structure file")
        .default("")

    val traitFile by option(help = "Phenotype file")
        .required()

    val traitFile by option(help = "Phenotype file")
        .required()

    override fun run() {

        // Create a map of key (i.e. chr10_154585427) to VCF file name (i.e. Zh-chr10_154585427-154627028.vcf)
        val vcfFileNamesPerRangeForAssemblies = File(vcfFilesPerRangeForAssembliesDir)
            .walk()
            .map { it.absolutePath }
            .filter { it.endsWith(".vcf") || it.endsWith(".vcf.gz") }
            .map { it.substringBeforeLast("-").substringAfter("Zh-") to it }
            .toMap()

        val phenotype = PhenotypeBuilder().fromFile(traitFile).build()[0]

        val populationStructure1 = populationStructureFile?.let { PhenotypeBuilder().fromFile(it).build()[0] }

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

}