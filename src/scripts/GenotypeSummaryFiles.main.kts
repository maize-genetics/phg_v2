#!/usr/bin/env kotlin

@file:DependsOn("net.maizegenetics:tassel:5.2.96")

import net.maizegenetics.analysis.data.GenotypeSummaryPlugin
import net.maizegenetics.dna.snp.ImportUtils
import net.maizegenetics.util.TableReportUtils
import java.io.File


val vcfInputDir = "/Users/terry/git/phgv2/data/test/smallseq"
val summaryOutputDir = "/Users/terry/git/phgv2/data/test/vcf-summary"

val vcfFiles = File(vcfInputDir)
    .walk()
    .map { it.absolutePath }
    .filter { it.endsWith(".vcf") || it.endsWith(".vcf.gz") }

vcfFiles
    .map { filename ->
        Pair(
            filename,
            ImportUtils.readFromVCF(filename, null, false, true)
        )
    }
    .filter { it.second != null }
    .forEach { (filename, genotypetable) ->
        println("Processing file: $filename")
        val outputFile = File(filename).name.substringBeforeLast(".vcf")
        val summaryTables = GenotypeSummaryPlugin(null, false).runPlugin(genotypetable)
        summaryTables.forEachIndexed { i, summaryTable ->
            val outputFileName = "$summaryOutputDir/${outputFile}_summary_$i.txt"
            TableReportUtils.saveDelimitedTableReport(summaryTable, File(outputFileName))
        }
    }