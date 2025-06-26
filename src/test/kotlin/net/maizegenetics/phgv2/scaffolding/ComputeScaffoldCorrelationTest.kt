package net.maizegenetics.phgv2.scaffolding

import biokotlin.util.bufferedReader
import org.jetbrains.letsPlot.core.spec.back.transform.bistro.util.map
import org.junit.jupiter.api.Test

class ComputeScaffoldCorrelationTest {

    @Test
    fun testComputeScaffoldCorrelation() {
        val collectedCountsDir = "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/CollectedCounts_2"
        val outputFile = "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/correlation_output_2.txt"

        val computeScaffoldCorrelation = ComputeScaffoldCorrelation()

        computeScaffoldCorrelation.computeScaffoldCorrelation(collectedCountsDir, outputFile, 0.9, CORRELATION_TYPE.PEARSON)
    }

    @Test
    fun testComputeScaffoldCorrelationSpearman() {
        val collectedCountsDir = "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/CollectedCounts_2"
        val outputFile = "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/correlation_output_spearman.txt"

        val computeScaffoldCorrelation = ComputeScaffoldCorrelation()

        computeScaffoldCorrelation.computeScaffoldCorrelation(collectedCountsDir, outputFile, 0.9, CORRELATION_TYPE.SPEARMAN)
    }

    @Test
    fun checkResults() {
        val file = "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/correlation_output_2.txt"

        val lines = bufferedReader(file).readLines().drop(1).map { it.split("\t") }
        println("Number of high Corr: ${lines.size}")

        var sameCount = 0
        var diffCount = 0

        val scaffoldSet = mutableSetOf<String>()


        for(line in lines) {
            val scaffold1 = line[0]
            val scaffold2 = line[1]
            val corr = line[2].toDouble()

            val sample1 = scaffold1.substringAfterLast("_")
            val sample2 = scaffold2.substringAfterLast("_")
            if (sample1 == sample2) {
                sameCount++
            } else {
                diffCount++
            }
            scaffoldSet.add(scaffold1)
            scaffoldSet.add(scaffold2)
        }

        println("Same sample count: $sameCount")
        println("Different sample count: $diffCount")

        println("Total unique scaffolds: ${scaffoldSet.size}")
        val groupedBySample = scaffoldSet.groupBy { it.substringAfterLast("_") }
        for ((sample, scaffolds) in groupedBySample) {
            println("Sample: $sample, Scaffolds: ${scaffolds.size}")
        }

    }

    @Test
    fun checkResultsSpearman() {
        val file = "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/correlation_output_spearman.txt"

        val lines = bufferedReader(file).readLines().drop(1).map { it.split("\t") }
        println("Number of high Corr: ${lines.size}")

        var sameCount = 0
        var diffCount = 0

        val scaffoldSet = mutableSetOf<String>()


        for(line in lines) {
            val scaffold1 = line[0]
            val scaffold2 = line[1]
            val corr = line[2].toDouble()

            val sample1 = scaffold1.substringAfterLast("_")
            val sample2 = scaffold2.substringAfterLast("_")
            if (sample1 == sample2) {
                sameCount++
            } else {
                diffCount++
            }
            scaffoldSet.add(scaffold1)
            scaffoldSet.add(scaffold2)
        }

        println("Same sample count: $sameCount")
        println("Different sample count: $diffCount")

        println("Total unique scaffolds: ${scaffoldSet.size}")
        val groupedBySample = scaffoldSet.groupBy { it.substringAfterLast("_") }
        for ((sample, scaffolds) in groupedBySample) {
            println("Sample: $sample, Scaffolds: ${scaffolds.size}")
        }

    }

}