package net.maizegenetics.phgv2.scaffolding

import org.junit.jupiter.api.Test

class ReportMemAlignmentsTest {

    @Test
    fun testReportMemAlignments() {
        val memFile = "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/mem_l100.lm12lm05_lm12.bed"
        val outputFile = "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/lm12lm05lm12_alignments_3.txt"

        val reportMemAlignments = ReportMemAlignments()

//        reportMemAlignments.reportMemAlignments(memFile, outputFile, 2)
    }

    @Test
    fun test4MemAlignments() {
        //mem_l100.lm12lm05_lm05.bed
        //mem_l100.lm12lm05_lm07.bed
        //mem_l100.lm12lm05_lm12.bed
        //mem_l100.lm12lm05_lm23.bed

        val memFiles = listOf(
            "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/mem_l100.lm12lm05_lm05.bed",
            "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/mem_l100.lm12lm05_lm07.bed",
            "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/mem_l100.lm12lm05_lm12.bed",
            "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/mem_l100.lm12lm05_lm23.bed"
        )
        val outputFiles = listOf(
            "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/lm12lm05_lm05_alignments_2.txt",
            "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/lm12lm05_lm07_alignments_2.txt",
            "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/lm12lm05_lm12_alignments_2.txt",
            "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/mem100s/lm12lm05/lm12lm05_lm23_alignments_2.txt"
        )
        val gbsNames = listOf("lm05", "lm07", "lm12", "lm23")

        val reportMemAlignments = ReportMemAlignments()

        for( i in memFiles.indices) {
            println("Processing ${memFiles[i]} with output ${outputFiles[i]}")
            reportMemAlignments.reportMemAlignments(memFiles[i], outputFiles[i], gbsNames[i],2)
        }


    }
}