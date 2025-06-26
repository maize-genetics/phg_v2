package net.maizegenetics.phgv2.scaffolding

import org.junit.jupiter.api.Test

class ProcessSamFileTest {

    @Test
    fun testProcessSamFile() {
        val samFile = "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/lm12_lm12_secondary.sam"
        val outputFile = "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/processSamFile/"

        val processSamFile = ProcessSamFile()

        processSamFile.processSAMFile(samFile, outputFile)
    }
}