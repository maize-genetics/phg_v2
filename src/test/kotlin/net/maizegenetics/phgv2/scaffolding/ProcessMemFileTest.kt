package net.maizegenetics.phgv2.scaffolding

import org.junit.jupiter.api.Test

class ProcessMemFileTest {
    @Test
    fun testProcessMemFile() {
        val memFile = "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/memAlignments/mem.index1_query1.bed"
        val outputFile = "/Users/zrm22/Desktop/JuneHackathon2025/scaffolding/memAlignments/processMem_index1_query_lm05/"

        val processMemFile = ProcessMemFile()

        processMemFile.processMemFile(memFile, outputFile)
    }
}