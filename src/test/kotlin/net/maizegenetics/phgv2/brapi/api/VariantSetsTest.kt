package net.maizegenetics.phgv2.brapi.api

import net.maizegenetics.phgv2.brapi.createSmallSeqTiledb
import net.maizegenetics.phgv2.brapi.resetDirs
import net.maizegenetics.phgv2.brapi.service.VariantSetsService
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.assertEquals

class VariantSetsTest {

    companion object {

        @JvmStatic
        @BeforeAll
        fun setup() {
            // delete, reset the directories
            resetDirs()

            // create the tiledb datasets, load them with from the vcf files
            // This will also create the AGC compressed file
            createSmallSeqTiledb()

            VariantSetsService
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(TestExtension.tempDir).deleteRecursively()
        }

    }

    @Test
    fun testVariantSetsService() {
        val filename = VariantSetsService.createAllSamplesHVCF()
        assertEquals(VariantSetsService.allSamplesHvcf, filename)
    }

}