package net.maizegenetics.phgv2.brapi.api

import net.maizegenetics.phgv2.brapi.createSmallSeqTiledb
import net.maizegenetics.phgv2.brapi.resetDirs
import net.maizegenetics.phgv2.brapi.service.VariantsService
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File

class VariantsTest {
    companion object {

        @JvmStatic
        @BeforeAll
        fun setup() {
            // delete, reset the directories
            resetDirs()

            // create the tiledb datasets, load them with from the vcf files
            // This will also create the AGC compressed file
            createSmallSeqTiledb()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            File(TestExtension.tempDir).deleteRecursively()
        }
    }

    @Test
    fun testCreateRefRangesForCache() {
        // THis test will verify the VarinetService.createRefRangesForCache() method
        // the createSmallSeqTiledb() should have created the tileDB folder, and copied
        // the bed file to it.  Use that for testing
        val bedFile = File("${TestExtension.testTileDBURI}/reference/").walk().filter { it.name.endsWith(".bed") }.toList()[0]
        val referenceRanges = VariantsService().createRefRangesForCache(bedFile)
        // Verify the number of lines in the file bedFile matches the number of reference ranges
        // in the referenceRanges list
        val bedFileLines = bedFile.readLines()
        assert(bedFileLines.size == referenceRanges.size)


    }
}