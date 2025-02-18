package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.fail

class RopeBwtChrIndexTest {

    @Test
    fun testCliktParams() {

        val ropeBWTChrIndex = RopeBwtChrIndex()
        val noKeyFile = ropeBWTChrIndex.test("--output-dir outputDir --index-file-prefix indexFilePrefix --threads 3")
        assertEquals(1, noKeyFile.statusCode)
        assertEquals("Usage: rope-bwt-chr-index [<options>]\n\n" +
                "Error: missing option --key-file\n", noKeyFile.stderr)


        val noOutputDir = ropeBWTChrIndex.test("--key-file keyFile --index-file-prefix indexFilePrefix --threads 3")
        assertEquals(1, noOutputDir.statusCode)
        assertEquals("Usage: rope-bwt-chr-index [<options>]\n\n" +
                "Error: missing option --output-dir\n", noOutputDir.stderr)

        val noIndexFilePrefix = ropeBWTChrIndex.test("--key-file keyFile --output-dir outputDir --threads 3")
        assertEquals(1, noIndexFilePrefix.statusCode)
        assertEquals("Usage: rope-bwt-chr-index [<options>]\n\n" +
                "Error: missing option --index-file-prefix\n", noIndexFilePrefix.stderr)
    }

    @Test
    fun testCreateChrIndex() {
        fail("Not yet implemented")
    }

    @Test
    fun testParseKeyFile() {
        fail("Not yet implemented")
    }

    @Test
    fun testProcessKeyFileRecord() {
        fail("Not yet implemented")
    }

    @Test
    fun testRenameFastaSeqs() {
        fail("Not yet implemented")
    }

    @Test
    fun testAddSeqToIndex() {
        fail("Not yet implemented")
    }

    @Test
    fun testBuildChrLengthFile() {
        fail("Not yet implemented")
    }


}