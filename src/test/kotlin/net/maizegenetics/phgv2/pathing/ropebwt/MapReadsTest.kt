package net.maizegenetics.phgv2.pathing.ropebwt

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.fail

class MapReadsTest {

    @Test
    fun testCliktParams() {
        val mapReads = MapReads()

        val noIndex = mapReads.test("--read-files test1.fq --output-dir testDir")
        assertEquals(1, noIndex.statusCode)
        assertEquals("Usage: map-reads [<options>]\n\n" +
                "Error: missing option --index\n", noIndex.stderr)

        val noReads = mapReads.test("--index testIndex --output-dir testDir")
        assertEquals(1, noReads.statusCode)
        assertEquals("Usage: map-reads [<options>]\n\n" +
                "Error: must provide one of --key-file, --read-files\n", noReads.stderr)

        val bothReadInputs = mapReads.test("--index testIndex --key-file testKeyFile --read-files test1.fq --output-dir testDir")
        assertEquals(1, bothReadInputs.statusCode)
        assertEquals("Usage: map-reads [<options>]\n\n" +
                "Error: option --key-file cannot be used with --read-files\n", bothReadInputs.stderr)

        val noOutputDir = mapReads.test("--index testIndex --read-files test1.fq")
        assertEquals(1, noOutputDir.statusCode)
        assertEquals("Usage: map-reads [<options>]\n\n" +
                "Error: missing option --output-dir\n", noOutputDir.stderr)
    }

    @Test
    fun testParseMem() {
        fail("Not yet implemented")
    }
}