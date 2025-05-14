package net.maizegenetics.phgv2.pathing.ropebwt

import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test
import com.github.ajalt.clikt.testing.test


class BuildSplineKnotsTest {
    @Test
    fun testCliktParams() {
        val buildSplineKnots = BuildSplineKnots()

        val noVcfDir = buildSplineKnots.test("--output-file dummyFile.json.gz")
        assertEquals(1, noVcfDir.statusCode)
        assertEquals("Usage: build-spline-knots [<options>]\n\n" +
                "Error: missing option --vcf-dir\n", noVcfDir.stderr)

        val noOutputFile = buildSplineKnots.test("--vcf-dir testDir")
        assertEquals(1, noOutputFile.statusCode)
        assertEquals("Usage: build-spline-knots [<options>]\n\n" +
                "Error: missing option --output-file\n", noOutputFile.stderr)
    }
}