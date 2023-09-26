package net.maizegenetics.phgv2.main

import biokotlin.featureTree.Genome
import org.junit.jupiter.api.Test
import kotlin.test.assertFails
import kotlin.test.assertEquals

class TestCreateRanges {
    @Test
    fun testIdMinMaxBounds() {
        val testGffPath = "data/testing/zm_b73v5_test.gff3.gz"
        val genes = Genome.fromGFF(testGffPath).genes()

        val cr = CreateRanges()

        val obsIdList01 = cr.idMinMaxBounds(genes, "gene", 0)
        val obsIdList02 = cr.idMinMaxBounds(genes, "cds", 0)
        val obsIdList03 = cr.idMinMaxBounds(genes, "gene", 100)

        val obsBedList01 = cr.generateBedRows(obsIdList01, genes)
        val obsBedList02 = cr.generateBedRows(obsIdList01, genes, ",")
        val obsBedList03 = cr.generateBedRows(obsIdList02, genes, featureId = "biotype")

        assertEquals(2, obsIdList01.size)
        assertEquals(obsIdList01[0], Pair(34616, 40203)) // should be 34616
        assertEquals(obsIdList01[1], Pair(41213, 46761))

        assertEquals(2, obsIdList02.size)
        assertEquals(obsIdList02[0], Pair(34721, 38365))
        assertEquals(obsIdList02[1], Pair(41526, 45912))

        assertEquals(2, obsIdList03.size)
        assertEquals(obsIdList03[0], Pair(34516, 40303))
        assertEquals(obsIdList03[1], Pair(41113, 46861))

        assertEquals(2, obsBedList01.size)
        assertEquals(obsBedList01[0], "chr1\t34616\t40203\tZm00001eb000010\t0\t+")
        assertEquals(obsBedList02[0], "chr1,34616,40203,Zm00001eb000010,0,+")
        assertEquals(obsBedList03[0], "chr1\t34721\t38365\tprotein_coding\t0\t+")

        assertFails {
            cr.idMinMaxBounds(genes, "geeeene", 0)
        }
    }
}
