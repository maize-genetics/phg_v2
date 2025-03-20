package net.maizegenetics.phgv2.utils

import biokotlin.seq.NucSeq
import net.maizegenetics.phgv2.cli.CreateMafVcf
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class VCFConversionUtilsTest {

    @Test
    fun testBuildRegionStrings() {
        //Create a list of variantContext records with refBlocks and SNPs to resize and create boundaries
        val variantContexts = listOf(
            createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
                Position("chr1",5), Position("chr1",10),
                Position("chr1",5), Position("chr1",10),"+"),
            createSNPVC("B97", Position("chr1",11),Position("chr1",11), Pair("A", "T"),
                Position("chr1",11), Position("chr1",11),"+"),
            createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
                Position("chr1",12), Position("chr1",15),
                Position("chr1",12), Position("chr1",15),"+"),
            createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
                Position("chr1",19), Position("chr1",20),
                Position("chr1",19), Position("chr1",20),"+"),
            createSNPVC("B97", Position("chr1",21),Position("chr1",21), Pair("A", "T"),
                Position("chr1",21), Position("chr1",21),"+"),
            createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
                Position("chr1",22), Position("chr1",25),
                Position("chr1",22), Position("chr1",25),"+"),
        )

        val regionStrings = VCFConversionUtils.buildNewAssemblyRegions(7,23,variantContexts)
            .map { "${it.first.contig}:${it.first.position}-${it.second.position}" }

        //List<Pair<Position,Position>>
        assertEquals(2, regionStrings.size)
        assertEquals("chr1:7-15", regionStrings[0])
        assertEquals("chr1:19-23", regionStrings[1])

        // test condition where one VariantContext spans the entire reference range

        val oneVariantContext = listOf(createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",60), Position("chr1",300),
            Position("chr1",160), Position("chr1",400),"+"))

        val containedRegionStrings = VCFConversionUtils.buildNewAssemblyRegions(180, 350, oneVariantContext)
            .map { "${it.first.contig}:${it.first.position}-${it.second.position}" }

        assertEquals(1, containedRegionStrings.size)
        assertEquals("chr1:180-350", containedRegionStrings[0])
    }

    @Test
    fun testConvertVariantContextToPositionRange() {
        val refBlockVariant = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",15), Position("chr1",20),"+")

        val range = VCFConversionUtils.convertVariantContextToPositionRange(refBlockVariant)

        assertEquals(Position("chr1",15), range.first)
        assertEquals(Position("chr1",20), range.second)


        //test a reverse complimented variant
        val refBlockVariant2 = createRefRangeVC(mapOf("chr1" to NucSeq("A".repeat(100))),"B97",
            Position("chr1",5), Position("chr1",10),
            Position("chr1",20), Position("chr1",15),"+")

        val range2 = VCFConversionUtils.convertVariantContextToPositionRange(refBlockVariant2)

        assertEquals(Position("chr1",20), range2.first)
        assertEquals(Position("chr1",15), range2.second)

    }

}