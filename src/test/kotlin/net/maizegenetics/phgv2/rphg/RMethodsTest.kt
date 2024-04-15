package net.maizegenetics.phgv2.rphg

import net.maizegenetics.phgv2.api.HaplotypeGraph
import org.junit.jupiter.api.Test
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.extension.ExtendWith
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class RMethodsTest {
    @Test
    fun testStringMatrixWithNames() {
        val testStringMatrix: RCharacterMatrix = MatrixWithNames(
            colNames = arrayOf("char_col1", "char_col2", "char_col3", "char_col4"),
            rowNames = arrayOf("row1", "row2", "row3"),
            matrixData = arrayOf(
                arrayOf("a", "b", "c"),
                arrayOf("d", "e", "f"),
                arrayOf("g", "h", "i"),
                arrayOf("j", "k", "l")
            )
        )

        assertEquals("phgv2_string_matrix", testStringMatrix.toString())
        assertEquals(4, testStringMatrix.colNames?.size)
        assertEquals(3, testStringMatrix.rowNames?.size)
        assertEquals(3, testStringMatrix.matrixData?.get(0)?.size)
        assertEquals(4, testStringMatrix.matrixData?.size)
        assertEquals(
            listOf("char_col1", "char_col2", "char_col3", "char_col4"),
            testStringMatrix.colNames?.toList()
        )
        assertEquals(
            listOf("row1", "row2", "row3"),
            testStringMatrix.rowNames?.toList()
        )
        assertEquals(
            "g",
            testStringMatrix.matrixData!![2][0]
        )
    }

    @Test
    fun testIntMatrixWithNames() {
        val testIntMatrix: RIntegerMatrix = MatrixWithNames (
            colNames = arrayOf("int_col1", "int_col2"),
            rowNames = arrayOf("row1", "row2", "row3"),
            matrixData = arrayOf(
                intArrayOf(1, 2, 3),
                intArrayOf(5, 6, 7)
            )
        )

        assertEquals("phgv2_int_matrix", testIntMatrix.toString())
        assertEquals(2, testIntMatrix.colNames?.size)
        assertEquals(3, testIntMatrix.rowNames?.size)
        assertEquals(3, testIntMatrix.matrixData?.get(0)?.size)
        assertEquals(2, testIntMatrix.matrixData?.size)
        assertEquals(
            listOf("int_col1", "int_col2"),
            testIntMatrix.colNames?.toList()
        )
        assertEquals(
            listOf("row1", "row2", "row3"),
            testIntMatrix.rowNames?.toList()
        )
        assertEquals(
            7,
            testIntMatrix.matrixData!![1][2]
        )
    }

    @Test
    fun testDoubleMatrixWithNames() {
        val testDblMatrix: RNumericMatrix = MatrixWithNames(
            arrayOf("dbl_col1", "dbl_col2"),
            arrayOf("row1", "row2", "row3", "row4"),
            arrayOf(
                doubleArrayOf(1.3, 2.1, 3.1434, 4.2),
                doubleArrayOf(5.1, 6.234, 7.0, 8.67)
            )
        )

        assertEquals("phgv2_dbl_matrix", testDblMatrix.toString())
        assertEquals(2, testDblMatrix.colNames?.size)
        assertEquals(4, testDblMatrix.rowNames?.size)
        assertEquals(4, testDblMatrix.matrixData?.get(0)?.size)
        assertEquals(2, testDblMatrix.matrixData?.size)
        assertEquals(
            listOf("dbl_col1", "dbl_col2"),
            testDblMatrix.colNames?.toList()
        )
        assertEquals(
            listOf("row1", "row2", "row3", "row4"),
            testDblMatrix.rowNames?.toList()
        )
        assertEquals(
            6.234,
            testDblMatrix.matrixData!![1][1]
        )
    }

    @Test
    fun testRList() {
        val testRList: RList = MatrixWithNames(
            arrayOf("rl_id1", "rl_id2"),
            arrayOf("row1", "row2", "row3"),
            arrayOf(
                arrayOf(1, 2, 3),
                arrayOf("this", "is", "categorical")
            )
        )

        assertEquals("phgv2_r_list", testRList.toString())
        assertEquals(2, testRList.colNames?.size)
        assertEquals(3, testRList.rowNames?.size)
        assertEquals(3, testRList.matrixData?.get(0)?.size)
        assertEquals(2, testRList.matrixData?.size)
        assertEquals(
            listOf("rl_id1", "rl_id2"),
            testRList.colNames?.toList()
        )
        assertEquals(
            listOf("row1", "row2", "row3"),
            testRList.rowNames?.toList()
        )
        assert(testRList.matrixData!![0].isArrayOf<Int>())
        assert(testRList.matrixData!![1].isArrayOf<String>())
    }

    @Test
    fun testGraphRefRangeRetrieval() {
        val graph = HaplotypeGraph(
            listOf(
                TestExtension.smallseqLineAHvcfFile,
                TestExtension.smallseqLineBHvcfFile,
                TestExtension.smallseqRefHvcfFile
            )
        )

        val rMethods = RMethods()
        val testRefRanges = rMethods.getRefRangesFromGraph(graph)

        assertEquals(
            listOf("seqname", "start", "end", "rr_id"),
            testRefRanges.colNames?.toList()
        )

        assertEquals("phgv2_r_list", testRefRanges.toString())
        assertEquals(40, testRefRanges.matrixData!![0].size)
        assertEquals(4, testRefRanges.matrixData!!.size)
        assert(testRefRanges.rowNames.isNullOrEmpty())
        assert(testRefRanges.matrixData!![0].isArrayOf<String>())
        assert(testRefRanges.matrixData!![1].isArrayOf<Int>())
        assert(testRefRanges.matrixData!![2].isArrayOf<Int>())
        assert(testRefRanges.matrixData!![3].isArrayOf<String>())
    }

    @Test
    fun testGraphAltHeaderRetrieval() {
        val graph = HaplotypeGraph(
            listOf(
                TestExtension.smallseqLineAHvcfFile,
                TestExtension.smallseqLineBHvcfFile,
                TestExtension.smallseqRefHvcfFile
            )
        )

        val rMethods = RMethods()
        val testAltHeader = rMethods.getAltHeadersFromGraph(graph)

        assertEquals(
            listOf(
                "hap_id",
                "sample_name",
                "description",
                "source",
                "checksum",
                "ref_range_hash"
            ),
            testAltHeader.colNames?.toList()
        )

        assertEquals("phgv2_r_list", testAltHeader.toString())
        assertEquals(116, testAltHeader.matrixData!![0].size)
        assertEquals(6, testAltHeader.matrixData!!.size)
        assert(testAltHeader.rowNames.isNullOrEmpty())
        (0..5).forEach {
            assert(testAltHeader.matrixData!![it].isArrayOf<String>())
        }
    }

    @Test
    fun testGraphAltHeaderPosRetrieval() {
        val graph = HaplotypeGraph(
            listOf(
                TestExtension.smallseqLineAHvcfFile,
                TestExtension.smallseqLineBHvcfFile,
                TestExtension.smallseqRefHvcfFile
            )
        )

        val rMethods = RMethods()
        val testAltPosHeader = rMethods.getAltHeaderPositionsFromGraph(graph)

        assertEquals(
            listOf(
                "hap_id",
                "contig_start",
                "contig_end",
                "start",
                "end"
            ),
            testAltPosHeader.colNames?.toList()
        )

        assertEquals("phgv2_r_list", testAltPosHeader.toString())
        assertEquals(116, testAltPosHeader.matrixData!![0].size)
        assertEquals(5, testAltPosHeader.matrixData!!.size)
        assert(testAltPosHeader.rowNames.isNullOrEmpty())

        (0..2).forEach {
            assert(testAltPosHeader.matrixData!![it].isArrayOf<String>())
        }

        (3..4).forEach {
            assert(testAltPosHeader.matrixData!![it].isArrayOf<Int>())
        }
    }

    @Test
    fun testGraphHapIdRetrieval() {
        val graph = HaplotypeGraph(
            listOf(
                TestExtension.smallseqLineAHvcfFile,
                TestExtension.smallseqLineBHvcfFile,
                TestExtension.smallseqRefHvcfFile
            )
        )

        val rMethods = RMethods()
        val testHapId = rMethods.getHapIdMatrixFromGraph(graph)

        assertEquals(
            listOf("LineA_G1", "LineB_G1", "Ref_G1"),
            testHapId.rowNames?.toList()
        )
        assertEquals("phgv2_string_matrix", testHapId.toString())
        assertEquals(40, testHapId.matrixData!![0].size)
        assertEquals(3, testHapId.matrixData!!.size)
        assertEquals("null", testHapId?.matrixData!![0][39])
        assert(testHapId?.matrixData!![0].isArrayOf<String>())
    }
}