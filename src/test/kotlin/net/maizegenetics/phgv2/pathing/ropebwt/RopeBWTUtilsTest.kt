package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.pathing.ropebwt.RopeBwtIndexTest.Companion.indexFilePrefix
import net.maizegenetics.phgv2.pathing.ropebwt.RopeBwtIndexTest.Companion.inputFasta
import org.junit.jupiter.api.Test
import java.io.File

class RopeBWTUtilsTest {
    //runBuildStep(inputFasta:String, indexFilePrefix:String, numThreads: Int, condaEnvPrefix:String)
    @Test
    fun testRunBuildStep() {
        val numThreads = 3
        RopeBWTUtils.runBuildStep(inputFasta, indexFilePrefix, numThreads, "")

        //verify that the output files exist
        val fmdFile = File("$indexFilePrefix.fmr")
        assert(fmdFile.exists())
    }

    //convertBWTIndex(indexFilePrefix: String, condaEnvPrefix: String)
    //deleteFMRIndex(indexFilePrefix: String)
    @Test
    fun testConvertAndDeleteBWTIndex() {
        val ropebwtIndex = RopeBwtIndex()
        val numThreads = 3

        RopeBWTUtils.runBuildStep(inputFasta, indexFilePrefix, numThreads, "")

        RopeBWTUtils.convertBWTIndex(indexFilePrefix, "")
        val fmdFile = File("$indexFilePrefix.fmd")
        assert(fmdFile.exists())

        RopeBWTUtils.deleteFMRIndex(indexFilePrefix)
        val fmrFile = File("$indexFilePrefix.fmr")
        assert(!fmrFile.exists())
    }

    //buildSuffixArray(indexFilePrefix: String, numThreads: Int, condaEnvPrefix: String)
    @Test
    fun testBuildSuffixArray() {
        val ropebwtIndex = RopeBwtIndex()
        val numThreads = 3

        RopeBWTUtils.runBuildStep(inputFasta, indexFilePrefix, numThreads, "")
        RopeBWTUtils.convertBWTIndex(indexFilePrefix, "")

        RopeBWTUtils.buildSuffixArray(indexFilePrefix, numThreads, "")
        val saFile = File("$indexFilePrefix.fmd.ssa")
        assert(saFile.exists())
    }
}