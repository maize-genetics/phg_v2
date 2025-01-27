package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.seqIO.NucSeqIO
import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import java.io.File
import kotlin.test.assertEquals

class RopeBwtIndexTest {


    companion object {
        val tempTestDir = "${TestExtension.tempDir}ropebwtTest/"
        val tempDBPathDir = "${TestExtension.testOutputFastaDir}dbPath/"
        val inputFasta = TestExtension.smallSeqInputDir+"/pangenome/pangenome.fa"
        val indexFilePrefix = tempTestDir+"testIndex"


        //Setup/download  files
        //Resetting on both setup and teardown just to be safe.
        @JvmStatic
        @BeforeAll
        fun setup() {
            resetDirs()
            setupDebugLogging()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            resetDirs()
        }

        private fun resetDirs() {

            File(TestExtension.tempDir).deleteRecursively()
            File(TestExtension.testOutputFastaDir).deleteRecursively()
            File(TestExtension.testOutputDir).deleteRecursively()
            File(tempTestDir).deleteRecursively()
            File(tempDBPathDir).deleteRecursively()

            File(TestExtension.tempDir).mkdirs()
            File(TestExtension.testOutputFastaDir).mkdirs()
            File(TestExtension.testOutputDir).mkdirs()
            File(tempTestDir).mkdirs()
            File(tempDBPathDir).mkdirs()


        }


    }

    @Test
    fun testCliktParams() {
        val ropebwtIndex = RopeBwtIndex()


        //pangenome-file
        //hvcf-dir
        //Leave off the input files
        val noInputs = ropebwtIndex.test("--db-path ${TestExtension.tempDir} --output-dir ${TestExtension.tempDir}ropebwtTest/ --index-file-prefix testIndex --num-threads 3 --delete-fmr-index")
        assertEquals(1, noInputs.statusCode)
        assertEquals("Usage: ropebwt-index [<options>]\n\n" +
                "Error: must provide one of --pangenome-file, --hvcf-dir\n", noInputs.stderr)

        val tooManyInputs = ropebwtIndex.test("--pangenome-file data/test/smallseq/pangenome.fa --hvcf-dir data/test/smallseq/hvcf/ --db-path ${TestExtension.tempDir} --output-dir ${TestExtension.tempDir}ropebwtTest/ --index-file-prefix testIndex --num-threads 3 --delete-fmr-index")
        assertEquals(1, tooManyInputs.statusCode)
        assertEquals("Usage: ropebwt-index [<options>]\n\n" +
                "Error: option --pangenome-file cannot be used with --hvcf-dir\n", tooManyInputs.stderr)

        //output-dir
        val noOutputDir = ropebwtIndex.test("--pangenome-file data/test/smallseq/pangenome.fa --db-path ${TestExtension.tempDir} --index-file-prefix testIndex --num-threads 3 --delete-fmr-index")
        assertEquals(1, noOutputDir.statusCode)
        assertEquals("Usage: ropebwt-index [<options>]\n\n" +
                "Error: missing option --output-dir\n", noOutputDir.stderr)

        //index-file-prefix
        val noIndexFilePrefix = ropebwtIndex.test("--pangenome-file data/test/smallseq/pangenome.fa --db-path ${TestExtension.tempDir} --output-dir ${TestExtension.tempDir}ropebwtTest/ --num-threads 3 --delete-fmr-index")
        assertEquals(1, noIndexFilePrefix.statusCode)
        assertEquals("Usage: ropebwt-index [<options>]\n\n" +
                "Error: missing option --index-file-prefix\n", noIndexFilePrefix.stderr)

    }

    //runBuildStep(inputFasta:String, indexFilePrefix:String, numThreads: Int, condaEnvPrefix:String)
    @Test
    fun testRunBuildStep() {
        val ropebwtIndex = RopeBwtIndex()
        val numThreads = 3
        ropebwtIndex.runBuildStep(inputFasta, indexFilePrefix, numThreads, "")

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

        ropebwtIndex.runBuildStep(inputFasta, indexFilePrefix, numThreads, "")

        ropebwtIndex.convertBWTIndex(indexFilePrefix, "")
        val fmdFile = File("$indexFilePrefix.fmd")
        assert(fmdFile.exists())

        ropebwtIndex.deleteFMRIndex(indexFilePrefix)
        val fmrFile = File("$indexFilePrefix.fmr")
        assert(!fmrFile.exists())
    }

    //buildSuffixArray(indexFilePrefix: String, numThreads: Int, condaEnvPrefix: String)
    @Test
    fun testBuildSuffixArray() {
        val ropebwtIndex = RopeBwtIndex()
        val numThreads = 3

        ropebwtIndex.runBuildStep(inputFasta, indexFilePrefix, numThreads, "")
        ropebwtIndex.convertBWTIndex(indexFilePrefix, "")

        ropebwtIndex.buildSuffixArray(indexFilePrefix, numThreads, "")
        val saFile = File("$indexFilePrefix.fmd.ssa")
        assert(saFile.exists())
    }
    //buildChrLengthFile(inputFasta: String, indexFilePrefix: String)
    @Test
    fun testBuildChrLengthFile() {
        val ropebwtIndex = RopeBwtIndex()
        ropebwtIndex.buildChrLengthFile(inputFasta, indexFilePrefix)

        //verify that the output file exists
        val chrLengthFileName = "$indexFilePrefix.fmd.len.gz"
        val chrLengthFile = File(chrLengthFileName)
        assert(chrLengthFile.exists())

        //verify that the output file has the correct number of lines
        val chrLengthLines = bufferedReader(chrLengthFileName).readLines()
        val numLines =chrLengthLines.size
        val nucSeq = NucSeqIO(inputFasta).readAll()
        assertEquals(nucSeq.keys.size, numLines)

        //Check the counts
        val chrLengths = chrLengthLines
            .map { it.split("\t") }
            .map { Pair(it[0],it[1].toInt()) }
            .toMap()

        for(key in nucSeq.keys) {
            assertEquals(nucSeq[key]!!.seq().length, chrLengths[key])
        }
    }

    @Test
    fun createInitialIndex() {
        resetDirs()
        val ropebwtIndex = RopeBwtIndex()
        val numThreads = 3
        ropebwtIndex.createInitialIndex(inputFasta, indexFilePrefix, numThreads, true,"")

        //verify that the output files exist
        val fmdFile = File("$indexFilePrefix.fmd")
        assert(fmdFile.exists())

        val ssaFile = File("$indexFilePrefix.fmd.ssa")
        assert(ssaFile.exists())

        val chrLengthFile = File("$indexFilePrefix.fmd.len.gz")
        assert(chrLengthFile.exists())

        //verify that the fmr file was deleted
        val fmrFile = File("$indexFilePrefix.fmr")
        assert(!fmrFile.exists())

        //check chr length file
        val chrLengthLines = bufferedReader("$indexFilePrefix.fmd.len.gz").readLines()
        val numLines =chrLengthLines.size
        val nucSeq = NucSeqIO(inputFasta).readAll()
        assertEquals(nucSeq.keys.size, numLines)

        //Check the counts
        val chrLengths = chrLengthLines
            .map { it.split("\t") }
            .map { Pair(it[0],it[1].toInt()) }
            .toMap()

        for(key in nucSeq.keys) {
            assertEquals(nucSeq[key]!!.seq().length, chrLengths[key])
        }

    }

}