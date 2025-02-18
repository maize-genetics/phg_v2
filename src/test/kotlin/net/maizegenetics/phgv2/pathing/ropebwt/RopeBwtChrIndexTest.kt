package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.seqIO.NucSeqIO
import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import java.io.File
import java.io.FileNotFoundException
import kotlin.test.assertEquals
import kotlin.test.fail

class RopeBwtChrIndexTest {

    companion object {
        val tempTestDir = "${TestExtension.tempDir}ropebwtTest/"

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
//            resetDirs()
        }

        private fun resetDirs() {

            File(TestExtension.tempDir).deleteRecursively()
            File(TestExtension.testOutputFastaDir).deleteRecursively()
            File(TestExtension.testOutputDir).deleteRecursively()
            File(tempTestDir).deleteRecursively()

            File(TestExtension.tempDir).mkdirs()
            File(TestExtension.testOutputFastaDir).mkdirs()
            File(TestExtension.testOutputDir).mkdirs()
            File(tempTestDir).mkdirs()
        }
    }

    @Test
    fun testCliktParams() {

        val ropeBWTChrIndex = RopeBwtChrIndex()
        val noKeyfile = ropeBWTChrIndex.test("--output-dir outputDir --index-file-prefix indexFilePrefix --threads 3")
        assertEquals(1, noKeyfile.statusCode)
        assertEquals("Usage: rope-bwt-chr-index [<options>]\n\n" +
                "Error: missing option --keyfile\n", noKeyfile.stderr)


        val noOutputDir = ropeBWTChrIndex.test("--keyfile keyFile --index-file-prefix indexFilePrefix --threads 3")
        assertEquals(1, noOutputDir.statusCode)
        assertEquals("Usage: rope-bwt-chr-index [<options>]\n\n" +
                "Error: missing option --output-dir\n", noOutputDir.stderr)

        val noIndexFilePrefix = ropeBWTChrIndex.test("--keyfile keyFile --output-dir outputDir --threads 3")
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
        val ropeBwtChrIndex = RopeBwtChrIndex()
        val keyfile = "data/test/ropebwt/asm_keyfile.txt"
        val keyFileParsed = ropeBwtChrIndex.parseKeyFile(keyfile)
        assertEquals(2, keyFileParsed.size)
        assertEquals(Pair("data/test/smallseq/Ref.fa", "Ref"), keyFileParsed[0])
        assertEquals(Pair("data/test/smallseq/LineA.fa", "LineA"), keyFileParsed[1])

        val keyFileBad = "data/test/ropebwt/asm_bad_keyfile.txt"
        assertThrows<IllegalArgumentException> { ropeBwtChrIndex.parseKeyFile(keyFileBad) }

        val noKeyFile = ""
        assertThrows<FileNotFoundException>{ ropeBwtChrIndex.parseKeyFile(noKeyFile) }
    }

    @Test
    fun testProcessKeyFileRecord() {
        val ropeBwtChrIndex = RopeBwtChrIndex()
        val fastaFile = "data/test/smallseq/Ref.fa"
        val sampleName = "sample1"
        val renameFastaDir = tempTestDir
        val (renameFastaFile, contigLengthPairs) = ropeBwtChrIndex.processKeyFileRecord(fastaFile, sampleName, renameFastaDir)

        //check the output file
        val originalNucSeq = NucSeqIO(fastaFile).readAll()
        NucSeqIO(renameFastaFile).readAll().forEach { nucSeq ->
            val outputContigName = nucSeq.key.split("_")
            val contigName = outputContigName[0]
            val originalSeq = originalNucSeq[contigName]
            assertEquals(originalSeq!!.id, contigName)
            assertEquals(originalSeq.seq(), nucSeq.value.seq())
        }

        //check the contig length pairs
        val contigLengths = NucSeqIO(fastaFile).readAll().map { Pair("${it.key}_${sampleName}", it.value.seq().length) }
        assertEquals(contigLengths, contigLengthPairs)
    }

    @Test
    fun testRenameFastaSeqs() {
        val ropeBwtChrIndex = RopeBwtChrIndex()
        val fastaFile = "data/test/smallseq/Ref.fa"
        val sampleName = "sample1"
        val outputFileName = "$tempTestDir/RefRename.fa"
        val contigLengthPairs = mutableListOf<Pair<String, Int>>()
        bufferedWriter( outputFileName ).use { writer ->
            ropeBwtChrIndex.renameFastaSeqs(fastaFile, sampleName, writer, contigLengthPairs)
        }
        val originalNucSeq = NucSeqIO(fastaFile).readAll()

        //Check the output file
        NucSeqIO(outputFileName).readAll().forEach { nucSeq ->
            val outputContigName = nucSeq.key.split("_")
            val contigName = outputContigName[0]
            val originalSeq = originalNucSeq[contigName]
            assertEquals(originalSeq!!.id, contigName)
            assertEquals(originalSeq.seq(), nucSeq.value.seq())
        }
    }

    @Test
    fun testAddSeqToIndex() {
        fail("Not yet implemented")
    }

    @Test
    fun testBuildChrLengthFile() {
        val ropeBwtChrIndex = RopeBwtChrIndex()
        val indexFilePrefix = "$tempTestDir/testIndex"
        val contigLengthPairs = listOf(Pair("seq1", 100), Pair("seq2", 200), Pair("seq3", 300))
        ropeBwtChrIndex.buildChrLengthFile(indexFilePrefix, contigLengthPairs)

        //Load in the file and check the contents
        val outputFile = "$tempTestDir/testIndex.fmd.len.gz"
        val lines = bufferedReader(outputFile).readLines()
        assertEquals(3, lines.size)
        assertEquals("seq1\t100", lines[0])
        assertEquals("seq2\t200", lines[1])
        assertEquals("seq3\t300", lines[2])

    }


}