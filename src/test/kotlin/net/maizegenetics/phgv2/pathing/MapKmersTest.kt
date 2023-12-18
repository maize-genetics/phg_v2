package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File

@ExtendWith(TestExtension::class)
class MapKmersTest {

    companion object {
        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.testOutputDir).mkdirs()
        }
        @JvmStatic
        @AfterAll
        fun tearDown() {
            File(TestExtension.testOutputDir).deleteRecursively()
        }
    }

    @Test
    fun testCliktParams() {
        val mapKmers = MapKmers()

        val resultMissingKmerIndex =
            mapKmers.test("--hvcf-dir ${TestExtension.testVCFDir} --read-files ${TestExtension.testReads} --output-dir ${TestExtension.testOutputDir}")
        assertEquals(resultMissingKmerIndex.statusCode, 1)
        assertEquals(
            "Usage: map-kmers [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --kmer-index: --kmer-index must not be blank\n", resultMissingKmerIndex.output
        )

        val resultMissingReadsAndKeyFile =
            mapKmers.test("--hvcf-dir ${TestExtension.testVCFDir} --kmer-index ${TestExtension.testKmerIndex} --output-dir ${TestExtension.testOutputDir}")
        assertEquals(resultMissingReadsAndKeyFile.statusCode, 1)
        assertEquals(
            "Usage: map-kmers [<options>]\n" +
                    "\n" +
                    "Error: must provide one of --key-file, --read-files\n", resultMissingReadsAndKeyFile.output
        )

        val resultHavingBothReadsAndKeyFile =
            mapKmers.test("--hvcf-dir ${TestExtension.testVCFDir} --kmer-index ${TestExtension.testKmerIndex} --output-dir ${TestExtension.testOutputDir} --key-file ${TestExtension.testReads} --read-files ${TestExtension.testReads}")

        assertEquals(resultHavingBothReadsAndKeyFile.statusCode, 1)
        //This returns the same error message regardless of ordering between key-file and read-files
        assertEquals(
            "Usage: map-kmers [<options>]\n" +
                    "\n" +
                    "Error: option --key-file cannot be used with --read-files\n", resultHavingBothReadsAndKeyFile.output
        )


        val resultMissingOutputDir =
            mapKmers.test("--hvcf-dir ${TestExtension.testVCFDir} --kmer-index ${TestExtension.testKmerIndex} --read-files ${TestExtension.testReads}")
        assertEquals(resultMissingOutputDir.statusCode, 1)
        assertEquals(
            "Usage: map-kmers [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --output-dir: --output-dir/-o must not be blank\n", resultMissingOutputDir.output
        )

        val testMissingHVCFDir = mapKmers.test("--kmer-index ${TestExtension.testKmerIndex} --read-files ${TestExtension.testReads} --output-dir ${TestExtension.testOutputDir}")
        assertEquals(testMissingHVCFDir.statusCode, 1)
        assertEquals(
            "Usage: map-kmers [<options>]\n" +
                    "\n" +
                    "Error: invalid value for --hvcf-dir: --hvcf-dir must not be blank\n", testMissingHVCFDir.output
        )
    }


    @Test
    fun testImportKmerMap() {
//        val kmerIndexFile = "data/test/kmerReadMapping/SimpleIndex.txt"
//        val kmerMapData = importKmerIndex(kmerIndexFile)
//
//        val kmerMap = kmerMapData.kmerHashToLongMap
//
//
//        println(kmerMap.keys)

        //We should try to round trip a kmerIndex.  Write it out and then read it back in and compare the two.



    }

    @Test
    fun testParsingInputFileClasses() {
        //load in testKeyFile to a ReadInputFile.KeyFile object and make sure it is parsed correctly
        val keyFile = ReadInputFile.KeyFile(TestExtension.testKeyFile)
        val keyFileData = keyFile.getReadFiles()
        assertEquals(keyFileData.size, 3)
        //sampleName  filename    filename2
        //sample1 read1.txt   read2.txt
        //sample2 read3.txt   read4.txt
        //sample3 read5.txt
        assertEquals(keyFileData[0].sampleName, "sample1")
        assertEquals(keyFileData[0].file1, "read1.txt")
        assertEquals(keyFileData[0].file2, "read2.txt")

        assertEquals(keyFileData[1].sampleName, "sample2")
        assertEquals(keyFileData[1].file1, "read3.txt")
        assertEquals(keyFileData[1].file2, "read4.txt")

        assertEquals(keyFileData[2].sampleName, "sample3")
        assertEquals(keyFileData[2].file1, "read5.txt")
        assertEquals(keyFileData[2].file2, "")


        //Check that it breaks with a file with missing header
        val keyFileMissingHeader = ReadInputFile.KeyFile(TestExtension.testKeyFileNoHeader)
        val thrownExceptionMissingHeader = assertThrows<IllegalStateException> { keyFileMissingHeader.getReadFiles() }
        assertEquals("Key file ${TestExtension.testKeyFileNoHeader} must have a column named sampleName.", thrownExceptionMissingHeader.message)

        val keyFileMissingFileName = ReadInputFile.KeyFile(TestExtension.testKeyFileMissingFileName)
        val thrownExceptionMissingFileName = assertThrows<IllegalStateException> { keyFileMissingFileName.getReadFiles() }
        assertEquals("Key file ${TestExtension.testKeyFileMissingFileName} must have a column named filename.", thrownExceptionMissingFileName.message)


        //Test parsing of ReadInputFile.ReadFiles
        val readString = "file1.txt,file2.txt"
        val readFiles = ReadInputFile.ReadFiles(readString)
        val readFilesData = readFiles.getReadFiles()
        assertEquals(readFilesData.size, 1)
        assertEquals(readFilesData[0].sampleName, "noSample")
        assertEquals(readFilesData[0].file1, "file1.txt")
        assertEquals(readFilesData[0].file2, "file2.txt")

        //Test parsing of ReadInputFile.ReadFiles with only one file
        val readString2 = "file1.txt"
        val readFiles2 = ReadInputFile.ReadFiles(readString2)
        val readFilesData2 = readFiles2.getReadFiles()
        assertEquals(readFilesData2.size, 1)
        assertEquals(readFilesData2[0].sampleName, "noSample")
        assertEquals(readFilesData2[0].file1, "file1.txt")
        assertEquals(readFilesData2[0].file2, "")

        //Test parsing of ReadInputFile.ReadFiles with no files
        val readString3 = ""
        val readFiles3 = ReadInputFile.ReadFiles(readString3)

        val thrownExceptionNoReadFiles = assertThrows<IllegalStateException> { readFiles3.getReadFiles() }
        assertEquals("--read-files must have at least one file.", thrownExceptionNoReadFiles.message)

        //Test parsing of ReadInputFile.ReadFiles with more than 2 files
        val readString4 = "file1.txt,file2.txt,file3.txt"
        val readFiles4 = ReadInputFile.ReadFiles(readString4)
        val thrownExceptionTooManyReadFiles = assertThrows<IllegalStateException> { readFiles4.getReadFiles() }
        assertEquals("--read-files must have 1 or 2 files separated by commas.  You provided: 3", thrownExceptionTooManyReadFiles.message)

        //test different delimiter
        val readString5 = "file1.txt|file2.txt"
        val readFiles5 = ReadInputFile.ReadFiles(readString5)
        val readFilesData5 = readFiles5.getReadFiles()
        assertEquals(readFilesData5.size, 1)
        assertEquals(readFilesData5[0].sampleName, "noSample")
        assertEquals(readFilesData5[0].file1, "file1.txt|file2.txt")
        assertEquals(readFilesData5[0].file2, "")

    }

    @Test
    fun testRoundTrippingReadMappings() {
        //fun exportReadMapping(outputFileName: String, hapIdMapping: Map<List<String>, Int>, taxon : String, fastqFiles: Pair<String,String>) {

        val hapIdMapping = mapOf(listOf("hap1", "hap2") to 15, listOf("hap3") to 20, listOf("hap1") to 25, listOf("hap1", "hap2", "hap3") to 30,
            listOf("hap10", "hap11", "hap12") to 35, listOf("hap10", "hap11") to 40, listOf("hap12") to 45)

        val sampleName = "sample1"
        val fastqFilesSingle = Pair(TestExtension.testReads, "")
        val fastqFilesPaired = Pair(TestExtension.testReads, TestExtension.testReads)

        //test single end mapping
        exportReadMapping(TestExtension.testOutputReadMappingSingleEnd, hapIdMapping, sampleName, fastqFilesSingle)

        val readMappingSingleEnd = importReadMapping(TestExtension.testOutputReadMappingSingleEnd)
        assertEquals(readMappingSingleEnd.size, 7)
        assertEquals(readMappingSingleEnd, hapIdMapping)


        //test paired end mapping
        exportReadMapping(TestExtension.testOutputReadMappingPairedEnd, hapIdMapping, sampleName, fastqFilesPaired)
        val readMappingPairedEnd = importReadMapping(TestExtension.testOutputReadMappingPairedEnd)
        assertEquals(readMappingPairedEnd.size, 7)
        assertEquals(readMappingPairedEnd, hapIdMapping)

        //Open up the files and make sure they have the right headers
        val readMappingSingleEndFile = File(TestExtension.testOutputReadMappingSingleEnd).readLines().filter { it.startsWith("#") }.toSet()

        val truthHeaderSingleEndLines = setOf("#sampleName=sample1", "#filename1=${TestExtension.testReads}")
        for(line in truthHeaderSingleEndLines) {
            assert(readMappingSingleEndFile.contains(line))
        }

        val readMappingPairedEndFile = File(TestExtension.testOutputReadMappingPairedEnd).readLines().filter { it.startsWith("#") }.toSet()
        val truthHeaderPairedEndLines = setOf("#sampleName=sample1", "#filename1=${TestExtension.testReads}", "#filename2=${TestExtension.testReads}")
        for(line in truthHeaderPairedEndLines) {
            assert(readMappingPairedEndFile.contains(line))
        }
    }

}