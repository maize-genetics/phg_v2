package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import com.google.common.io.Files
import net.maizegenetics.phgv2.utils.bgzipAndIndexGVCFfile
import net.maizegenetics.phgv2.utils.verifyURI
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.assertThrows
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.Test
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class LoadVCFTest {
    companion object {
        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"

        @JvmStatic
        @BeforeAll
        fun setup() {
            File(tempDir).mkdirs()
        }

        // Comment out the tearDown()if you need to look at the logs files created by ProcessBuilder
        // commands.
        @JvmStatic
        @AfterAll
        fun teardown() {
            File(tempDir).deleteRecursively()
        }
    }

    @Test
    fun testCliktParams() {
        val loadVCF = LoadVcf()

        // Test missing vcf-dir parameter
        val resultMissingVCFDir = loadVCF.test("--db-path ${TestExtension.testTileDBURI} ")
        assertEquals(resultMissingVCFDir.statusCode, 1)
        assertEquals("Usage: load-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --vcf-dir: --vcf-dir must not be blank\n",resultMissingVCFDir.output)

        // Test missing db-path parameter
        val resultMissingDB = loadVCF.test("--vcf-dir ${TestExtension.testVCFDir} ")
        assertEquals(resultMissingDB.statusCode, 1)
        assertEquals("Usage: load-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --db-path: --db-path must not be blank\n",resultMissingDB.output)

    }

    @Test
    fun testBadURIs() {

        val badPath = "${TestExtension.tempDir}badPATH"
        val URI = "${TestExtension.tempDir}hvcf_dataset"
        //Test non-existant dbPath
        assertThrows<IllegalStateException> {
            //Check that an error is thrown when the dbPath folder does not exist
            verifyURI(badPath, URI)
        }

        // Test tiledb URI is a file, not a tiledb array
        val goodPath = "${TestExtension.tempDir}"
        val fileURI = "hvcf_dataset"
        val filePath = "${TestExtension.tempDir}${fileURI}"

        File(filePath).bufferedWriter().use {
            it.write("this is just a text file")
        }
        assertThrows<IllegalArgumentException> {
            //Check that an error is thrown when the dbPath is good, but the URI is a regular file, not a tiledb array
            verifyURI(goodPath, fileURI)
        }

        // Test tiledbURI is a directory, but not a tiledb array
        // delete the filePath created above
        File(filePath).delete()

        // Recreate it as a directory
        val dirURI = "${TestExtension.tempDir}hvcf_dataset"
        File(dirURI).mkdirs()
        assertThrows<IllegalArgumentException> {
            //Check that an error is thrown when the dbPath is good, but the URI is a regular directory file, not a tiledb array
            verifyURI(goodPath, fileURI)
        }
    }

    @Test
    fun testGetFileLists() {

        // copy the files from data/test/smallseq to the tempDir
        // call bgzip to get the files we need for the test
        var origGvcfFile = "data/test/smallseq/LineA.g.vcf"
        var testGvcfFile = "${TestExtension.testVCFDir}LineA.g.vcf"
        Files.copy(File(origGvcfFile), File(testGvcfFile))

        val origHvcfFile = "data/test/smallseq/Ref.h.vcf"
        val testHvcfFile = "${TestExtension.testVCFDir}Ref.h.vcf"
        Files.copy(File(origHvcfFile), File(testHvcfFile))

        // call bgzipAndIndexGVCFfile to zip and index the file
        bgzipAndIndexGVCFfile(testGvcfFile)
        bgzipAndIndexGVCFfile(testHvcfFile)

        // Folder with both g.vcf.gz and h.vcf.gz files, and other files
        // verify the g.vcf.gz is returned in lists.first, and h.vcf.gz is returned in lists.second
        var vcfDir = TestExtension.testVCFDir
        val lists = LoadVcf().getFileLists(vcfDir)
        assertEquals(1,lists.first.size)
        assertEquals(1,lists.second.size)
        assertEquals(lists.first[0], "${TestExtension.testVCFDir}LineA.g.vcf.gz")
        assertEquals(lists.second[0], "${TestExtension.testVCFDir}Ref.h.vcf.gz")

        // Check a folder with no g.vcf.gz or h.vcf.gz files
        // THe docs folder has no vcf files, so the lists should be empty
        vcfDir = "docs"
        val emptyLists = LoadVcf().getFileLists(vcfDir)
        assertEquals(0, emptyLists.first.size)
        assertEquals(0, emptyLists.second.size)

    }

    @Test
    fun testLoadVCF() {
        // first create the tiledb arrays by calling Initdb
        val dbPath = "${TestExtension.testTileDBURI}"
        // make the dbPath directory if it does not exist
        File(dbPath).mkdirs()
        Initdb().createDataSets(dbPath)

        // verify the dbPath directory exists with subdirectories hvcf_dataset and gvcf_dataset
        val dbDir = File(dbPath)
        assertEquals(true, dbDir.exists())
        assertEquals(true, dbDir.isDirectory)
        assertEquals(3, dbDir.listFiles().size) // includes a temp dir created for log file output
        assertTrue(dbDir.listFiles().map { it.name }.contains("hvcf_dataset"))
        assertTrue(dbDir.listFiles().map { it.name }.contains("gvcf_dataset"))

        // Test the code in LoadVCF.run() that exits when there are no files to load
        val loadVCF = LoadVcf()
        var vcfDir = "docs" // this folder has no vcf files

        // This call hits code that returns with an error message there are no vcf files to load
        assertThrows<IllegalArgumentException> {
            //Check that an error is thrown when the dbPath is good, but the URI is a regular directory file, not a tiledb array
            loadVCF.test("--vcf-dir ${vcfDir} --db-path ${dbPath} ")
        }

        // copy the files from data/test/smallseq to the tempDir
        // call bgzip to get the files we need for the test
        var origGvcfFile = "data/test/smallseq/LineA.g.vcf"
        var testGvcfFile = "${TestExtension.testVCFDir}LineA.g.vcf"
        Files.copy(File(origGvcfFile), File(testGvcfFile))

        val origHvcfFile = "data/test/smallseq/Ref.h.vcf"
        val testHvcfFile = "${TestExtension.testVCFDir}Ref.h.vcf"
        Files.copy(File(origHvcfFile), File(testHvcfFile))

        // call bgzipAndIndexGVCFfile to zip and index the file
        bgzipAndIndexGVCFfile(testGvcfFile)
        bgzipAndIndexGVCFfile(testHvcfFile)

        // load the vcf files stored in the data/test/smallseq folder
        // This will load LineA.gvcf to the gvcf_dataset and Ref.hvcf to the hvcf_dataset
        vcfDir = TestExtension.testVCFDir
        var result = loadVCF.test("--vcf-dir ${vcfDir} --db-path ${dbPath} ")
        assertEquals(0,result.statusCode)

        // get the list of samples from tiledb gvcf_dataset
        var uri = "${dbPath}/gvcf_dataset"
        val sampleList = LoadVcf().getTileDBSampleLists(uri)

        // verify the samples are in the tiledb array
        assertEquals(1,sampleList.size)
        assertEquals("LineA", sampleList[0])

        // get the list of samples from the tiledb hvcf_dataset
        uri = "${dbPath}/hvcf_dataset"
        val sampleList2 = LoadVcf().getTileDBSampleLists(uri)

        //verify the samples are in the tiledb array
        assertEquals(1,sampleList2.size)
        assertEquals("Ref",sampleList2[0])

        // Verify the hvcf file was copied to dbPath/hvcf_files folder
        val hvcfFile = "${dbPath}/hvcf_files/Ref.h.vcf.gz"
        assertEquals(true, File(hvcfFile).exists())

        // Export to hvcf and verify the number of lines in the vcf file
        val exportHvcf = ExportVcf()
        // datatype defaults to hvcf
        val exportResultHvcf = exportHvcf.test("--db-path ${dbPath} --sample-names Ref -o ${TestExtension.tempDir}")
        assertEquals(0, exportResultHvcf.statusCode)

        // Tiledb writes the file to the output folder with a name of ${sampleName}.vcf
        // Verify the tiledb exported file has the same number of lines as does the
        // original vcf file.
        val exportFileHvcf = "${TestExtension.tempDir}Ref.vcf"
        val origFileH = "data/test/smallseq/Ref.h.vcf"
        val origFileLinesH = File(origFileH).readLines().size
        val exportFileLines = File(exportFileHvcf).readLines().size
        assertEquals(origFileLinesH, exportFileLines)

        // Export the gvcf file, verify the number of lines in the exported vcf file matches the number in the original vcf file
        val exportResultGvcf = exportHvcf.test("--db-path ${dbPath} --sample-names LineA --dataset-type gvcf -o ${TestExtension.tempDir}")
        assertEquals(0, exportResultGvcf.statusCode )

        val exportFileGvcf = "${TestExtension.tempDir}LineA.vcf"
        val origFileG = "data/test/smallseq/LineA.g.vcf"
        val origFileLinesG = File(origFileG).readLines().size
        val exportFileLinesG = File(exportFileGvcf).readLines().size
        assertEquals(origFileLinesG, exportFileLinesG)

        // Load the same files again, verify no error message,
        result = loadVCF.test(  "--vcf-dir ${vcfDir} --db-path ${dbPath} ")
        assertEquals(0, result.statusCode )

        // Export the hvcf and gvcf files again, verify the line numbers
        // are still the same as the original vcf files (ie, not increased due to duplicate entries)
        val exportResultH2 = exportHvcf.test("--db-path ${dbPath} --sample-names Ref --dataset-type hvcf -o ${TestExtension.tempDir}")
        assertEquals(exportResultH2.statusCode, 0)
        val exportFileHvcf2 = "${TestExtension.tempDir}Ref.vcf"
        val exportFileLinesH2 = File(exportFileHvcf2).readLines().size
        assertEquals(origFileLinesH, exportFileLinesH2) // same original hvcf file for comparision

        // Check the gvcf file lines
        val exportResultG2 = exportHvcf.test("--db-path ${dbPath} --sample-names LineA --dataset-type gvcf -o ${TestExtension.tempDir}")
        assertEquals(0, exportResultG2.statusCode )

        val exportFileGvcf2 = "${TestExtension.tempDir}LineA.vcf"
        val exportFileLinesG2 = File(exportFileGvcf2).readLines().size
        assertEquals(origFileLinesG, exportFileLinesG2) // same original gvcf file for comparision, number of lines should not have increased


    }

}