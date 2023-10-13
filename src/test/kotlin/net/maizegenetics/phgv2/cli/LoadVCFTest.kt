package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import com.google.common.io.Files
import net.maizegenetics.phgv2.utils.VariantLoadingUtilsTest
import net.maizegenetics.phgv2.utils.bgzipAndIndexGVCFfile
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
        val resultMissingVCFDir = loadVCF.test("--db-path ${TestExtension.testTileDBURI} --temp-dir ${TestExtension.tempDir}")
        assertEquals(resultMissingVCFDir.statusCode, 1)
        assertEquals("Usage: load-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --vcf-dir: --vcf-dir must not be blank\n",resultMissingVCFDir.output)

        // Test missing db-path parameter
        val resultMissingDB = loadVCF.test("--vcf-dir ${TestExtension.testVCFDir} --temp-dir ${TestExtension.tempDir}")
        assertEquals(resultMissingDB.statusCode, 1)
        assertEquals("Usage: load-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --db-path: --db-path must not be blank\n",resultMissingDB.output)

        // Test missing temp-dir parameter
        val resultMissingTempDir = loadVCF.test("--vcf-dir ${TestExtension.testVCFDir} --db-path ${TestExtension.testTileDBURI}")
        assertEquals(resultMissingTempDir.statusCode, 1)
        assertEquals("Usage: load-vcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --temp-dir: --temp-dir must not be blank\n",resultMissingTempDir.output)
    }

    @Test
    fun testBadURIs() {

        val badPath = "${TestExtension.tempDir}badPATH"
        val URI = "${TestExtension.tempDir}hvcf_dataset"
        //Test non-existant dbPath
        assertThrows<IllegalStateException> {
            //Check that an error is thrown when the dbPath folder does not exist
            LoadVcf().verifyURI(badPath, URI)
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
            LoadVcf().verifyURI(goodPath, fileURI)
        }

        // Test tiledbURI is a directory, but not a tiledb array
        // delete the filePath created above
        File(filePath).delete()

        // Recreate it as a directory
        val dirURI = "${TestExtension.tempDir}hvcf_dataset"
        File(dirURI).mkdirs()
        assertThrows<IllegalStateException> {
            //Check that an error is thrown when the dbPath is good, but the URI is a regular directory file, not a tiledb array
            LoadVcf().verifyURI(goodPath, fileURI)
        }
    }

    @Test
    fun testGetFileLists() {
        // Change to TestExtension when that is updated to reference the files in data/test/smallseq

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
        assertEquals(lists.first.size, 1)
        assertEquals(lists.second.size, 1)
        assertEquals(lists.first[0], "${TestExtension.testVCFDir}LineA.g.vcf.gz")
        assertEquals(lists.second[0], "${TestExtension.testVCFDir}Ref.h.vcf.gz")

        // Check a folder with no g.vcf.gz or h.vcf.gz files
        // THe docs folder has no vcf files, so the lists should be empty
        vcfDir = "docs"
        val emptyLists = LoadVcf().getFileLists(vcfDir)
        assertEquals(emptyLists.first.size, 0)
        assertEquals(emptyLists.second.size, 0)

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
        assertEquals(dbDir.exists(), true)
        assertEquals(dbDir.isDirectory, true)
        assertEquals(dbDir.listFiles().size, 3) // includes a temp dir created for log file output
        assertTrue(dbDir.listFiles().map { it.name }.contains("hvcf_dataset"))
        assertTrue(dbDir.listFiles().map { it.name }.contains("gvcf_dataset"))

        // Test the code in LoadVCF.run() that exits when there are no files to load
        val loadVCF = LoadVcf()
        var vcfDir = "docs" // this folder has no vcf files

        // This call hits code that returns with an error message there are no vcf files to load
        assertThrows<IllegalArgumentException> {
            //Check that an error is thrown when the dbPath is good, but the URI is a regular directory file, not a tiledb array
            loadVCF.test("--vcf-dir ${vcfDir} --db-path ${dbPath} --temp-dir ${TestExtension.tempDir}")
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

        // now load the vcf files stored in the data/test/smallseq folder
        vcfDir = TestExtension.testVCFDir
        val result = loadVCF.test("--vcf-dir ${vcfDir} --db-path ${dbPath} --temp-dir ${TestExtension.tempDir}")

        // to verify the load worked, need to query tiledb for sample names and other data
        // This isn't going to work until we have means to query tiledb arrays.
        // This will be added in the next Pull Request

    }

}