package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class Hvcf2GvcfTest {
    companion object {

        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.testVCFDir).mkdirs()
            File(TestExtension.testTileDBURI).mkdirs()
            Initdb().createDataSets(TestExtension.testTileDBURI,"")
        }

        @JvmStatic
        @AfterAll
        fun tearDown() {
            File(TestExtension.tempDir).deleteRecursively()
            File(TestExtension.testVCFDir).deleteRecursively()
            File(TestExtension.testTileDBURI).deleteRecursively()

        }
    }

    @Test
    fun testSimpleHvcf2Gvcf() {
        // This is a basic test.  We copy an hvcf file created from CreateMafVcf to a new location
        // and run that through hvcf2gvcf.  We then verify that the output file exists.

        val dbPath = TestExtension.testTileDBURI
        val refName = "Ref"
        val refUrl = TestExtension.refURL

        val ranges = "data/test/smallseq/anchors.bed"
        val refFasta = "data/test/smallseq/Ref.fa"


        //Run InitDB
        println("createSmallSeqTiledb - calling Initdb")
        val initdb = Initdb()
        initdb.test("--db-path ${dbPath}")

        // Should run PrepareAssemblies ??  Or are they fine for smallSeq?
        // SmallSeq fastas already contain sampleName=<sampleName> in the header.

        // Tiledb datasets were created during initialization (@beforeAll)
        println("running agcCompress")
        val agcCompress = AgcCompress()
        var agcResult = agcCompress.test("--fasta-list ${TestExtension.smallseqAssembliesListFile} --db-path ${dbPath} --reference-file ${TestExtension.smallseqRefFile}")
        println(agcResult.output)

        println("running CreateRefVcf")
        var result = CreateRefVcf().test("--bed $ranges --reference-name $refName --reference-file $refFasta --reference-url ${refUrl} --db-path $dbPath")
        assertEquals(0, result.statusCode )

        // Verify the tiledbUri/reference folder exists and contains the ranges file
        val referenceDir = "${dbPath}/reference/"

        // RUn alignAssemblies test to get MAF files.
        // Run createMAFVCf on the assemblies LineA and LIneB to get
        // data into the db.

        println("running AlignAssemblies")
        val alignAssemblies = AlignAssemblies()

        result = alignAssemblies.test(
            "--gff ${TestExtension.smallseqAnchorsGffFile} --reference-file ${TestExtension.smallseqRefFile} " +
                    "--assembly-file-list ${TestExtension.smallseqAssembliesListFile} -o ${TestExtension.tempDir} --total-threads 1 --in-parallel 1"
        )

        println("testRunningAlignAssemblies: result output: ${result.output}")
        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        // Load assemblies using CreateMafVcf - creates and loads gvcf and hvcf
        // remember - there is no gvcf for the ref
        println("running CreateMafVcf")
        val createMafVcf = CreateMafVcf()
        result = createMafVcf.test("--db-path ${dbPath} --bed data/test/smallSeq/anchors.bed --reference-file ${refFasta} --maf-dir ${TestExtension.tempDir} -o ${TestExtension.testVCFDir}")
        println(result.output)

        // Need to load the vcf now!
        // Load the vcf files into the tiledb dataset
        println("running LoadVcf")
        val loadVcf = LoadVcf()
        result = loadVcf.test("--db-path ${dbPath} --vcf-dir ${TestExtension.testVCFDir}")

        // verify the gvcf and hvcf files were created, and that the hvcf file was written somewhere -
        // are we doing that by default?  if so, where?  Is it undewr the reference folder?
        // RUn to here, stop, verify the output files.  It is written to the tempDir/vcfDir (${TestExtension.testVCFDir})
        // so copy LineB.h.vcf.gz to LineBPath.h.vcf.gz
        // then run hvcf2gvcf on LineBPath.h.vcf.gz

        // Copy the $TestExtension.testVCFDir/LineB.h.vcf.gz to $TestExtension.testVCFDir/LineBPath.h.vcf.gz
        // Make directory ${TestExtension.testVCFDir}/testOutputGVCFDir

        println("NOW ... running hvcf2gvcf")
        val testGVCFdir = "${TestExtension.testVCFDir}/testOutputGVCFDir"
        File(testGVCFdir).mkdirs()
        val lineBPathHvcf = "${testGVCFdir}/LineBPath.h.vcf.gz"
        val lineBHvcf = "${TestExtension.testVCFDir}/LineB.h.vcf.gz"
        File(lineBHvcf).copyTo(File(lineBPathHvcf))

        // Run hvcf2gvcf on the copied file
        val hvcf2gvcf = Hvcf2Gvcf()
        result = hvcf2gvcf.test("--db-path ${dbPath} --hvcf-dir $testGVCFdir --output-dir ${testGVCFdir} --reference-file ${refFasta}")

        println("done !!")


    }

}