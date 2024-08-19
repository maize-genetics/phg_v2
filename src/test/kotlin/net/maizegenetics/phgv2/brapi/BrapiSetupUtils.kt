package net.maizegenetics.phgv2.brapi

import com.github.ajalt.clikt.testing.test
import com.google.common.io.Files
import net.maizegenetics.phgv2.cli.*
import net.maizegenetics.phgv2.utils.bgzipAndIndexGVCFfile
import java.io.File
import kotlin.test.assertEquals

/**
 * This function creates a tiledb dataset from the smallseq assemblies and loads the smallseq vcf files into it.
 * This is needed for testing various brapi endpoint.
 * It will also create the AGC compressed files from the smallSeq fastas.
 */
fun createSmallSeqTiledb(tiledbURI:String = TestExtension.testTileDBURI) {
    val setupEnv = SetupEnvironment()
    println("createSmallSeqTiledb: calling setupEnv")
    setupEnv.test("--output-dir ${TestExtension.tempDir}")

    //Run InitDB
    println("createSmallSeqTiledb - calling Initdb")
    val initdb = Initdb()
    initdb.test("--db-path ${tiledbURI}")

    //Create the agc file:
    println("createSmallSeqTiledb - calling AgcCompress")
    val agcCompress = AgcCompress()
    var agcResult = agcCompress.test("--fasta-list ${TestExtension.smallseqAssembliesListFile} --db-path ${tiledbURI} --reference-file ${TestExtension.smallseqRefFile}")
    println(agcResult.output)

    // Run CreateRefVcf
    println("createSmallSeqTiledb - calling CreateRefVcf")
    val createRefVcf = CreateRefVcf()
    val refName = "Ref"
    val refUrl = TestExtension.refURL

    val ranges = "data/test/smallseq/anchors.bed"
    val genome = "data/test/smallseq/Ref.fa"

    val result = CreateRefVcf().test("--bed $ranges --reference-name $refName --reference-file $genome --reference-url ${refUrl} --db-path $tiledbURI")
    assertEquals(0, result.statusCode )

    //Load All HVCFs into Tile DB
    println("createSmallSeqTiledb - calling LoadVcf")
    loadVcfFiles(tiledbURI)

}

fun resetDirs() {
    File(TestExtension.tempDir).deleteRecursively()
    File(TestExtension.asmDir).deleteRecursively()
    File(TestExtension.testVCFDir).deleteRecursively()
    File(TestExtension.testMafDir).deleteRecursively()
    File(TestExtension.testInputFastaDir).deleteRecursively()
    File(TestExtension.testOutputFastaDir).deleteRecursively()
    File(TestExtension.testOutputGVCFDir).deleteRecursively()
    File(TestExtension.testTileDBURI).deleteRecursively()



    File(TestExtension.tempDir).mkdirs()
    File(TestExtension.asmDir).mkdirs()
    File(TestExtension.testVCFDir).mkdirs()
    File(TestExtension.testMafDir).mkdirs()
    File(TestExtension.testInputFastaDir).mkdirs()
    File(TestExtension.testOutputFastaDir).mkdirs()
    File(TestExtension.testOutputGVCFDir).mkdirs()
    File(TestExtension.testTileDBURI).mkdirs()
}

fun loadVcfFiles(tiledbURI:String = TestExtension.testTileDBURI) {
    // copy the files from data/test/smallseq to the tempDir
    // call bgzip to get the files we need for the test

    var origHvcfFile = "data/test/smallseq/Ref.h.vcf"
    var testHvcfFile = "${TestExtension.testVCFDir}Ref.h.vcf"
    Files.copy(File(origHvcfFile), File(testHvcfFile))

    // call bgzipAndIndexGVCFfile to zip and index the file
    bgzipAndIndexGVCFfile(testHvcfFile, "", false)

    // Copy LineA h.vcf
    origHvcfFile = "data/test/smallseq/LineA.h.vcf"
    testHvcfFile = "${TestExtension.testVCFDir}LineA.h.vcf"
    Files.copy(File(origHvcfFile), File(testHvcfFile))

    // call bgzipAndIndexGVCFfile to zip and index the file
    bgzipAndIndexGVCFfile(testHvcfFile, "", false)

    // Copy LineB h.vcf
    origHvcfFile = "data/test/smallseq/LineB.h.vcf"
    testHvcfFile = "${TestExtension.testVCFDir}LineB.h.vcf"
    Files.copy(File(origHvcfFile), File(testHvcfFile))

    // call bgzipAndIndexGVCFfile to zip and index the file
    bgzipAndIndexGVCFfile(testHvcfFile, "", false)

    // load the vcf files stored in the data/test/smallseq folder
    // This will load LineA.gvcf to the gvcf_dataset and Ref.hvcf to the hvcf_dataset
    println("loadVcfFiles - loading to tiledb")
    val loadVCF = LoadVcf()
    val vcfDir = TestExtension.testVCFDir
    var result = loadVCF.test("--vcf-dir ${vcfDir} --db-path ${tiledbURI} ")
    assertEquals(result.statusCode, 0)

    // Copy the reference fasta and bed file to folder TestExtension.testTileDBURI/reference
    // create the folder if it does not exist
    val refDir = "${TestExtension.testTileDBURI}/reference"
    File(refDir).mkdirs()
    var origRefFile = "data/test/smallseq/Ref.fa"
    var testRefFile = "${refDir}/Ref.fa"
    Files.copy(File(origRefFile), File(testRefFile))

    // Copy the bed file from data/test/smallseq/anchors.bed to the testTileDBURI/reference folder
    val origBedFile = "data/test/smallseq/anchors.bed"
    val testBedFile = "${refDir}/anchors.bed"
    Files.copy(File(origBedFile), File(testBedFile))


}