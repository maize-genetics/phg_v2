package net.maizegenetics.phgv2.cli

import org.junit.jupiter.api.extension.BeforeAllCallback
import org.junit.jupiter.api.extension.ExtensionContext
import java.io.File

class TestExtension : BeforeAllCallback {

    companion object {

        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"
        val testVCFDir = "${tempDir}vcfDir/"
        val testMafDir = "${tempDir}mafDir/"
        val testOutputFastaDir = "${tempDir}outputFastaDir/"
        val testOutputGVCFDIr = "${tempDir}outputGVCFDir/"

        const val smallseqLineAFile = "data/test/smallseq/LineA.fa"
        const val smallseqLineAMafFile = "data/test/smallseq/LineA.maf"
        const val smallseqLineBFile = "data/test/smallseq/LineB.fa"
        const val smallseqLineBMafFile = "data/test/smallseq/LineB.maf"
        const val smallseqRefFile = "data/test/smallseq/Ref.fa"
        const val smallseqAnchorsBedFile = "data/test/smallseq/anchors.bed"
        const val smallseqAnchorsGffFile = "data/test/smallseq/anchors.gff"
        const val smallseqGvcfFile = "data/test/smallseq/sample.gvcf"
        const val smallseqAssembliesListFile = "data/test/smallseq/assembliesList.txt"

        const val refLineName = "Ref"
        const val refFastaName = "Ref.fa"
        const val refURL = "https://s3.amazonaws.com/maizegenetics/phg/phgV2Test/Ref.fa" // this is a dummy URL
        val asmList = listOf("LineA.fa", "LineB.fa")


        val testTileDBURI = "${tempDir}testTileDBURI"
        const val testGFFFile = ""
        val testBEDFile = "${tempDir}testBEDFile.bed"
        val asmDir = "${tempDir}asmDir/"
        val testRefFasta = "${asmDir}${refFastaName}"

        val testAGCFile = "${tempDir}testAGCFile.agc"


        val testOutputRefFasta = "${testOutputFastaDir}${refFastaName}"

    }

    override fun beforeAll(context: ExtensionContext) {
        //Setup the test document environments

        File(tempDir).mkdirs()
        File(asmDir).mkdirs()
        File(testVCFDir).mkdirs()
        File(testMafDir).mkdirs()
        File(testOutputFastaDir).mkdirs()
        File(testOutputGVCFDIr).mkdirs()

    }
}