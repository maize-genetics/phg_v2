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


        val refFastaName = "Ref.fa"
        val asmList = listOf("LineA.fa", "LineB.fa")


        val testTileDBURI = "${tempDir}testTileDBURI"
        val testGFFFile = ""
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
    }




}