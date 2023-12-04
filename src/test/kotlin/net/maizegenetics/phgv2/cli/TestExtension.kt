package net.maizegenetics.phgv2.cli

import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.junit.jupiter.api.extension.BeforeAllCallback
import org.junit.jupiter.api.extension.ExtensionContext
import java.io.File

class TestExtension : BeforeAllCallback {

    companion object {

        val tempDir = "${System.getProperty("user.home")}/temp/phgv2Tests/tempDir/"
        val testVCFDir = "${tempDir}vcfDir/"
        val testMafDir = "${tempDir}mafDir/"
        val testInputFastaDir = "${tempDir}inputFastaDir/"
        val testOutputFastaDir = "${tempDir}outputFastaDir/"
        val testOutputGVCFDIr = "${tempDir}outputGVCFDir/"

        const val smallSeqInputDir = "data/test/smallseq/"
        const val smallseqLineAFile = "${smallSeqInputDir}LineA.fa"
        const val smallseqLineAMafFile = "${smallSeqInputDir}LineA.maf"
        const val smallSeqLineAGvcfFile = "${smallSeqInputDir}LineA.g.vcf"
        const val smallseqLineBMafFile = "${smallSeqInputDir}LineB.maf"
        const val smallseqRefFile = "${smallSeqInputDir}Ref.fa"
        const val smallseqAnchorsBedFile = "${smallSeqInputDir}anchors.bed"
        const val smallseqAnchorsGffFile = "${smallSeqInputDir}anchors.gff"
        const val smallseqGvcfFile = "${smallSeqInputDir}sample.gvcf"
        const val smallseqAssembliesListFile = "${smallSeqInputDir}assembliesList.txt"
        const val smallseqRefHvcfFile = "${smallSeqInputDir}Ref.h.vcf"
        const val smallseqLineAHvcfFile = "${smallSeqInputDir}LineA.h.vcf"
        const val smallseqLineAHvcfFileBadAltTag = "${smallSeqInputDir}LineA_old_BadALTHeader.h.vcf"
        const val smallseqLineBHvcfFile = "${smallSeqInputDir}LineB.h.vcf"
        const val smallseqLineBHvcfFileBadAltTag = "${smallSeqInputDir}LineB_old_BadALTHeader.h.vcf"

        const val refLineName = "Ref"
        const val refFastaName = "Ref.fa"
        const val refURL = "https://s3.amazonaws.com/maizegenetics/phg/phgV2Test/Ref.fa" // this is a dummy URL
        val asmList = listOf("LineA", "LineB")


        val testTileDBURI = "${tempDir}testTileDBURI/"
        const val testGFFFile = ""
        val testBEDFile = "${tempDir}testBEDFile.bed"
        val asmDir = "${tempDir}asmDir/"
        val testRefFasta = "${asmDir}${refFastaName}"

        val testAGCFile = "${tempDir}testAGCFile.agc"


        val testOutputRefFasta = "${testOutputFastaDir}${refFastaName}"
        val testOutputHaplotypeFasta = "${testOutputFastaDir}haplotypes.fa"
        val testOutputCompositeFasta = "${testOutputFastaDir}composite.fa"

    }

    override fun beforeAll(context: ExtensionContext) {

        // Setup the test document environments

        setupDebugLogging()

        File(tempDir).mkdirs()
        File(asmDir).mkdirs()
        File(testVCFDir).mkdirs()
        File(testMafDir).mkdirs()
        File(testInputFastaDir).mkdirs()
        File(testOutputFastaDir).mkdirs()
        File(testOutputGVCFDIr).mkdirs()
        File(testTileDBURI).mkdirs()

    }
}