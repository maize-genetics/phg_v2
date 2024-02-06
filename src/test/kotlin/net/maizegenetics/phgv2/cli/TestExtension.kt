package net.maizegenetics.phgv2.cli

import net.maizegenetics.phgv2.utils.setupDebugLogging
import org.junit.jupiter.api.extension.BeforeAllCallback
import org.junit.jupiter.api.extension.ExtensionContext
import java.io.File

class TestExtension : BeforeAllCallback {

    companion object {
        //the next line converts Windows \ to linux / in the user home path
        val userHome = System.getProperty("user.home").replace('\\', '/')
        val tempDir = "${userHome}/temp/phgv2Tests/tempDir/"
        val testVCFDir = "${tempDir}vcfDir/"
        val testMafDir = "${tempDir}mafDir/"
        val testInputFastaDir = "${tempDir}inputFastaDir/"
        val testOutputFastaDir = "${tempDir}outputFastaDir/"
        val testOutputGVCFDIr = "${tempDir}outputGVCFDir/"


        const val smallSeqInputDir = "data/test/smallseq/"

        const val smallseqLineAFile = "${smallSeqInputDir}LineA.fa"
        const val smallseqLineAMafFile = "${smallSeqInputDir}LineA.maf"
        const val smallSeqLineAGvcfFile = "${smallSeqInputDir}LineA.g.vcf"
        const val smallseqLineBFile = "${smallSeqInputDir}LineB.fa"
        const val smallseqLineBMafFile = "${smallSeqInputDir}LineB.maf"
        const val smallSeqLineBGvcfFile = "${smallSeqInputDir}LineB.g.vcf"
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

        const val exportGraphDir = "data/test/exportGraph/"
        const val exportGraphSingleSample = "${exportGraphDir}testSingleSampleHaplotypeGraph.vcf"
        const val exportGraphMultiSample = "${exportGraphDir}testMultipleFilesHaplotypeGraph.vcf"

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

        //Imputation test files
        val readMappingDir = "data/test/kmerReadMapping/"
        val testKmerIndex = "${tempDir}testKmerIndex.txt"
        val testReads = "${tempDir}testReads.fastq"
        val testOutputDir = "${tempDir}testOutputDir/"
        val testOutputReadMappingSingleEnd = "${testOutputDir}readMapping_single.txt"
        val testOutputReadMappingPairedEnd = "${testOutputDir}readMapping_paired.txt"
        val testKeyFile = "${readMappingDir}keyFile.txt"
        val testKeyFileNoHeader = "${readMappingDir}keyFileNoHeader.txt"
        val testKeyFileMissingFileName = "${readMappingDir}keyFileMissingFileName.txt"

        val smallSeqSimReads = "data/test/kmerReadMapping/simulatedReads/"
        val testLineASimReadsPrefix = "${smallSeqSimReads}LineA"
        val testLineBSimReadsPrefix = "${smallSeqSimReads}LineB"
        val testLineABSimReadsPrefix = "${smallSeqSimReads}LineA_LineB"

        val outputReadMappingLineA = "${testOutputDir}LineA_1_readMapping.txt"
        val outputReadMappingLineASingle = "${testOutputDir}LineA_2_readMapping.txt"
        val outputReadMappingLineB = "${testOutputDir}LineB_1_readMapping.txt"
        val outputReadMappingLineBSingle = "${testOutputDir}LineB_2_readMapping.txt"
        val outputReadMappingLineAB = "${testOutputDir}LineAB_1_readMapping.txt"


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