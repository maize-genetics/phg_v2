package net.maizegenetics.phgv2.cli

import biokotlin.seqIO.NucSeqIO
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension.Companion.asmList
import net.maizegenetics.phgv2.utils.getChecksumForString
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue


@ExtendWith(TestExtension::class)
class FullPipelineIT {

    companion object {
        //Setup/download  files
        //The tempfileDir is already created by the IntegrationTestExtension

        @JvmStatic
        @BeforeAll
        fun setup() {
            resetDirs()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
//            //Delete the tempDir
//            resetDirs()
        }

        fun resetDirs() {
            File(TestExtension.tempDir).deleteRecursively()
            File(TestExtension.asmDir).deleteRecursively()
            File(TestExtension.testVCFDir).deleteRecursively()
            File(TestExtension.testMafDir).deleteRecursively()
            File(TestExtension.testInputFastaDir).deleteRecursively()
            File(TestExtension.testOutputFastaDir).deleteRecursively()
            File(TestExtension.testOutputGVCFDIr).deleteRecursively()
            File(TestExtension.testTileDBURI).deleteRecursively()



            File(TestExtension.tempDir).mkdirs()
            File(TestExtension.asmDir).mkdirs()
            File(TestExtension.testVCFDir).mkdirs()
            File(TestExtension.testMafDir).mkdirs()
            File(TestExtension.testInputFastaDir).mkdirs()
            File(TestExtension.testOutputFastaDir).mkdirs()
            File(TestExtension.testOutputGVCFDIr).mkdirs()
            File(TestExtension.testTileDBURI).mkdirs()
        }
    }
    @Test
    fun testFullPipeline() {
        //Run the full pipeline
        //Create environment

        val setupEnv = SetupEnvironment()
        setupEnv.test("--output-dir ${TestExtension.tempDir}")

        //Run InitDB
        val initdb = Initdb()
        initdb.test("--db-path ${TestExtension.testTileDBURI}")

        //Run CreateRanges
        val createRanges = CreateRanges()
        val createRangesResult = createRanges.test("--gff ${TestExtension.smallseqAnchorsGffFile} --output ${TestExtension.testBEDFile} --reference-file ${TestExtension.smallseqRefFile}")

        //Create the agc record:
        val agcCompress = AgcCompress()
        var agcResult = agcCompress.test("--fasta-list ${TestExtension.smallseqAssembliesListFile} --db-path ${TestExtension.testTileDBURI} --reference-file ${TestExtension.smallseqRefFile}")
        println(agcResult.output)

        println(createRangesResult.output)
        //Run BuildRefVCF
        val createRefVcf = CreateRefVcf()
        val createRefVcfResult = createRefVcf.test("--bed ${TestExtension.testBEDFile} --reference-file ${TestExtension.smallseqRefFile} --reference-name Ref -o ${TestExtension.testVCFDir}")
        println(createRefVcfResult.output)

        //Run Anchorwave
        val alignAssemblies = AlignAssemblies()
        val alignAssembliesResult = alignAssemblies.test(
            "--gff ${TestExtension.smallseqAnchorsGffFile} --reference-file ${TestExtension.smallseqRefFile} " +
                    "-a ${TestExtension.smallseqAssembliesListFile} -o ${TestExtension.testMafDir}"
        )
        println(alignAssembliesResult.output)



        //Run BuildMafVCF
        val createMafVCF = CreateMafVcf()
        val createMAFVCFResult = createMafVCF.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} " +
                "--reference-file ${TestExtension.smallseqRefFile} --maf-dir ${TestExtension.testMafDir} -o ${TestExtension.testVCFDir}")
        println(createMAFVCFResult.output)

        //Load All HVCFs into Tile DB
        val loadVCF = LoadVcf()
        val loadVCFResult = loadVCF.test("--vcf-dir ${TestExtension.testVCFDir} --db-path ${TestExtension.testTileDBURI}")
        println(loadVCFResult.output)

        //Pull out the HVCF from TileDB
        val exportHVCF = ExportHvcf()
        val exportHVCFRefResult = exportHVCF.test("--db-path ${TestExtension.testTileDBURI} --sample-names Ref,LineA,LineB -o ${TestExtension.testOutputGVCFDIr}")
        println(exportHVCFRefResult.output)

        //Run Fasta Generator for REF
//        val createFastaFromHvcf = CreateFastaFromHvcf()
//        val createFastaFromHvcfResult = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} --agc-file ${TestExtension.testAGCFile} --sample-name Ref -o ${TestExtension.testOutputRefFasta}")
//        println(createFastaFromHvcfResult.output)

//        buildAllFastas()
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val createFastaFromHvcfRefResult = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} --fasta-type haplotype --hvcf-dir ${TestExtension.testOutputGVCFDIr} -o ${TestExtension.testOutputHaplotypeFasta}")
        println(createFastaFromHvcfRefResult.output)

        //Open the output fasta
        //loop through the fasta file and take the sequence and compare its hash with the id in the header
        NucSeqIO(TestExtension.testOutputHaplotypeFasta).readAll().forEach { chr, seq ->
            val header = seq.id
            val headerParts = header.split(" ")
            val id = headerParts[0]
            val seqHash = getChecksumForString(seq.seq())
            assertEquals(id, seqHash, "Hashes do not match for $id")
        }


        //build a composite genome from the HVCFs

        createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} --fasta-type composite --hvcf-file ${TestExtension.testOutputGVCFDIr}/Ref.vcf -o ${TestExtension.testOutputFastaDir}/Ref_composite.fa")
        createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} --fasta-type composite --hvcf-file ${TestExtension.testOutputGVCFDIr}/LineA.vcf -o ${TestExtension.testOutputFastaDir}/LineA_composite.fa")
        createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} --fasta-type composite --hvcf-file ${TestExtension.testOutputGVCFDIr}/LineB.vcf -o ${TestExtension.testOutputFastaDir}/LineB_composite.fa")

        //Load in the initial HVCFs and check that each of their ids are represented in the composite genome


//        //Compare ref to input
        val refDiff =  compareFastaSeqs(TestExtension.smallseqRefFile, "${TestExtension.testOutputFastaDir}/Ref_composite.fa")
        assertTrue(refDiff < 0.0001, "Ref Fasta is not the same as input")
//        //For each asm
        for(asmName in asmList) {
            //Run Fasta Generator for ASM
//            createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} --agc-file ${TestExtension.testAGCFile} --sample-name ${asmName} -o ${TestExtension.testOutputFastaDir}${asmName}")
            //print the fastaNames
            println("Comparing ${asmName} fasta")
            //Compare asm to input
            val asmDiff =  compareFastaSeqs("${TestExtension.smallSeqInputDir}${asmName}.fa", "${TestExtension.testOutputFastaDir}${asmName}_composite.fa")
            assertTrue(asmDiff < 0.0001, "${asmName} Fasta is not the same as input")
        }


    }

    fun compareFastaSeqs(inputFasta: String, outputFasta: String): Double {
        //Compare the input and output fasta files
        //Return the percent difference

        //load the input fasta into a NucSeq
        val truthSeqChrMap = NucSeqIO(inputFasta).readAll()
        //load the output fasta into a NucSeq
        val generatedSeqChrMap = NucSeqIO(outputFasta).readAll()
        var diffCount = 0
        var length = 0
        //compare the two NucSeqs
        for(chr in truthSeqChrMap.keys) {
            val truthSeq = truthSeqChrMap[chr]
            val generatedSeq = generatedSeqChrMap[chr]
            if(truthSeq != null && generatedSeq != null) {

//                for(i in 0 until truthSeq.seq().length) {
                for(i in 0 until generatedSeq.seq().length) { //Need to use the generatedSeq here as we lose some info coming from MAF
                    if(truthSeq[i] != generatedSeq[i]) {
                        diffCount++
                    }
                    length++
                }
            }
        }

        return diffCount.toDouble()/length.toDouble()
    }

}