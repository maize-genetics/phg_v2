package net.maizegenetics.phgv2.cli

import biokotlin.seqIO.NucSeqIO
import org.junit.jupiter.api.extension.ExtendWith
import kotlin.test.Test
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension.Companion.asmList
import kotlin.test.Ignore
import kotlin.test.assertTrue


@ExtendWith(TestExtension::class)
class FullPipelineIT {

    companion object {
        //Setup/download  files
        //The tempfileDir is already created by the IntegrationTestExtension

    }

    @Ignore
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
        val result = createRanges.test("--gff ${TestExtension.smallseqAnchorsGffFile} --output ${TestExtension.testBEDFile}")
        println(result.output)
        //Run BuildRefVCF
        val createRefVcf = CreateRefVcf()
        createRefVcf.test("--bed ${TestExtension.testBEDFile} --reference ${TestExtension.testRefFasta} -o ${TestExtension.testVCFDir}")

        //Run Anchorwave
        val alignAssemblies = AlignAssemblies()
        alignAssemblies.test(
            "--gff ${TestExtension.smallseqAnchorsGffFile} --ref ${TestExtension.smallseqRefFile} " +
                    "-a ${TestExtension.smallseqAssembliesListFile} -o ${TestExtension.testMafDir}"
        )

        //Run BuildMafVCF
        val createMafVCF = CreateMafVcf()
        createMafVCF.test("--db-path ${TestExtension.testTileDBURI} --bed ${TestExtension.testBEDFile} " +
                "--reference ${TestExtension.testRefFasta} --maf-dir ${TestExtension.testMafDir} -o ${TestExtension.testVCFDir}")
        //Load All HVCFs into Tile DB
        val loadVCF = LoadVcf()
        loadVCF.test("--vcf-dir ${TestExtension.testVCFDir} --db-path ${TestExtension.testTileDBURI}")


        //Pull out the HVCF from TileDB
        TODO("Pull out the HVCF from TileDB")

        //Run Fasta Generator for REF
        val fastaGenerator = CreateFastaFromHvcf()
        fastaGenerator.test("--db-path ${TestExtension.testTileDBURI} --agc-file ${TestExtension.testAGCFile} --sample-name Ref -o ${TestExtension.testOutputRefFasta}")
        //Compare ref to input
        val refDiff =  compareFastaSeqs(TestExtension.testRefFasta, TestExtension.testOutputRefFasta)
        assertTrue(refDiff < 0.0001, "Ref Fasta is not the same as input")
        //For each asm
        for(asmName in asmList) {
            //Run Fasta Generator for ASM
            fastaGenerator.test("--db-path ${TestExtension.testTileDBURI} --agc-file ${TestExtension.testAGCFile} --sample-name ${asmName} -o ${TestExtension.testOutputFastaDir}${asmName}")
            //Compare asm to input
            val asmDiff =  compareFastaSeqs("${TestExtension.asmDir}${asmName}", "${TestExtension.testOutputFastaDir}${asmName}")
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
                //walk through each bp and compare
                for(i in 0 until truthSeq.size()) {
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