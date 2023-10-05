package net.maizegenetics.phgv2.cli

import biokotlin.seqIO.NucSeqIO
import org.junit.jupiter.api.extension.ExtendWith
import kotlin.test.Test
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.IntegrationTestExtension.Companion.asmList
import kotlin.test.Ignore
import kotlin.test.assertTrue


@ExtendWith(IntegrationTestExtension::class)
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

        TODO("Call CreateEnvironment command")

        //Run InitDB
        val initdb = Initdb()
        initdb.test("--db-path ${IntegrationTestExtension.testTileDBURI}")
        //Run CreateRanges
        val createRanges = CreateRanges()
        createRanges.test("--gff ${IntegrationTestExtension.testGFFFile} --output ${IntegrationTestExtension.testBEDFile}")
        //Run BuildRefVCF
        val buildRefVCF = BuildRefVcf()
        buildRefVCF.test("--bed ${IntegrationTestExtension.testBEDFile} --reference ${IntegrationTestExtension.testRefFasta} -o ${IntegrationTestExtension.testVCFDir}")

        //Run Anchorwave
        TODO("Run Anchorwave")

        //Run BuildMafVCF
        val buildMafVCF = BuildMafVcf()
        buildMafVCF.test("--maf-dir ${IntegrationTestExtension.testMafDir} -o ${IntegrationTestExtension.testVCFDir}")
        //Load All HVCFs into Tile DB
        val loadVCF = LoadVcf()
        loadVCF.test("--vcf-dir ${IntegrationTestExtension.testVCFDir} --db-path ${IntegrationTestExtension.testTileDBURI}")


        //Pull out the HVCF from TileDB
        TODO("Pull out the HVCF from TileDB")

        //Run Fasta Generator for REF
        val fastaGenerator = FastaGenerator()
        fastaGenerator.test("--db-path ${IntegrationTestExtension.testTileDBURI} --agc-file ${IntegrationTestExtension.testAGCFile} --sample-name Ref -o ${IntegrationTestExtension.testOutputRefFasta}")
        //Compare ref to input
        val refDiff =  compareFastaSeqs(IntegrationTestExtension.testRefFasta, IntegrationTestExtension.testOutputRefFasta)
        assertTrue(refDiff < 0.0001, "Ref Fasta is not the same as input")
        //For each asm
        for(asmName in asmList) {
            //Run Fasta Generator for ASM
            fastaGenerator.test("--db-path ${IntegrationTestExtension.testTileDBURI} --agc-file ${IntegrationTestExtension.testAGCFile} --sample-name ${asmName} -o ${IntegrationTestExtension.testOutputFastaDir}${asmName}")
            //Compare asm to input
            val asmDiff =  compareFastaSeqs("${IntegrationTestExtension.asmDir}${asmName}", "${IntegrationTestExtension.testOutputFastaDir}${asmName}")
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