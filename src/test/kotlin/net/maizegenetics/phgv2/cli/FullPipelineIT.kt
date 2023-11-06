package net.maizegenetics.phgv2.cli

import biokotlin.seqIO.NucSeqIO
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension.Companion.asmList
import net.maizegenetics.phgv2.utils.getChecksumForString
import org.junit.jupiter.api.AfterAll
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
        @AfterAll
        fun teardown() {
            //Delete the tempDir
            //File(TestExtension.tempDir).deleteRecursively()

            File(TestExtension.testTileDBURI).deleteRecursively()

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
        val createRangesResult = createRanges.test("--gff ${TestExtension.smallseqAnchorsGffFile} --output ${TestExtension.testBEDFile}")

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


        //Load in the initial HVCFs and check that each of their ids are represented in the composite genome


//        //Compare ref to input
//        val refDiff =  compareFastaSeqs(TestExtension.testRefFasta, TestExtension.testOutputRefFasta)
//        assertTrue(refDiff < 0.0001, "Ref Fasta is not the same as input")
//        //For each asm
//        for(asmName in asmList) {
//            //Run Fasta Generator for ASM
////            createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} --agc-file ${TestExtension.testAGCFile} --sample-name ${asmName} -o ${TestExtension.testOutputFastaDir}${asmName}")
//            //Compare asm to input
//            val asmDiff =  compareFastaSeqs("${TestExtension.asmDir}${asmName}", "${TestExtension.testOutputFastaDir}${asmName}")
//            assertTrue(asmDiff < 0.0001, "${asmName} Fasta is not the same as input")
//        }


    }

    fun exportAllHvcf() {
        val exportHVCF = ExportHvcf()
        val exportHVCFRefResult = exportHVCF.test("--db-path ${TestExtension.testTileDBURI} --sample-names Ref -o ${TestExtension.testOutputGVCFDIr}")
        println(exportHVCFRefResult.output)
        val exportHVCFLineAResult = exportHVCF.test("--db-path ${TestExtension.testTileDBURI} --sample-names LineA -o ${TestExtension.testOutputGVCFDIr}")
        println(exportHVCFLineAResult.output)
        val exportHVCFLineBResult = exportHVCF.test("--db-path ${TestExtension.testTileDBURI} --sample-names LineB -o ${TestExtension.testOutputGVCFDIr}")
        println(exportHVCFLineBResult.output)
    }

    fun buildAllFastas() {
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val createFastaFromHvcfRefResult = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} --fasta-type haplotype --hvcf-file ${TestExtension.testOutputGVCFDIr}/Ref.vcf -o ${TestExtension.testOutputRefFasta}")
        println(createFastaFromHvcfRefResult.output)
        val createFastaFromHvcfLineAResult = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} --fasta-type haplotype --hvcf-file ${TestExtension.testOutputGVCFDIr}/LineA.vcf -o ${TestExtension.testOutputFastaDir}LineA.fa")
        println(createFastaFromHvcfLineAResult.output)
        val createFastaFromHvcfLineBResult = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} --fasta-type haplotype --hvcf-file ${TestExtension.testOutputGVCFDIr}/LineB.vcf -o ${TestExtension.testOutputFastaDir}LineB.fa")
        println(createFastaFromHvcfLineBResult.output)

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