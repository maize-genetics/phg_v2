package net.maizegenetics.phgv2.cli

import biokotlin.seqIO.NucSeqIO
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.phgv2.cli.TestExtension.Companion.asmList
import net.maizegenetics.phgv2.pathing.BuildKmerIndex
import net.maizegenetics.phgv2.pathing.FindPaths
import net.maizegenetics.phgv2.pathing.MapKmers
import net.maizegenetics.phgv2.utils.getBufferedWriter
import net.maizegenetics.phgv2.utils.getChecksumForString
import net.maizegenetics.phgv2.utils.retrieveAgcGenomes
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue


/**
 * This test will do a full end to end Integration test of the pipeline using the command line parameter inputs.
 * As of PHGv2.1 this will test the output sequences that they match the input assemblies(minus the parts which are not anchorwave aligned).
 */
@ExtendWith(TestExtension::class)
class FullPipelineIT {

    companion object {
        //Setup/download  files
        //Resetting on both setup and teardown just to be safe.
        @JvmStatic
        @BeforeAll
        fun setup() {
            resetDirs()
        }

        @JvmStatic
        @AfterAll
        fun teardown() {
            resetDirs()
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
//            File(TestExtension.testOutputDir).deleteRecursively()


            File(TestExtension.tempDir).mkdirs()
            File(TestExtension.asmDir).mkdirs()
            File(TestExtension.testVCFDir).mkdirs()
            File(TestExtension.testMafDir).mkdirs()
            File(TestExtension.testInputFastaDir).mkdirs()
            File(TestExtension.testOutputFastaDir).mkdirs()
            File(TestExtension.testOutputGVCFDIr).mkdirs()
            File(TestExtension.testTileDBURI).mkdirs()
            File(TestExtension.testOutputDir).mkdirs()
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
        val createRefVcfResult = createRefVcf.test("--bed ${TestExtension.testBEDFile} --reference-file ${TestExtension.smallseqRefFile} --reference-name Ref --db-path ${TestExtension.testTileDBURI}")
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
        val exportHVCF = ExportVcf()
        val exportHVCFRefResult = exportHVCF.test("--db-path ${TestExtension.testTileDBURI} --sample-names Ref,LineA,LineB -o ${TestExtension.testOutputGVCFDIr}")
        println(exportHVCFRefResult.output)

        //Create a fasta from the HVCF
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
        
        //Compare the outputs.
        val refDiff =  compareFastaSeqs(TestExtension.smallseqRefFile, "${TestExtension.testOutputFastaDir}/Ref_composite.fa")
        assertTrue(refDiff < 0.00001, "Ref Fasta is not the same as input")
        for(asmName in asmList) {
            println("Comparing ${asmName} fasta")
            //Compare asm to input
            val asmDiff =  compareFastaSeqs("${TestExtension.smallSeqInputDir}${asmName}.fa", "${TestExtension.testOutputFastaDir}${asmName}_composite.fa")
            assertTrue(asmDiff < 0.00001, "${asmName} Fasta is not the same as input")
        }

        //build a kmer index
        val buildKmerIndexArgs = "--agc-path ${TestExtension.testTileDBURI} --hvcf-dir ${TestExtension.testVCFDir}"
        val indexResult = BuildKmerIndex().test(buildKmerIndexArgs)
        assertEquals(0, indexResult.statusCode, "Kmer Indexing failed")
        println(indexResult.output)

        //create some reads with fairly high coverage to make sure mapping results are consistent so that they can be tested.
        val lineAFastqFilename = TestExtension.testInputFastaDir + "readsLineA.fastq"
        writeFastq(lineAFastqFilename, createHaploidReadsLineA())

        //map reads
        val mappingArgs = "--hvcf-dir ${TestExtension.testVCFDir} --kmer-index ${TestExtension.testVCFDir}kmerIndex.txt $ --read-files $lineAFastqFilename"
        val mapResult = MapKmers().test(mappingArgs)
        assertEquals(0, indexResult.statusCode, "Kmer Indexing failed")
        println(indexResult.output)

        //write a keyfile for FindPaths
        val pathKeyfile = "${TestExtension.testOutputDir}path_keyfile.txt"
        getBufferedWriter(pathKeyfile).use {myWriter ->
            myWriter.write("sampleName\treadMappingFiles\n")
            myWriter.write("TestSample\t$lineAFastqFilename\n")
        }

        //impute paths
        val pathArgs = "--path-keyfile $pathKeyfile --hvcf-dir ${TestExtension.testVCFDir} " +
                "--reference-genome ${TestExtension.smallseqRefFile} --output-dir ${TestExtension.testOutputDir} " +
                "--path-type haploid --prob-same-gamete 0.95"
        val pathResult = FindPaths().test(pathArgs)
        assertEquals(0, indexResult.statusCode, "Kmer Indexing failed")
        println(indexResult.output)

        //check paths


    }

    /**
     * Function to compare between two fasta sequences.
     * Because the alignment can lose information we cannot do a simple full identity check.
     * This counts number of differences/total number of bps in the chromosomes.
     *
     * The denominator is limited by the generated chromosome length because we use proali and it will lose the ends of the chromosomes until it hits an anchor.
     */
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

    private fun createHaploidReadsLineA(): List<String> {
        val coverage = 0.2
        val readLength = 100

        //create reads from lineA
        val readList = mutableListOf<String>()
        val seqMap = retrieveAgcGenomes(TestExtension.testTileDBURI, listOf("LineA"))
        seqMap.entries.forEach { (nameChr, seq) ->
            val chrlen = seq.size()
            val numberOfReads = (coverage * chrlen).toInt()
            repeat(numberOfReads) {
                val start = Random.nextInt(chrlen)
                val end = start + readLength
                readList.add(seq[start, end].seq())
            }
        }
        return readList
    }

    private fun writeFastq(filename: String, reads: List<String>) {
        getBufferedWriter(filename).use { myWriter ->
            var readCount = 1
            for (read in reads) {
                myWriter.write("@read$readCount\n")
                myWriter.write(read)
                myWriter.write("\n")
                myWriter.write("+\n")
                myWriter.write("x".repeat(read.length))
                myWriter.write("\n")
                readCount++
            }
        }
    }

}