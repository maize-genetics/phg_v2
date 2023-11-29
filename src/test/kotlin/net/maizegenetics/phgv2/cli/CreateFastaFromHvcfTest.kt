package net.maizegenetics.phgv2.cli

import biokotlin.seqIO.NucSeqIO
import com.github.ajalt.clikt.testing.test
import htsjdk.variant.vcf.VCFAltHeaderLine
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeader
import htsjdk.variant.vcf.VCFHeaderVersion
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.getChecksumForString
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
class CreateFastaFromHvcfTest {
    companion object {
        @JvmStatic
        @BeforeAll
        fun setup() {
            //create an AGC record with the Ref in it
            val altFileListFile = TestExtension.testOutputFastaDir+"/agc_altList.txt"
            BufferedWriter(FileWriter(altFileListFile)).use { writer ->
                writer.write("data/test/smallseq/LineA.fa\n")
                writer.write("data/test/smallseq/LineB.fa\n")
                writer.write("data/test/smallseq/LineC.fa\n")
            }

            val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"
            File(dbPath).mkdirs()

            //Call AGCCompress to create the AGC file
            val agcCompress = AgcCompress()
            agcCompress.processAGCFiles(dbPath,altFileListFile,"data/test/smallseq/Ref.fa")
        }
    }

    @Test
    fun testCliktParams() {
        val createFastaFromHvcf = CreateFastaFromHvcf()

        val resultMissingDB = createFastaFromHvcf.test("-o ${TestExtension.testOutputRefFasta} --fasta-type composite --hvcf-file /test_file.h.vcf")
        assertEquals(resultMissingDB.statusCode, 1)
        assertEquals("Usage: create-fasta-from-hvcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --db-path: --db-path must not be blank\n",resultMissingDB.output)

        val resultMissingOutput = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} --fasta-type haplotype --hvcf-file /test_file.h.vcf")
        assertEquals(resultMissingOutput.statusCode, 1)
        assertEquals("Usage: create-fasta-from-hvcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --output: --output/-o must not be blank\n",resultMissingOutput.output)

        val resultMissingFastaType = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} -o ${TestExtension.testOutputRefFasta} --hvcf-file /test_file.h.vcf")
        assertEquals(resultMissingFastaType.statusCode, 1)
        assertEquals("Usage: create-fasta-from-hvcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --fasta-type: --fasta-type must be either composite or haplotype\n",resultMissingFastaType.output)

        val resultBadFastaType = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} -o ${TestExtension.testOutputRefFasta} --fasta-type bad")
        assertEquals(resultBadFastaType.statusCode, 1)
        assertEquals("Usage: create-fasta-from-hvcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --fasta-type: --fasta-type must be either composite or haplotype\n",resultBadFastaType.output)
    }

    @Test
    fun testParseALTHeader() {
        val refHVCFFile = File("data/test/smallseq/Ref.h.vcf")
        val vcfReader = VCFFileReader(refHVCFFile, false)
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val altHeaders= createFastaFromHvcf.parseALTHeader(vcfReader.header)

        assertEquals(altHeaders.size, 40)

        //check a few to make sure they are correct
//        ##ALT=<ID=2b4590f722ef9229c15d29e0b4e51a0e,Description="haplotype data for line: Ref",Number=6,Source="data/test/smallseq/Ref.fa",Contig=1,Start=11001,End=12000,Checksum=Md5,RefRange=2b4590f722ef9229c15d29e0b4e51a0e>
        assertTrue(altHeaders.containsKey("2b4590f722ef9229c15d29e0b4e51a0e"))
        val currentHeader2b = altHeaders["2b4590f722ef9229c15d29e0b4e51a0e"]
        assertEquals(currentHeader2b?.id, "2b4590f722ef9229c15d29e0b4e51a0e")
        assertEquals(currentHeader2b?.description, "\"haplotype data for line: Ref\"")
        assertEquals(currentHeader2b?.number, "6")
        assertEquals(currentHeader2b?.source, "\"data/test/smallseq/Ref.fa\"")
        assertEquals(currentHeader2b?.regions?.get(0), Pair(Position("1",11001), Position("1",12000)))
        assertEquals(currentHeader2b?.checksum, "Md5")
        assertEquals(currentHeader2b?.refRange, "2b4590f722ef9229c15d29e0b4e51a0e")
//        ##ALT=<ID=db22dfc14799b1aa666eb7d571cf04ec,Description="haplotype data for line: Ref",Number=6,Source="data/test/smallseq/Ref.fa",Contig=2,Start=16501,End=17500,Checksum=Md5,RefRange=db22dfc14799b1aa666eb7d571cf04ec>
        assertTrue(altHeaders.containsKey("db22dfc14799b1aa666eb7d571cf04ec"))
        val currentHeaderdb = altHeaders["db22dfc14799b1aa666eb7d571cf04ec"]
        assertEquals(currentHeaderdb?.id, "db22dfc14799b1aa666eb7d571cf04ec")
        assertEquals(currentHeaderdb?.description, "\"haplotype data for line: Ref\"")
        assertEquals(currentHeaderdb?.number, "6")
        assertEquals(currentHeaderdb?.source, "\"data/test/smallseq/Ref.fa\"")
        assertEquals(currentHeaderdb?.regions?.get(0), Pair(Position("2",16501), Position("2",17500)))
        assertEquals(currentHeaderdb?.checksum, "Md5")
        assertEquals(currentHeaderdb?.refRange, "db22dfc14799b1aa666eb7d571cf04ec")

    //        ##ALT=<ID=5812acb1aff74866003656316c4539a6,Description="haplotype data for line: Ref",Number=6,Source="data/test/smallseq/Ref.fa",Contig=2,Start=1,End=1000,Checksum=Md5,RefRange=5812acb1aff74866003656316c4539a6>
        assertTrue(altHeaders.containsKey("5812acb1aff74866003656316c4539a6"))
        val currentHeader581 = altHeaders["5812acb1aff74866003656316c4539a6"]
        assertEquals(currentHeader581?.id, "5812acb1aff74866003656316c4539a6")
        assertEquals(currentHeader581?.description, "\"haplotype data for line: Ref\"")
        assertEquals(currentHeader581?.number, "6")
        assertEquals(currentHeader581?.source, "\"data/test/smallseq/Ref.fa\"")
        assertEquals(currentHeader581?.regions?.get(0), Pair(Position("2",1), Position("2",1000)))
        assertEquals(currentHeader581?.checksum, "Md5")
        assertEquals(currentHeader581?.refRange, "5812acb1aff74866003656316c4539a6")

        //Note, ID and description are both required by VCFAltHeaderLine so we do not need to check them.
        assertFailsWith<IllegalStateException>(
            message = "No exception found when Testing Number",
            block = {
                //,Regions=1:34001-38500
                createFastaFromHvcf.parseALTHeader(VCFHeader(setOf(VCFAltHeaderLine("<ID=id," +
                        "Description=\"haplotype data for line: testSample\">," +
                        "Number_bad=6," +
                        "Source=\"archive.agc\"," +
                        "Regions=\"1:200-300\"," +
                        "Checksum=\"Md5\"," +
                        "RefRange=\"id\">", VCFHeaderVersion.VCF4_2))))
            }
        )

        assertFailsWith<IllegalStateException>(
            message = "No exception found when Testing Source",
            block = {
                createFastaFromHvcf.parseALTHeader(VCFHeader(setOf(VCFAltHeaderLine("<ID=id," +
                        "Description=\"haplotype data for line: testSample\">," +
                        "Number=6," +
                        "Source_bad=\"archive.agc\"," +
                        "Regions=\"1:200-300\"," +
                        "Checksum=\"Md5\"," +
                        "RefRange=\"id\">", VCFHeaderVersion.VCF4_2))))
            }
        )

        assertFailsWith<IllegalStateException>(
            message = "No exception found when Testing checksum",
            block = {
                createFastaFromHvcf.parseALTHeader(VCFHeader(setOf(VCFAltHeaderLine("<ID=id," +
                        "Description=\"haplotype data for line: testSample\">," +
                        "Number=6," +
                        "Source=\"archive.agc\"," +
                        "Regions=\"1:200-300\"," +
                        "Checksum_bad=\"Md5\"," +
                        "RefRange=\"id\">", VCFHeaderVersion.VCF4_2))))
            }
        )
        assertFailsWith<IllegalStateException>(
            message = "No exception found when Testing refRange",
            block = {
                createFastaFromHvcf.parseALTHeader(VCFHeader(setOf(VCFAltHeaderLine("<ID=id," +
                        "Description=\"haplotype data for line: testSample\">," +
                        "Number=6," +
                        "Source=\"archive.agc\"," +
                        "Regions=\"1:200-300\"," +
                        "Checksum=\"Md5\"," +
                        "RefRange_bad=\"id\">", VCFHeaderVersion.VCF4_2))))
            }
        )
    }


    @Test
    fun testWriteCompositeSequence() {
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val outputFile = "${TestExtension.testOutputFastaDir}/testWriteCompositeSequence.fa"

        //Build 5 random sequences of 300bps each
        val rand = Random(12345)
        val alleles = listOf("A","C","G","T")

        val seqs = listOf(
            (1..300).map { alleles[rand.nextInt(4)] }.joinToString(""),
            (1..300).map { alleles[rand.nextInt(4)] }.joinToString(""),
            (1..300).map { alleles[rand.nextInt(4)] }.joinToString(""),
            (1..300).map { alleles[rand.nextInt(4)] }.joinToString(""),
            (1..300).map { alleles[rand.nextInt(4)] }.joinToString("")
        )

        //Build some simple Haplotype Sequences
        val haplotypeSequences =  listOf(
            HaplotypeSequence(getChecksumForString(seqs[0]), seqs[0], getChecksumForString(seqs[0]), "1", 1, 300, listOf(Pair(Position("1", 1), Position("1", 300)))),
            HaplotypeSequence(getChecksumForString(seqs[1]), seqs[1], getChecksumForString(seqs[1]), "1", 301, 600, listOf(Pair(Position("1", 301), Position("1", 600)))),
            HaplotypeSequence(getChecksumForString(seqs[2]), seqs[2], getChecksumForString(seqs[2]), "1", 601, 900, listOf(Pair(Position("1", 601), Position("1",900)))),
            HaplotypeSequence(getChecksumForString(seqs[3]), seqs[3], getChecksumForString(seqs[3]), "2", 1, 300, listOf(Pair(Position("2", 1), Position("1",300)))),
            HaplotypeSequence(getChecksumForString(seqs[4]), seqs[4], getChecksumForString(seqs[4]), "2", 301, 600, listOf(Pair(Position("2", 301), Position("1", 600))))
        )

        BufferedWriter(FileWriter(outputFile)).use { writer ->
            createFastaFromHvcf.writeCompositeSequence(writer, haplotypeSequences)
        }

        //Check that the file was created
        val file = File(outputFile)
        assertTrue(file.exists())

        //load in the file and check that the sequences are correct
        val importedSeqs = NucSeqIO(outputFile).readAll()

        assertEquals("${seqs[0]}${seqs[1]}${seqs[2]}", importedSeqs["1"]!!.seq())
        assertEquals("${seqs[3]}${seqs[4]}", importedSeqs["2"]!!.seq())
    }

    @Test
    fun testWriteHaplotypeSequence() {
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val outputFile = "${TestExtension.testOutputFastaDir}/testWriteHaplotypeSequence.fa"

        //Build 5 random sequences of 300bps each
        val rand = Random(12345)
        val alleles = listOf("A","C","G","T")

        val seqs = listOf(
            (1..300).map { alleles[rand.nextInt(4)] }.joinToString(""),
            (1..300).map { alleles[rand.nextInt(4)] }.joinToString(""),
            (1..300).map { alleles[rand.nextInt(4)] }.joinToString(""),
            (1..300).map { alleles[rand.nextInt(4)] }.joinToString(""),
            (1..300).map { alleles[rand.nextInt(4)] }.joinToString("")
        )

        //Build some simple Haplotype Sequences
        val haplotypeSequences =  listOf(
            HaplotypeSequence(getChecksumForString(seqs[0]), seqs[0], getChecksumForString(seqs[0]), "1", 1, 300, listOf(Pair(Position("1", 1), Position("1", 300)))),
            HaplotypeSequence(getChecksumForString(seqs[1]), seqs[1], getChecksumForString(seqs[1]), "1", 301, 600, listOf(Pair(Position("1", 301), Position("1", 600)))),
            HaplotypeSequence(getChecksumForString(seqs[2]), seqs[2], getChecksumForString(seqs[2]), "1", 601, 900, listOf(Pair(Position("1", 601), Position("1",900)))),
            HaplotypeSequence(getChecksumForString(seqs[3]), seqs[3], getChecksumForString(seqs[3]), "2", 1, 300, listOf(Pair(Position("2", 1), Position("1",300)))),
            HaplotypeSequence(getChecksumForString(seqs[4]), seqs[4], getChecksumForString(seqs[4]), "2", 301, 600, listOf(Pair(Position("2", 301), Position("1", 600))))
        )

        BufferedWriter(FileWriter(outputFile)).use { writer ->
            createFastaFromHvcf.writeHaplotypeSequence(writer, haplotypeSequences)
        }

        //Check that the file was created
        val file = File(outputFile)
        assertTrue(file.exists())

        //load in the file and check that the sequences are correct
        val importedSeqs = NucSeqIO(outputFile).readAll()

        assertEquals(seqs[0], importedSeqs[getChecksumForString(seqs[0])]!!.seq())
        assertEquals(seqs[1], importedSeqs[getChecksumForString(seqs[1])]!!.seq())
        assertEquals(seqs[2], importedSeqs[getChecksumForString(seqs[2])]!!.seq())
        assertEquals(seqs[3], importedSeqs[getChecksumForString(seqs[3])]!!.seq())
        assertEquals(seqs[4], importedSeqs[getChecksumForString(seqs[4])]!!.seq())

        //Check to make sure the metadata is correct
        compareFastaDescriptions(importedSeqs[getChecksumForString(seqs[0])]!!.description?:"", haplotypeSequences[0])
        compareFastaDescriptions(importedSeqs[getChecksumForString(seqs[1])]!!.description?:"", haplotypeSequences[1])
        compareFastaDescriptions(importedSeqs[getChecksumForString(seqs[2])]!!.description?:"", haplotypeSequences[2])
        compareFastaDescriptions(importedSeqs[getChecksumForString(seqs[3])]!!.description?:"", haplotypeSequences[3])
        compareFastaDescriptions(importedSeqs[getChecksumForString(seqs[4])]!!.description?:"", haplotypeSequences[4])


    }
    fun compareFastaDescriptions(description : String, hapSequence: HaplotypeSequence) {
        val descriptionParsed = description.replace(">","ID=")
            .split(" ")
            .map{ it.split("=") }
            .associate{ it[0] to it[1] }

        assertEquals(hapSequence.id, descriptionParsed["ID"])
        assertEquals(hapSequence.refRangeId, descriptionParsed["Ref_Range_Id"])
        assertEquals(hapSequence.refContig, descriptionParsed["Ref_Contig"])
        assertEquals(hapSequence.refStart, descriptionParsed["Ref_Start"]?.toInt()?:-1)
        assertEquals(hapSequence.refEnd, descriptionParsed["Ref_End"]?.toInt()?:-1)
        assertEquals("${hapSequence.asmRegions.first().first.contig}:${hapSequence.asmRegions.first().first.position}-${hapSequence.asmRegions.first().second.position}", descriptionParsed["Asm_Regions"])

    }

    @Test
    fun testCreateHaplotypeSequences() {
        val refHVCFFile = File("data/test/smallseq/Ref.h.vcf")
        val vcfReader = VCFFileReader(refHVCFFile, false)
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val altHeaders= createFastaFromHvcf.parseALTHeader(vcfReader.header)
        
        val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"

        val hapSequence = createFastaFromHvcf.createHaplotypeSequences(dbPath, "Ref", vcfReader.iterator().asSequence().toList(), altHeaders)

        assertEquals(40, hapSequence.size)
        val truthHashes = altHeaders.values.map { it.id }.toSet()

        //This verifies that we do indeed extract out the correct sequences
        hapSequence.forEach{
            assertTrue(truthHashes.contains(it.id))
        }
    }

    @Test
    fun testBuildFastaFromHVCF() {
        //buildFastaFromHVCF(dbPath: String, outputFile: String, fastaType:String, hvcfFile : String)
        val refHVCFFileName = "data/test/smallseq/Ref.h.vcf"
        val vcfReader = VCFFileReader(File(refHVCFFileName), false)
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val altHeaders= createFastaFromHvcf.parseALTHeader(vcfReader.header)

        val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"

        createFastaFromHvcf.buildFastaFromHVCF(dbPath, "${TestExtension.testOutputFastaDir}/Ref_Test_output.fa", "composite","",refHVCFFileName)

        //Compare the composite against the truth input
        val truthFasta = NucSeqIO("data/test/smallseq/Ref.fa").readAll()
        val outputFastaComposite = NucSeqIO("${TestExtension.testOutputFastaDir}/Ref_Test_output.fa").readAll()
        for(chr in truthFasta.keys) {
            assertEquals(getChecksumForString(truthFasta[chr]!!.seq()), getChecksumForString(outputFastaComposite[chr]!!.seq()))
        }

        //build a haplotype fasta as well
        createFastaFromHvcf.buildFastaFromHVCF(dbPath, "${TestExtension.testOutputFastaDir}/Ref_Test_output_hap.fa", "haplotype","",refHVCFFileName)

        //get truth hashes:
        val truthHashes = altHeaders.values.map { it.id }.toSet()
        //load in the haplotypes and build hashes
        val outputFastaHaplotypes = NucSeqIO("${TestExtension.testOutputFastaDir}/Ref_Test_output_hap.fa").readAll()
        outputFastaHaplotypes.values.map { getChecksumForString(it.seq()) }.toSet().forEach{
            assertTrue(truthHashes.contains(it))
        }




    }
}