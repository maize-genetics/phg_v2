package net.maizegenetics.phgv2.cli

import biokotlin.seq.NucSeq
import biokotlin.seqIO.NucSeqIO
import com.github.ajalt.clikt.testing.test
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFAltHeaderLine
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeader
import htsjdk.variant.vcf.VCFHeaderVersion
import net.maizegenetics.phgv2.api.SampleGamete
import net.maizegenetics.phgv2.utils.*
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
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
            agcCompress.processAGCFiles(dbPath,altFileListFile,"data/test/smallseq/Ref.fa","")
        }
    }

    @Test
    fun testCliktParams() {
        val createFastaFromHvcf = CreateFastaFromHvcf()

        val resultMissingOutput = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} --fasta-type haplotype --hvcf-file /test_file.h.vcf")
        assertEquals(resultMissingOutput.statusCode, 1)
        assertEquals("Usage: create-fasta-from-hvcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --output-dir: --output-dir/-o must not be blank\n",resultMissingOutput.output)

        val resultMissingFastaType = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} -o ${TestExtension.tempDir} --hvcf-file /test_file.h.vcf")
        assertEquals(resultMissingFastaType.statusCode, 1)
        assertEquals("Usage: create-fasta-from-hvcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --fasta-type: --fasta-type must be either composite or haplotype\n",resultMissingFastaType.output)

        val resultBadFastaType = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} -o ${TestExtension.tempDir} --fasta-type bad")
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
        val altHeaders= parseALTHeader(vcfReader.header)

        assertEquals(altHeaders.size, 40)

        //check a few to make sure they are correct
//        ##ALT=<ID=2b4590f722ef9229c15d29e0b4e51a0e,Description="haplotype data for line: Ref",Number=6,Source="data/test/smallseq/Ref.fa",Contig=1,Start=11001,End=12000,Checksum=Md5,RefRange=2b4590f722ef9229c15d29e0b4e51a0e>
        assertTrue(altHeaders.containsKey("2b4590f722ef9229c15d29e0b4e51a0e"))
        val currentHeader2b = altHeaders["2b4590f722ef9229c15d29e0b4e51a0e"]
        assertEquals(currentHeader2b?.id, "2b4590f722ef9229c15d29e0b4e51a0e")
        assertEquals(currentHeader2b?.description, "haplotype data for line: Ref")
        assertEquals(currentHeader2b?.source, "data/test/smallseq/Ref.fa")
        assertEquals(currentHeader2b?.regions?.get(0), Pair(Position("1",11001), Position("1",12000)))
        assertEquals(currentHeader2b?.checksum, "Md5")
        assertEquals(currentHeader2b?.refRange, "2b4590f722ef9229c15d29e0b4e51a0e")
//        ##ALT=<ID=db22dfc14799b1aa666eb7d571cf04ec,Description="haplotype data for line: Ref",Number=6,Source="data/test/smallseq/Ref.fa",Contig=2,Start=16501,End=17500,Checksum=Md5,RefRange=db22dfc14799b1aa666eb7d571cf04ec>
        assertTrue(altHeaders.containsKey("db22dfc14799b1aa666eb7d571cf04ec"))
        val currentHeaderdb = altHeaders["db22dfc14799b1aa666eb7d571cf04ec"]
        assertEquals(currentHeaderdb?.id, "db22dfc14799b1aa666eb7d571cf04ec")
        assertEquals(currentHeaderdb?.description, "haplotype data for line: Ref")
        assertEquals(currentHeaderdb?.source, "data/test/smallseq/Ref.fa")
        assertEquals(currentHeaderdb?.regions?.get(0), Pair(Position("2",16501), Position("2",17500)))
        assertEquals(currentHeaderdb?.checksum, "Md5")
        assertEquals(currentHeaderdb?.refRange, "db22dfc14799b1aa666eb7d571cf04ec")

    //        ##ALT=<ID=5812acb1aff74866003656316c4539a6,Description="haplotype data for line: Ref",Number=6,Source="data/test/smallseq/Ref.fa",Contig=2,Start=1,End=1000,Checksum=Md5,RefRange=5812acb1aff74866003656316c4539a6>
        assertTrue(altHeaders.containsKey("5812acb1aff74866003656316c4539a6"))
        val currentHeader581 = altHeaders["5812acb1aff74866003656316c4539a6"]
        assertEquals(currentHeader581?.id, "5812acb1aff74866003656316c4539a6")
        assertEquals(currentHeader581?.description, "haplotype data for line: Ref")
        assertEquals(currentHeader581?.source, "data/test/smallseq/Ref.fa")
        assertEquals(currentHeader581?.regions?.get(0), Pair(Position("2",1), Position("2",1000)))
        assertEquals(currentHeader581?.checksum, "Md5")
        assertEquals(currentHeader581?.refRange, "5812acb1aff74866003656316c4539a6")


        assertFailsWith<IllegalStateException>(
            message = "No exception found when Testing Source",
            block = {
                parseALTHeader(VCFHeader(setOf(VCFAltHeaderLine("<ID=id," +
                        "Description=\"haplotype data for line: testSample\">," +
                        "Source_bad=\"archive.agc\"," +
                        "Regions=\"1:200-300\"," +
                        "Checksum=\"Md5\"," +
                        "RefRange=\"id\">", VCFHeaderVersion.VCF4_2))))
            }
        )
        assertFailsWith<IllegalStateException>(
            message = "No exception found when Testing checksum",
            block = {
                parseALTHeader(VCFHeader(setOf(VCFAltHeaderLine("<ID=id," +
                        "Description=\"haplotype data for line: testSample\">," +
                        "Source=\"archive.agc\"," +
                        "Regions=\"1:200-300\"," +
                        "Checksum_bad=\"Md5\"," +
                        "RefRange=\"id\">", VCFHeaderVersion.VCF4_2))))
            }
        )
        assertFailsWith<IllegalStateException>(
            message = "No exception found when Testing refRange",
            block = {
                parseALTHeader(VCFHeader(setOf(VCFAltHeaderLine("<ID=id," +
                        "Description=\"haplotype data for line: testSample\">," +
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
        val altHeaders= parseALTHeader(vcfReader.header)
        
        val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"

        val hapSequence = createFastaFromHvcf.createHaplotypeSequences(dbPath, "Ref", vcfReader.iterator().asSequence().toList(), altHeaders,"")

        assertEquals(40, hapSequence.size)
        val truthHashes = altHeaders.values.map { it.id }.toSet()

        //This verifies that we do indeed extract out the correct sequences
        hapSequence.forEach{
            assertTrue(truthHashes.contains(it.id))
        }
    }

    @Test
    fun testCreateHaplotypeSequencesFromImputedVCf() {
        // This file has a samplename of "LineImpute".
        // All of the haplotypes for chrom1 are from LineA, all of the haplotypes for chrom2 are from LineB
        val refHVCFFile = File("data/test/smallseq/LineImpute.h.vcf")
        val vcfReader = VCFFileReader(refHVCFFile, false)
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val altHeaders= parseALTHeader(vcfReader.header)

        val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"

        val hapSequence = createFastaFromHvcf.createHaplotypeSequences(dbPath, "LineImpute", vcfReader.iterator().asSequence().toList(), altHeaders,"")

        assertEquals(37, hapSequence.size)
        val truthHashes = altHeaders.values.map { it.id }.toSet()

        //This verifies that we do indeed extract out the correct sequences
        hapSequence.forEach{
            assertTrue(truthHashes.contains(it.id))
        }
    }

    @Test
    fun testHaplotypeSequenceMissingHaps() {
        // This file has a samplename of "LineImputeMissingHaps".
        // All of the haplotypes for chrom1 are from LineA, all of the haplotypes for chrom2 are from LineB
        // In addition, there are 2 missing chrom1 haplotypes in the file
        // This is a copy of the LineImpute.h.vcf file with the 2nd and 3rd haplotypes marked as missing
        val refHVCFFile = File("data/test/smallseq/LineImputeMissingHaps.h.vcf")
        val vcfReader = VCFFileReader(refHVCFFile, false)
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val altHeaders= parseALTHeader(vcfReader.header)

        val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"

        val hapSequence = createFastaFromHvcf.createHaplotypeSequences(dbPath, "LineImputeMissingHaps", vcfReader.iterator().asSequence().toList(), altHeaders,"")

        // There are 37 haplotypes in the file, but 2 are missing, so we should only get 35 haplotypes
        assertEquals(35, hapSequence.size)
        val truthHashes = altHeaders.values.map { it.id }.toSet()

        //This verifies that we do indeed extract out the correct sequences
        hapSequence.forEach{
            assertTrue(truthHashes.contains(it.id))
        }
    }

    @Test
    fun testMissingAltHeaders() {
        // This file has a samplename of "LineImputeMissingAH".
        // All of the haplotypes for chrom1 are from LineA, all of the haplotypes for chrom2 are from LineB
        // This is a copy of the LineImpute.h.vcf file with the 4th haplotype having a bad hapid, and
        // not represented in the ALT Headers
        val refHVCFFile = File("data/test/smallseq/LineImputeMissingAH.h.vcf")
        val vcfReader = VCFFileReader(refHVCFFile, false)
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val altHeaders= parseALTHeader(vcfReader.header)

        val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"

        assertThrows<IllegalStateException> {
            //Check that an error is thrown when the variants contain a hapId that is not in the ALT Headers
            createFastaFromHvcf.createHaplotypeSequences(dbPath, "LineImputeMissingAH", vcfReader.iterator().asSequence().toList(), altHeaders,"")

        }
    }

    @Test
    fun testBuildFastaFromHVCF_fileInput() {
        //buildFastaFromHVCF(dbPath: String, outputFile: String, fastaType:String, hvcfFile : String)
        val refHVCFFileName = "data/test/smallseq/Ref.h.vcf"
        val vcfReader = VCFFileReader(File(refHVCFFileName), false)
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val altHeaders= parseALTHeader(vcfReader.header)

        val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"

        createFastaFromHvcf.buildFastaFromHVCF(dbPath, "${TestExtension.testOutputFastaDir}", "composite","",refHVCFFileName,"")

        //Compare the composite against the truth input
        val truthFasta = NucSeqIO("data/test/smallseq/Ref.fa").readAll()
        val outputFastaComposite = NucSeqIO("${TestExtension.testOutputFastaDir}/Ref.fa").readAll()
        for(chr in truthFasta.keys) {
            assertEquals(getChecksumForString(truthFasta[chr]!!.seq()), getChecksumForString(outputFastaComposite[chr]!!.seq()))
        }

        //build a haplotype fasta as well
        createFastaFromHvcf.buildFastaFromHVCF(dbPath, "${TestExtension.testOutputFastaDir}", "haplotype","",refHVCFFileName,"")

        //get truth hashes:
        val truthHashes = altHeaders.values.map { it.id }.toSet()
        //load in the haplotypes and build hashes
        val outputFastaHaplotypes = NucSeqIO("${TestExtension.testOutputFastaDir}/Ref.fa").readAll()
        outputFastaHaplotypes.values.map { getChecksumForString(it.seq()) }.toSet().forEach{
            assertTrue(truthHashes.contains(it))
        }

    }

    @Test
    fun testBuildFastaFromHVCF_dirInput() {
        // This test builds fastas for multiple hvcf files in a directory
        // It verifies that the output files are created and that the sequences are correct,
        // and there are individual fasta files for each hvcf file.

        // The first hvcf file is the Ref hvcf file.  Ref.h.vcf stored to test/data/smallseq was created by
        // splitting the genome based on the anchors.bed file, then creating haplotypes for those regions.
        // When creating a fasta based on that h.vcf file, the output fasta should be the same as the original
        //
        // The second hvcf file is the LineB hvcf file.  LineB.h.vcf stored to test/data/smallseq was created by
        // aligning LineB to Ref.fa and creating haplotypes. When creating a fasta based on
        // that h.vcf file, we can't guarantee the results will be identical to the original LineB.fa file.
        // THis is because some pieces of sequence may be left out of the alignment, or insertions may be moved
        // to a different place (included as part of other haplotype), etc.
        //
        // Because of this, we create our own LineB.h.vcf file by running the LineB.fa file through
        // CreateRefVcf.createRefHvcf().  This will create a LineB.h.vcf file the same
        // way the Ref.h.vcf file was created.  The haplotype sequence in the LineB.h.vcf file should be the same
        // sequence as the original LineB.fa file.

        val refHVCFFileName = "data/test/smallseq/Ref.h.vcf"
        val lineBFasta = "data/test/smallseq/LineB.fa"

        // Write the Ref hvcf file to the output directory: TestExtension.testOutputHVCFDir/hvcf_files
        // We add the subfolder "hcf_files" as this is a subfolder to which createRefHvcf will write the hvcf files.
        val hvcfDir = "${TestExtension.testOutputHVCFDir}/hvcf_files"
        val outputFastaDir = TestExtension.testOutputFastaDir
        File(hvcfDir).mkdirs()
        File(outputFastaDir).mkdirs()
        File(TestExtension.testOutputHVCFDir+"/hvcf_files/Ref.h.vcf").writeBytes(File(refHVCFFileName).readBytes())

        // Create LineB.h.vcf file
        // We don't use the LIneB.h.vcf stored in test/smallseq because we aren't sure how it was
        // created, but probably by aligning LineB to Ref.fa.  The haplotype sequences may be slightly different
        // from what is stored in LineB.fa originally.  Instead create a new LineB.h.vcf file that has the
        // haploytpes created the same way the ref haplotypes were created.

        val createRefVcf = CreateRefVcf()
        val ranges = "data/test/smallseq/anchors.bed"

        createRefVcf.createRefHvcf(ranges,lineBFasta,"LineB","",TestExtension.testOutputHVCFDir)

        // createRefHvcf will compress the output file.
        val lineBHvcfFileName = TestExtension.testOutputHVCFDir+"/hvcf_files/LineB.h.vcf.gz"

        val refVcfReader = VCFFileReader(File(refHVCFFileName), false)
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val refAltHeaders= parseALTHeader(refVcfReader.header)

        val lineBVcfReader = VCFFileReader(File(lineBHvcfFileName), false)
        val lineBAltHeaders= parseALTHeader(lineBVcfReader.header)

        val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"

        createFastaFromHvcf.buildFastaFromHVCF(dbPath, "${TestExtension.testOutputFastaDir}", "composite",hvcfDir,"","")

        //Compare the composite against the truth input
        val refTruthFasta = NucSeqIO("data/test/smallseq/Ref.fa").readAll()
        val refOutputFastaComposite = NucSeqIO("${TestExtension.testOutputFastaDir}/Ref.fa").readAll()
        for(chr in refTruthFasta.keys) {
            assertEquals(getChecksumForString(refTruthFasta[chr]!!.seq()), getChecksumForString(refOutputFastaComposite[chr]!!.seq()))
        }

        // Compare the LineB composite against the truth input
        val lineBTruthFasta = NucSeqIO("data/test/smallseq/LineB.fa").readAll()
        val lineBOutputFastaComposite = NucSeqIO("${TestExtension.testOutputFastaDir}/LineB.fa").readAll()
        for(chr in lineBTruthFasta.keys) {
            println("Compare Composite lineB chr: $chr")
            assertEquals(getChecksumForString(lineBTruthFasta[chr]!!.seq()), getChecksumForString(lineBOutputFastaComposite[chr]!!.seq()))
        }

        //build a haplotype fasta as well
        createFastaFromHvcf.buildFastaFromHVCF(dbPath, "${TestExtension.testOutputFastaDir}", "haplotype",hvcfDir,"","")

        //get truth hashes:
        val truthHashes = refAltHeaders.values.map { it.id }.toSet()
        //load in the haplotypes and build hashes
        val outputFastaHaplotypes = NucSeqIO("${TestExtension.testOutputFastaDir}/Ref.fa").readAll()
        outputFastaHaplotypes.values.map { getChecksumForString(it.seq()) }.toSet().forEach{
            assertTrue(truthHashes.contains(it))
        }

        // Get LineB truth hashes:
        val lineBTruthHashes = lineBAltHeaders.values.map { it.id }.toSet()
        // Load in the LineB haplotypes and build hashes
        val lineBOutputFastaHaplotypes = NucSeqIO("${TestExtension.testOutputFastaDir}/LineB.fa").readAll()
        lineBOutputFastaHaplotypes.values.map { getChecksumForString(it.seq()) }.toSet().forEach{
            assertTrue(lineBTruthHashes.contains(it))
        }

    }


    @Test
    fun testExtractInversions() {
        val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"
        val sample = SampleGamete("Ref")
        val seqs = NucSeqIO("data/test/smallseq/Ref.fa").readAll()

        //create a list of a single haplotype variant with inverted regions
        val hvcfRecord = createHVCFRecord("Ref", Position("1",1),Position("1",100),  Pair("A","2b4590f722ef9229c15d29e0b4e51a0e"))
        val variantList = listOf<VariantContext>(hvcfRecord)
        //create a map of altHeaders
        val altHeader = AltHeaderMetaData("2b4590f722ef9229c15d29e0b4e51a0e","\"haplotype data for line: Ref\"","\"data/test/smallseq/Ref.fa\"",sample,
            listOf(Pair(Position("1",1), Position("1",50)), Pair(Position("1",60), Position("1", 100))), "Md5", "2b4590f722ef9229c15d29e0b4e51a0e")

        val altHeaders = mapOf<String, AltHeaderMetaData>("2b4590f722ef9229c15d29e0b4e51a0e" to altHeader)

        val createFastaFromHvcf = CreateFastaFromHvcf()

        val sequences = createFastaFromHvcf.createHaplotypeSequences(dbPath, sample.name, variantList, altHeaders,"")

        val firstSeq = sequences.first()
        assertEquals(2, firstSeq.asmRegions.size)
        assertEquals(Position("1",1), firstSeq.asmRegions[0].first)
        assertEquals(Position("1",50), firstSeq.asmRegions[0].second)
        assertEquals(Position("1",60), firstSeq.asmRegions[1].first)
        assertEquals(Position("1",100), firstSeq.asmRegions[1].second)

        assertEquals("2b4590f722ef9229c15d29e0b4e51a0e", firstSeq.id)
        assertEquals("${seqs["1"]!![0..49]}${seqs["1"]!![59..99]}", firstSeq.sequence)


        //Test a variant with a single inversion
        val hvcfRecord2 = createHVCFRecord("Ref", Position("1",1),Position("1",100),  Pair("A","db22dfc14799b1aa666eb7d571cf04ec"))
        val variantList2 = listOf<VariantContext>(hvcfRecord2)
        //create a map of altHeaders
        val altHeader2 =AltHeaderMetaData("db22dfc14799b1aa666eb7d571cf04ec","\"haplotype data for line: Ref\"","\"data/test/smallseq/Ref.fa\"",sample,
            listOf(Pair(Position("1",50), Position("1",1)), Pair(Position("1",100), Position("1", 60))), "Md5", "db22dfc14799b1aa666eb7d571cf04ec")

        val altHeaders2 = mapOf<String, AltHeaderMetaData>("db22dfc14799b1aa666eb7d571cf04ec" to altHeader2)

        val sequences2 = createFastaFromHvcf.createHaplotypeSequences(dbPath, sample.name, variantList2, altHeaders2,"")

        val firstSeq2 = sequences2.first()
        assertEquals(2, firstSeq2.asmRegions.size)
        assertEquals(Position("1",50), firstSeq2.asmRegions[0].first)
        assertEquals(Position("1",1), firstSeq2.asmRegions[0].second)
        assertEquals(Position("1",100), firstSeq2.asmRegions[1].first)
        assertEquals(Position("1",60), firstSeq2.asmRegions[1].second)

        assertEquals("db22dfc14799b1aa666eb7d571cf04ec", firstSeq2.id)
        assertEquals("${seqs["1"]!![0..49].reverse_complement()}${seqs["1"]!![59..99].reverse_complement()}", firstSeq2.sequence)


        //Test a variant with mixed normal and inversion
        val hvcfRecord3 = createHVCFRecord("Ref", Position("1",1),Position("1",100),  Pair("A","5812acb1aff74866003656316c4539a6"))
        val variantList3 = listOf<VariantContext>(hvcfRecord3)
        //create a map of altHeaders

        val altHeader3 =AltHeaderMetaData("5812acb1aff74866003656316c4539a6","\"haplotype data for line: Ref\"","\"data/test/smallseq/Ref.fa\"",sample,
            listOf(Pair(Position("1",1), Position("1",50)), Pair(Position("1",100), Position("1", 60))), "Md5", "5812acb1aff74866003656316c4539a6")

        val altHeaders3 = mapOf<String, AltHeaderMetaData>("5812acb1aff74866003656316c4539a6" to altHeader3)

        val sequences3 = createFastaFromHvcf.createHaplotypeSequences(dbPath, sample.name, variantList3, altHeaders3,"")

        val firstSeq3 = sequences3.first()
        assertEquals(2, firstSeq3.asmRegions.size)
        assertEquals(Position("1",1), firstSeq3.asmRegions[0].first)
        assertEquals(Position("1",50), firstSeq3.asmRegions[0].second)
        assertEquals(Position("1",100), firstSeq3.asmRegions[1].first)
        assertEquals(Position("1",60), firstSeq3.asmRegions[1].second)

        assertEquals("5812acb1aff74866003656316c4539a6", firstSeq3.id)
        assertEquals("${seqs["1"]!![0..49]}${seqs["1"]!![59..99].reverse_complement()}", firstSeq3.sequence)

    }

    @Test
    fun testHapSequenceRetrieval() {
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val seqs = NucSeqIO("data/test/smallseq/Ref.fa").readAll()

        //make a simple hapSeq Objects
        val firstHapSeq = HaplotypeSequence("id1","", "id1", "1", 1, 100,
            listOf(Pair(Position("1", 1), Position("1", 50)),
                    Pair(Position("1", 60), Position("1", 100))))

        val region = listOf("Ref@1:0-49","Ref@1:59-99")

        val extractedSeqs = mapOf<Pair<String,String>,NucSeq>(Pair("Ref","1:0-49") to seqs["1"]!![0..49],
            Pair("Ref","1:59-99") to seqs["1"]!![59 .. 99])

        val hapSeq = createFastaFromHvcf.buildHapSeq(extractedSeqs, region, firstHapSeq)
        assertEquals(91, hapSeq.length)
        assertEquals(extractedSeqs.values.joinToString("") { it.seq() }, hapSeq)

        //test inversion
        val secondHapSeq = HaplotypeSequence("id2","", "id2", "1", 1, 100,
            listOf(Pair(Position("1", 50), Position("1", 1)),
                Pair(Position("1", 100), Position("1", 60))))

        val region2 = listOf("Ref@1:49-0","Ref@1:59-99")

        val extractedSeqs2 = mapOf<Pair<String,String>,NucSeq>(Pair("Ref","1:49-0") to seqs["1"]!![0..49],
            Pair("Ref","1:59-99") to seqs["1"]!![59 .. 99])

        val hapSeq2 = createFastaFromHvcf.buildHapSeq(extractedSeqs2, region2, secondHapSeq)
        assertEquals(91, hapSeq2.length)
        assertEquals(extractedSeqs2.values.joinToString("") { it.reverse_complement().seq() }, hapSeq2)


        // test mixed normal and inversion
        val thirdHapSeq = HaplotypeSequence("id3","", "id3", "1", 1, 100,
            listOf(Pair(Position("1", 1), Position("1", 50)),
                Pair(Position("1", 100), Position("1", 60))))

        val region3 = listOf("Ref@1:0-49","Ref@1:59-99")

        val extractedSeqs3 = mapOf<Pair<String,String>,NucSeq>(Pair("Ref","1:0-49") to seqs["1"]!![0..49],
            Pair("Ref","1:59-99") to seqs["1"]!![59 .. 99])

        val hapSeq3 = createFastaFromHvcf.buildHapSeq(extractedSeqs3, region3, thirdHapSeq)
        assertEquals(91, hapSeq3.length)
        val expectedSeq = extractedSeqs3[Pair("Ref","1:0-49")]!!.seq() + extractedSeqs3[Pair("Ref","1:59-99")]!!.reverse_complement().seq()
        assertEquals(expectedSeq, hapSeq3)

    }
}