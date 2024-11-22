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
                "Error: missing option --output-dir\n",resultMissingOutput.output)

        val resultMissingFastaType = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} -o ${TestExtension.tempDir} --hvcf-file /test_file.h.vcf")
        assertEquals(resultMissingFastaType.statusCode, 1)
        assertEquals("Usage: create-fasta-from-hvcf [<options>]\n" +
                "\n" +
                "Error: missing option --fasta-type\n",resultMissingFastaType.output)

        val resultBadFastaType = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} -o ${TestExtension.tempDir} --fasta-type bad --hvcf-file /test_file.h.vcf")
        assertEquals(resultBadFastaType.statusCode, 1)
        assertEquals("Usage: create-fasta-from-hvcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --fasta-type: invalid choice: bad. (choose from composite, haplotype, pangenomeHaplotype)\n",resultBadFastaType.output)

        // Only one of --hvcf-file or --hvcf-dir can be used.  THis test has both sent
        val resultHvcfFileAndDir = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} -o ${TestExtension.tempDir} --fasta-type haplotype --hvcf-file /test_file.h.vcf --hvcf-dir /test_dir")
        assertEquals(resultHvcfFileAndDir.statusCode, 1)
        assertEquals("Usage: create-fasta-from-hvcf [<options>]\n" +
                "\n" +
                "Error: option --hvcf-file cannot be used with --hvcf-dir\n",resultHvcfFileAndDir.output)

        val resultNoHvcfInput = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} -o ${TestExtension.tempDir} --fasta-type haplotype ")
        assertEquals(resultNoHvcfInput.statusCode, 1)
        assertEquals("Usage: create-fasta-from-hvcf [<options>]\n" +
                "\n" +
                "Error: must provide one of --hvcf-file, --hvcf-dir\n",resultNoHvcfInput.output)
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
        // THere is a duplicate here - it will be removed when writing the file
        val haplotypeSequences =  listOf(
            HaplotypeSequence(getChecksumForString(seqs[0]), seqs[0], getChecksumForString(seqs[0]), "1", 1, 300, listOf(Pair(Position("1", 1), Position("1", 300)))),
            HaplotypeSequence(getChecksumForString(seqs[1]), seqs[1], getChecksumForString(seqs[1]), "1", 301, 600, listOf(Pair(Position("1", 301), Position("1", 600)))),
            HaplotypeSequence(getChecksumForString(seqs[3]), seqs[3], getChecksumForString(seqs[3]), "2", 1, 300, listOf(Pair(Position("2", 1), Position("1",300)))),
            HaplotypeSequence(getChecksumForString(seqs[2]), seqs[2], getChecksumForString(seqs[2]), "1", 601, 900, listOf(Pair(Position("1", 601), Position("1",900)))),
            HaplotypeSequence(getChecksumForString(seqs[3]), seqs[3], getChecksumForString(seqs[3]), "2", 1, 300, listOf(Pair(Position("2", 1), Position("1",300)))),
            HaplotypeSequence(getChecksumForString(seqs[4]), seqs[4], getChecksumForString(seqs[4]), "2", 301, 600, listOf(Pair(Position("2", 301), Position("1", 600))))
        )

        // Write the haplotype fasta
        BufferedWriter(FileWriter(outputFile)).use { writer ->
            createFastaFromHvcf.writeHaplotypeSequence(writer, haplotypeSequences)
        }

        //Check that the file was created
        val file = File(outputFile)
        assertTrue(file.exists())

        //load in the file and check that the sequences are correct
        val importedSeqs = NucSeqIO(outputFile).readAll()

        // Verify there are only 5 sequences in the file (the duplicate was removed)
        assertEquals(5, importedSeqs.size)

        assertEquals(seqs[0], importedSeqs[getChecksumForString(seqs[0])]!!.seq())
        assertEquals(seqs[1], importedSeqs[getChecksumForString(seqs[1])]!!.seq())
        assertEquals(seqs[2], importedSeqs[getChecksumForString(seqs[2])]!!.seq())
        assertEquals(seqs[3], importedSeqs[getChecksumForString(seqs[3])]!!.seq())
        assertEquals(seqs[4], importedSeqs[getChecksumForString(seqs[4])]!!.seq())

        // writeHaplotypeSequences processed the list as a set, and orders are not maintained in a set.
        // Verify that each value exists at some place in the file.
        compareFastaDescriptions(importedSeqs[getChecksumForString(seqs[0])]!!.description?:"", haplotypeSequences.firstOrNull { it.id == getChecksumForString(seqs[0]) }!!)
        compareFastaDescriptions(importedSeqs[getChecksumForString(seqs[1])]!!.description?:"", haplotypeSequences.firstOrNull { it.id == getChecksumForString(seqs[1]) }!!)
        compareFastaDescriptions(importedSeqs[getChecksumForString(seqs[2])]!!.description?:"", haplotypeSequences.firstOrNull { it.id == getChecksumForString(seqs[2]) }!!)
        compareFastaDescriptions(importedSeqs[getChecksumForString(seqs[3])]!!.description?:"", haplotypeSequences.firstOrNull { it.id == getChecksumForString(seqs[3]) }!!)
        compareFastaDescriptions(importedSeqs[getChecksumForString(seqs[4])]!!.description?:"", haplotypeSequences.firstOrNull { it.id == getChecksumForString(seqs[4]) }!!)

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

        val hvcfInput = HvcfInput.HvcfFile(refHVCFFileName)
        createFastaFromHvcf.buildFastaFromHVCF(dbPath, "${TestExtension.testOutputFastaDir}", CreateFastaFromHvcf.FastaType.composite, hvcfInput,"")

        //Compare the composite against the truth input
        val truthFasta = NucSeqIO("data/test/smallseq/Ref.fa").readAll()
        val outputFastaComposite = NucSeqIO("${TestExtension.testOutputFastaDir}/Ref_composite.fa").readAll()
        for(chr in truthFasta.keys) {
            assertEquals(getChecksumForString(truthFasta[chr]!!.seq()), getChecksumForString(outputFastaComposite[chr]!!.seq()))
        }

        //build a haplotype fasta as well
        createFastaFromHvcf.buildFastaFromHVCF(dbPath, "${TestExtension.testOutputFastaDir}", CreateFastaFromHvcf.FastaType.haplotype, hvcfInput,"")

        //get truth hashes:
        val truthHashes = altHeaders.values.map { it.id }.toSet()
        //load in the haplotypes and build hashes
        val outputFastaHaplotypes = NucSeqIO("${TestExtension.testOutputFastaDir}/Ref_haplotype.fa").readAll()
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
        // sequence as the original LineB.fa file.   This file has been created and stored to data/test/createFastaFromHvcf/LineB.h.vcf

        val refHVCFFileName = "data/test/smallseq/Ref.h.vcf"
        val lineBHvcfFileName = "data/test/createFastaFromHvcf/LineB.h.vcf"

        // Write the Ref and LineB hvcf files to the output directory: TestExtension.testOutputHVCFDir
        val hvcfDir = TestExtension.testOutputHVCFDir
        val hvcfInput = HvcfInput.HvcfDir(TestExtension.testOutputHVCFDir)
        File(hvcfDir).mkdirs()
        // File(TestExtension.smallSeqLineAGvcfFile).copyTo(lineAGvcf, true)
        File(refHVCFFileName).copyTo(File(TestExtension.testOutputHVCFDir+"/Ref.h.vcf"), true)
        File(lineBHvcfFileName).copyTo(File(TestExtension.testOutputHVCFDir+"/LineB.h.vcf"),true)

        val outputFastaDir = TestExtension.testOutputFastaDir
        File(outputFastaDir).mkdirs()

        // Get the Ref and LineB ALT headers
        val refVcfReader = VCFFileReader(File(refHVCFFileName), false)
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val refAltHeaders= parseALTHeader(refVcfReader.header)

        val lineBVcfReader = VCFFileReader(File(lineBHvcfFileName), false)
        val lineBAltHeaders= parseALTHeader(lineBVcfReader.header)

        // Build fastas from the hvcf files
        val dbPath = "${TestExtension.testOutputFastaDir}/dbPath"
        createFastaFromHvcf.buildFastaFromHVCF(dbPath, "${TestExtension.testOutputFastaDir}", CreateFastaFromHvcf.FastaType.composite, hvcfInput,"")

        //Compare the ref composite against the truth input
        val refTruthFasta = NucSeqIO("data/test/smallseq/Ref.fa").readAll()
        val refOutputFastaComposite = NucSeqIO("${TestExtension.testOutputFastaDir}/Ref_composite.fa").readAll()
        for(chr in refTruthFasta.keys) {
            assertEquals(getChecksumForString(refTruthFasta[chr]!!.seq()), getChecksumForString(refOutputFastaComposite[chr]!!.seq()))
        }

        // Compare the LineB composite against the truth input
        val lineBTruthFasta = NucSeqIO("data/test/smallseq/LineB.fa").readAll()
        val lineBOutputFastaComposite = NucSeqIO("${TestExtension.testOutputFastaDir}/LineB_composite.fa").readAll()
        for(chr in lineBTruthFasta.keys) {
            println("Compare Composite lineB chr: $chr")
            assertEquals(getChecksumForString(lineBTruthFasta[chr]!!.seq()), getChecksumForString(lineBOutputFastaComposite[chr]!!.seq()))
        }

        //build the haplotype fastas as well
        createFastaFromHvcf.buildFastaFromHVCF(dbPath, "${TestExtension.testOutputFastaDir}", CreateFastaFromHvcf.FastaType.haplotype, hvcfInput,"")

        //get ref truth hashes:
        val truthHashes = refAltHeaders.values.map { it.id }.toSet()
        //load in the haplotypes and compare hashes
        val outputFastaHaplotypes = NucSeqIO("${TestExtension.testOutputFastaDir}/Ref_haplotype.fa").readAll()
        outputFastaHaplotypes.values.map { getChecksumForString(it.seq()) }.toSet().forEach{
            assertTrue(truthHashes.contains(it))
        }

        // Get LineB truth hashes:
        val lineBTruthHashes = lineBAltHeaders.values.map { it.id }.toSet()
        // Load in the LineB haplotypes and compare hashes
        val lineBOutputFastaHaplotypes = NucSeqIO("${TestExtension.testOutputFastaDir}/LineB_haplotype.fa").readAll()
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

        val hapSeq = buildHapSeq(extractedSeqs, region, firstHapSeq)
        assertEquals(91, hapSeq.length)
        assertEquals(extractedSeqs.values.joinToString("") { it.seq() }, hapSeq)

        //test inversion
        val secondHapSeq = HaplotypeSequence("id2","", "id2", "1", 1, 100,
            listOf(Pair(Position("1", 50), Position("1", 1)),
                Pair(Position("1", 100), Position("1", 60))))

        val region2 = listOf("Ref@1:49-0","Ref@1:59-99")

        val extractedSeqs2 = mapOf<Pair<String,String>,NucSeq>(Pair("Ref","1:49-0") to seqs["1"]!![0..49],
            Pair("Ref","1:59-99") to seqs["1"]!![59 .. 99])

        val hapSeq2 = buildHapSeq(extractedSeqs2, region2, secondHapSeq)
        assertEquals(91, hapSeq2.length)
        assertEquals(extractedSeqs2.values.joinToString("") { it.reverse_complement().seq() }, hapSeq2)


        // test mixed normal and inversion
        val thirdHapSeq = HaplotypeSequence("id3","", "id3", "1", 1, 100,
            listOf(Pair(Position("1", 1), Position("1", 50)),
                Pair(Position("1", 100), Position("1", 60))))

        val region3 = listOf("Ref@1:0-49","Ref@1:59-99")

        val extractedSeqs3 = mapOf<Pair<String,String>,NucSeq>(Pair("Ref","1:0-49") to seqs["1"]!![0..49],
            Pair("Ref","1:59-99") to seqs["1"]!![59 .. 99])

        val hapSeq3 = buildHapSeq(extractedSeqs3, region3, thirdHapSeq)
        assertEquals(91, hapSeq3.length)
        val expectedSeq = extractedSeqs3[Pair("Ref","1:0-49")]!!.seq() + extractedSeqs3[Pair("Ref","1:59-99")]!!.reverse_complement().seq()
        assertEquals(expectedSeq, hapSeq3)

    }
}