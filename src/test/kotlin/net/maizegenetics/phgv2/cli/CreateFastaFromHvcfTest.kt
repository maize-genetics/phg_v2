package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
import htsjdk.variant.vcf.VCFFileReader
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
class CreateFastaFromHvcfTest {
    companion object {

    }

    @Test
    fun testCliktParams() {
        val createFastaFromHvcf = CreateFastaFromHvcf()
        val resultGood = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} -o ${TestExtension.testOutputRefFasta} --fasta-type composite")

        val resultMissingDB = createFastaFromHvcf.test("-o ${TestExtension.testOutputRefFasta} --fasta-type composite")
        assertEquals(resultMissingDB.statusCode, 1)
        assertEquals("Usage: create-fasta-from-hvcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --db-path: --db-path must not be blank\n",resultMissingDB.output)

        val resultMissingOutput = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} --fasta-type haplotype")
        assertEquals(resultMissingOutput.statusCode, 1)
        assertEquals("Usage: create-fasta-from-hvcf [<options>]\n" +
                "\n" +
                "Error: invalid value for --output: --output/-o must not be blank\n",resultMissingOutput.output)

        val resultMissingFastaType = createFastaFromHvcf.test("--db-path ${TestExtension.testTileDBURI} -o ${TestExtension.testOutputRefFasta}")
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
        assertEquals(currentHeader2b?.contig, "1")
        assertEquals(currentHeader2b?.start, 11001)
        assertEquals(currentHeader2b?.end, 12000)
        assertEquals(currentHeader2b?.checksum, "Md5")
        assertEquals(currentHeader2b?.refRange, "2b4590f722ef9229c15d29e0b4e51a0e")
//        ##ALT=<ID=db22dfc14799b1aa666eb7d571cf04ec,Description="haplotype data for line: Ref",Number=6,Source="data/test/smallseq/Ref.fa",Contig=2,Start=16501,End=17500,Checksum=Md5,RefRange=db22dfc14799b1aa666eb7d571cf04ec>
        assertTrue(altHeaders.containsKey("db22dfc14799b1aa666eb7d571cf04ec"))
        val currentHeaderdb = altHeaders["db22dfc14799b1aa666eb7d571cf04ec"]
        assertEquals(currentHeaderdb?.id, "db22dfc14799b1aa666eb7d571cf04ec")
        assertEquals(currentHeaderdb?.description, "\"haplotype data for line: Ref\"")
        assertEquals(currentHeaderdb?.number, "6")
        assertEquals(currentHeaderdb?.source, "\"data/test/smallseq/Ref.fa\"")
        assertEquals(currentHeaderdb?.contig, "2")
        assertEquals(currentHeaderdb?.start, 16501)
        assertEquals(currentHeaderdb?.end, 17500)
        assertEquals(currentHeaderdb?.checksum, "Md5")
        assertEquals(currentHeaderdb?.refRange, "db22dfc14799b1aa666eb7d571cf04ec")

    //        ##ALT=<ID=5812acb1aff74866003656316c4539a6,Description="haplotype data for line: Ref",Number=6,Source="data/test/smallseq/Ref.fa",Contig=2,Start=1,End=1000,Checksum=Md5,RefRange=5812acb1aff74866003656316c4539a6>
        assertTrue(altHeaders.containsKey("5812acb1aff74866003656316c4539a6"))
        val currentHeader581 = altHeaders["5812acb1aff74866003656316c4539a6"]
        assertEquals(currentHeader581?.id, "5812acb1aff74866003656316c4539a6")
        assertEquals(currentHeader581?.description, "\"haplotype data for line: Ref\"")
        assertEquals(currentHeader581?.number, "6")
        assertEquals(currentHeader581?.source, "\"data/test/smallseq/Ref.fa\"")
        assertEquals(currentHeader581?.contig, "2")
        assertEquals(currentHeader581?.start, 1)
        assertEquals(currentHeader581?.end, 1000)
        assertEquals(currentHeader581?.checksum, "Md5")
        assertEquals(currentHeader581?.refRange, "5812acb1aff74866003656316c4539a6")

    }


}