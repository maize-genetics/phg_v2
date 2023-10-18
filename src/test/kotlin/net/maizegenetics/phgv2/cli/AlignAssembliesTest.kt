package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.testing.test
// import net.maizegenetics.phgv2.utils.testMergingMAF
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.test.assertEquals
import kotlin.test.assertTrue

@ExtendWith(TestExtension::class)
class AlignAssembliesTest {

    @Test
    fun testRunningAlignAssemblies() {

        // phg align-assemblies --gff /workdir/tmc46/AlignAssemblies/smallSeq_data/anchors.gff
        // --ref /workdir/tmc46/AlignAssemblies/smallSeq_data/Ref.fa
        // -a /workdir/tmc46/AlignAssemblies/assembliesList.txt
        // -o /workdir/tmc46/AlignAssemblies/temp

        val alignAssemblies = AlignAssemblies()

        val result = alignAssemblies.test(
            "--gff ${TestExtension.smallseqAnchorsGffFile} --ref ${TestExtension.smallseqRefFile} " +
                    "-a ${TestExtension.smallseqAssembliesListFile} -o ${TestExtension.tempDir}"
        )

        println("testRunningAlignAssemblies: result output: ${result.output}")

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        val lineAMAF = TestExtension.tempDir + "LineA.maf"
        assertTrue(File(lineAMAF).exists(), "File $lineAMAF does not exist")

        val lineBMAF = TestExtension.tempDir + "LineB.maf"
        assertTrue(File(lineBMAF).exists(), "File $lineBMAF does not exist")

        // val mafOutput = TestExtension.tempDir + "LineA_unsplit.maf"
        // println("mafOutput: $mafOutput")
        // testMergingMAF(TestExtension.smallseqLineAMafFile, mafOutput)
        // testMergingMAF(lineAMAF, TestExtension.tempDir + "LineA_unsplit2.maf")

        // val mafBlocks = getMAFblocks("${TestExtension.tempDir}/LineA.maf")

        // mafBlocks[0].forEach { println(it) }

    }

}