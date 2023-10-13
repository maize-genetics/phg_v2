package net.maizegenetics.phgv2.cli

import biokotlin.genome.getMAFblocks
import com.github.ajalt.clikt.testing.test
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

        assertEquals(result.statusCode, 0, "status code not 0: ${result.statusCode}")

        assertTrue(
            File(TestExtension.smallseqLineAMafFile).exists(),
            "File ${TestExtension.smallseqLineAMafFile} does not exist"
        )

        assertTrue(
            File(TestExtension.smallseqLineBMafFile).exists(),
            "File ${TestExtension.smallseqLineBMafFile} does not exist"
        )

        // val mafBlocks = getMAFblocks("${TestExtension.tempDir}/LineA.maf")

        // mafBlocks[0].forEach { println(it) }

    }

}