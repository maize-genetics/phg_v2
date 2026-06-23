package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.cli.TestExtension
import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader
import org.junit.jupiter.api.AfterAll
import org.junit.jupiter.api.BeforeAll
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import java.io.File
import kotlin.random.Random
import kotlin.test.assertContentEquals

@ExtendWith(TestExtension::class)
class MostLikelyPs4gParentsTest {

    companion object {
        @JvmStatic
        @BeforeAll
        fun setup() {
            File(TestExtension.testOutputDir).deleteRecursively()
            File(TestExtension.testOutputDir).mkdirs()
        }

        @JvmStatic
        @AfterAll
        fun tearDown() {
            //comment out the following line to inspect the test results after the tests have been run
//            File(TestExtension.testOutputDir).deleteRecursively()
        }
    }

    val ps4gFilename = "${TestExtension.testOutputDir}likely-parent.ps4g"

    @Test
    fun testLikelyParents() {
        writePs4gTestFile()
        val mlp = MostLikelyPs4gParents(Ps4gFileReader(ps4gFilename), setOf("chr1"))
        val bestParentList = mlp.bestParents(2)
        println(bestParentList)
        assertContentEquals(listOf(0,1), bestParentList)
    }

    fun writePs4gTestFile() {
        File(ps4gFilename).bufferedWriter(Charsets.UTF_8).use { writer ->
            writer.write("#PS4G\n")
            writer.write("#version=2.0\n")
            writer.write("#Command: PHGv2 Command:\n")
            writer.write("#TotalUniqueCounts: 50\n" +
                    "#gamete\tgameteIndex\tcount\n")

            //generate some reads
            val numberOfGametes = 5
            val gameteNames = listOf("gamete1", "gamete2", "gamete3", "gamete4", "gamete5")
            val contig = "chr1"
            val numberOfPositions = 25
            val numberOfReadsAtAPosition = 2
            val probNonParent = 0.5

            val readMappingString = mutableListOf<String>()
            val gameteCounts = intArrayOf(0,0,0,0,0)

            for (pos in 0 until numberOfPositions) {

                for (readNumber in 0 until 2) {
                    val tempList = mutableListOf<Int>()
                    if (Random.nextDouble() < 0.5) {
                        tempList.add(0)
                        gameteCounts[0]++
                        if (Random.nextDouble() < probNonParent) {
                            tempList.add(1)
                            gameteCounts[1]++
                        }
                    } else {
                        if (Random.nextDouble() < probNonParent) {
                            tempList.add(0)
                            gameteCounts[0]++
                        }
                        tempList.add(1)
                        gameteCounts[1]++
                    }
                    if (Random.nextDouble() < probNonParent) {
                        tempList.add(2)
                        gameteCounts[2]++
                    }
                    if (Random.nextDouble() < probNonParent) {
                        tempList.add(3)
                        gameteCounts[3]++
                    }
                    if (Random.nextDouble() < probNonParent) {
                        tempList.add(4)
                        gameteCounts[4]++
                    }
                    readMappingString.add("${tempList.joinToString(",")}\tchr1\t$pos\t1")
                }


            }

            //write the gametes
            for (ndx in 0 until numberOfGametes) {
                writer.write("#${gameteNames[ndx]}\t$ndx\t${gameteCounts[ndx]}\n")
            }

            //write the set counts
            writer.write("gameteSet\trefContig\trefPosBinned\tcount\n")
            readMappingString.forEach { writer.write("$it\n") }

        }
    }
}