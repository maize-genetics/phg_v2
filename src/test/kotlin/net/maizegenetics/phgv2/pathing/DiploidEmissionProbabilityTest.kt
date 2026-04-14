package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.extension.ExtendWith
import kotlin.math.ln
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class DiploidEmissionProbabilityTest {

    @Test
    fun testSmallFactorials() {
        val maxdiff = 1e-10
        val expected = listOf(0.0, ln(1.0), ln(2.0), ln(6.0), ln(24.0))
        for (ndx in 0..4) {
            println("expected = ${expected[ndx]}, actual = ${DiploidEmissionProbability.smallFactorials[ndx]}")
            assertEquals(expected[ndx], DiploidEmissionProbability.smallFactorials[ndx], maxdiff)
        }
    }

    @Test
    fun testmultinomialProbabilities() {

        val maxdiff = 1e-5
        //This tests lnMultinomialProbability(counts: IntArray, probabilities: DoubleArray): Double
        //The following values were calculated using R:
        //counts    probabilities   lnProbability
        //(1,1,2,1)     (.3,.3,.3,.1)   -3.024132
        //(2,0,1,1)     (.3,.3,.3,.1)   -3.429597
        //(2,0,3,1)     (.3,.3,.3,.1)   -4.228105
        //(3,2,15,2)     (.3,.3,.3,.1)  -11.29077
        //(5,7,25,5)     (.3,.3,.3,.1)  -14.39179
        var lnprob1 = DiploidEmissionProbability.lnMultinomialProbability(intArrayOf(1,1,2,1), doubleArrayOf(.3,.3,.3,.1))
        print("expected = -3.024132, actual = $lnprob1\n")
        assertEquals(-3.024132, lnprob1, maxdiff)

        lnprob1 = DiploidEmissionProbability.lnMultinomialProbability(intArrayOf(2,0,1,1), doubleArrayOf(.3,.3,.3,.1))
        print("expected = -3.429597, actual = $lnprob1\n")
        assertEquals(-3.429597, lnprob1, maxdiff)

        lnprob1 = DiploidEmissionProbability.lnMultinomialProbability(intArrayOf(2,0,3,1), doubleArrayOf(.3,.3,.3,.1))
        print("expected = -4.228105, actual = $lnprob1\n")
        assertEquals(-4.228105, lnprob1, maxdiff)

        lnprob1 = DiploidEmissionProbability.lnMultinomialProbability(intArrayOf(3,2,15,2), doubleArrayOf(.3,.3,.3,.1))
        print("expected = -11.29077, actual = $lnprob1\n")
        assertEquals(-11.29077, lnprob1, maxdiff)

        lnprob1 = DiploidEmissionProbability.lnMultinomialProbability(intArrayOf(5,7,25,5), doubleArrayOf(.3,.3,.3,.1))
        print("expected = -14.39179, actual = $lnprob1\n")
        assertEquals(-14.39179, lnprob1, maxdiff)

    }
}