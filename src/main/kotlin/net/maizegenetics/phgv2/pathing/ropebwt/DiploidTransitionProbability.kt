package net.maizegenetics.phgv2.pathing.ropebwt

import kotlin.math.ln

class DiploidTransitionProbability(val pNoSwitch: Double, val inbreedingCoef: Double, numberOfGenomes: Int) {
    val pSwitch = (1.0 - pNoSwitch) / (numberOfGenomes - 1)
    val pss = pSwitch * pSwitch
    val psn = pSwitch * pNoSwitch
    val pnn = pNoSwitch * pNoSwitch

    fun calculate(from: Pair<Int,Int>, to: Pair<Int,Int>): Double {
        val transitionP =  when (inbreedingCoef) {
            0.0 -> probabilityForF0(from, to)
            1.0 -> probabilityForF1(from, to)
            else -> (1.0 - inbreedingCoef) * probabilityForF0(from, to) + inbreedingCoef * probabilityForF1(from, to)
        }

        return ln(transitionP)
    }

    fun probabilityForF0(from: Pair<Int,Int>, to: Pair<Int,Int>):Double {
        return if (from.first == to.first) {
            if (from.second == to.second) pnn else psn
        } else {
            if (from.second == to.second) psn else pss
        }
    }

    fun probabilityForF1(from: Pair<Int,Int>, to: Pair<Int,Int>):Double {
        return if (from.first == from.second)  {  //from is homozygous
            if (to.first == to.second) {
                if (from.first == to.first) pNoSwitch else pSwitch
            } else 0.0
        } else { //from is heterozygous
            if (to.first == to.second) {
                if (from.first == to.first || from.second == to.second) pSwitch
                else pss
            } else 0.0
        }

    }

}