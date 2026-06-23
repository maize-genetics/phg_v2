package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader
import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader.Ps4gGameteSet
import org.junit.jupiter.api.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals

class DiploidPS4GEmissionProbabilityTest {

    @Test
    fun testIndexCountsForOneIndex() {
        val gameteSetList = mutableListOf<Ps4gGameteSet>()
        gameteSetList.add(Ps4gGameteSet(intArrayOf(0,1), 4))
        gameteSetList.add(Ps4gGameteSet(intArrayOf(2,3), 2))
        gameteSetList.add(Ps4gGameteSet(intArrayOf(0,1,2), 3))
        gameteSetList.add(Ps4gGameteSet(intArrayOf(0,1,3), 1))
        gameteSetList.add(Ps4gGameteSet(intArrayOf(1,3), 5))

        val readMap = mapOf(Pair(1, gameteSetList))
        val diploidProb = DiploidPS4GEmissionProbability(readMap, (0..3).toSet(), 0.95)

        var sampleCount = diploidProb.indexCountsForOneIndex(gameteSetList, 0)
        assertContentEquals(intArrayOf(8,15), sampleCount)
        sampleCount = diploidProb.indexCountsForOneIndex(gameteSetList, 1)
        assertContentEquals(intArrayOf(13,15), sampleCount)
        sampleCount = diploidProb.indexCountsForOneIndex(gameteSetList, 2)
        assertContentEquals(intArrayOf(5,15), sampleCount)
        sampleCount = diploidProb.indexCountsForOneIndex(gameteSetList, 3)
        assertContentEquals(intArrayOf(8,15), sampleCount)

    }

    @Test
    fun testIndexCountsForTwoIndexes() {
        val gameteSetList = mutableListOf<Ps4gGameteSet>()
        gameteSetList.add(Ps4gGameteSet(intArrayOf(0,1), 4))
        gameteSetList.add(Ps4gGameteSet(intArrayOf(2,3), 2))
        gameteSetList.add(Ps4gGameteSet(intArrayOf(0,1,2), 3))
        gameteSetList.add(Ps4gGameteSet(intArrayOf(0,1,3), 1))
        gameteSetList.add(Ps4gGameteSet(intArrayOf(1,3), 5))

        val readMap = mapOf(Pair(1, gameteSetList))
        val diploidProb = DiploidPS4GEmissionProbability(readMap, (0..3).toSet(), 0.95)

        var sampleCount = diploidProb.indexCountsForTwoIndices(gameteSetList, 0,0)
        assertContentEquals(intArrayOf(0,0,8,7), sampleCount)
        sampleCount = diploidProb.indexCountsForTwoIndices(gameteSetList, 0,1)
        assertContentEquals(intArrayOf(0,5,8,2), sampleCount)
        sampleCount = diploidProb.indexCountsForTwoIndices(gameteSetList, 0,2)
        assertContentEquals(intArrayOf(5,2,3,5), sampleCount)
        sampleCount = diploidProb.indexCountsForTwoIndices(gameteSetList, 0,3)
        assertContentEquals(intArrayOf(7,7,1,0), sampleCount)
        sampleCount = diploidProb.indexCountsForTwoIndices(gameteSetList, 1,1)
        assertContentEquals(intArrayOf(0,0,13,2), sampleCount)
    }
}