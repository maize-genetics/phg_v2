package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader

class MostLikelyPs4gParents(val ps4gReader: Ps4gFileReader, val chromosomeSet: Set<String>) {

    fun bestParents(numberOfParents: Int) : Set<Int> {
        require(numberOfParents > 1) {"Number of best parents must be greater than 1."}
        val bestParentSet = mutableSetOf<Int>()
        val gameteSets = mutableListOf<Ps4gFileReader.Ps4gGameteSet>()
        for (contig in chromosomeSet) {
            val contigReadMap = ps4gReader.readMapForContig(contig)
            check(contigReadMap != null) { "contigReadMap is null for contig = $contig" +
                    " caused by an incorrect contig name passed to MostLikelyParents." }

            gameteSets.addAll(contigReadMap.entries.flatMap { it.value })
        }
        //1. count hits to each gamete
        //2. choose parent with the highest count
        //3. add parent to best parent list and remove from count map

        val gameteCounts = gameteCountsFromGameteSets(gameteSets).toMutableMap()
        val chosenParent =  gameteCounts.maxByOrNull { it.value }?.key
        check(chosenParent != null) { "chosen Parent is null in MostLikelyParents." }
        bestParentSet.add(chosenParent)
        gameteCounts.remove(chosenParent)

//        4. count # of reads hitting each parent, filtering out reads hitting chosen parent
//        5. choose parent with the highest count
//        6. add this parent to best parent list and remove it from gameteCounts
        val nextParent = bestParentFromFilteredGameteSets(gameteSets, chosenParent)
        check(nextParent != null) { "next Parent is null in MostLikelyParents." }
        bestParentSet.add(nextParent)
        gameteCounts.remove(nextParent)

        //repeat process until enough parents have been chosen
        while (bestParentSet.size < numberOfParents) {
            val highestCountParent = gameteCounts.maxByOrNull { it.value }?.key
            check(highestCountParent != null) {"Could only choose ${bestParentSet.size} best parents."}
            bestParentSet.add(highestCountParent)
            gameteCounts.remove(highestCountParent)

            if (bestParentSet.size < numberOfParents) {
                val complimentParent = bestParentFromFilteredGameteSets(gameteSets, chosenParent)
                check(complimentParent != null) {"Could only choose ${bestParentSet.size} best parents."}
                bestParentSet.add(complimentParent)
                gameteCounts.remove(complimentParent)
            }

        }

        return bestParentSet
    }

    fun bestParentFromFilteredGameteSets(gameteSets: List<Ps4gFileReader.Ps4gGameteSet>, excludedParent: Int): Int? {
        val filteredGameteSets = gameteSets.filter { !it.gameteIndices.contains(excludedParent) }
        val gameteCounts = gameteCountsFromGameteSets(filteredGameteSets)
        return gameteCounts.maxByOrNull { it.value }?.key
    }

    fun gameteCountsFromGameteSets(gameteSets: List<Ps4gFileReader.Ps4gGameteSet>): Map<Int, Int> {
        return gameteSets.flatMap { gameteSet ->
            gameteSet.gameteIndices.map { index ->
                Pair(
                    index,
                    gameteSet.count
                )
            }
        }.groupingBy { it.first }.eachCount().toMutableMap()
    }

}