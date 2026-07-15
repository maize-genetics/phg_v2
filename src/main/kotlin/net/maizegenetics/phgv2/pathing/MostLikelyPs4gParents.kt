package net.maizegenetics.phgv2.pathing

import net.maizegenetics.phgv2.pathing.ropebwt.Ps4gFileReader

/**
 * Given the read mappings in a PS4G file, finds the most likely parents (gametes) in a stepwise fashion.
 * The gamete hit by the most reads is selected as a likely parent. Then the subset of reads that do not hit that
 * parent is created, and the gamete with the most reads in that subset is added as the second likely parent. That
 * second step selects a parent complementary to the first, which is useful for diploid path finding where the two
 * parents are expected to cover different reads. Additional parents are added in order of decreasing read count
 * until numberOfParents have been chosen.
 *
 * Parents are identified by their gamete index, as used in the PS4G file. [Ps4gFileReader.gameteIndexMap] maps
 * those indices back to gamete names.
 *
 * @param ps4gReader    a reader for the PS4G file holding the read mappings to be counted
 * @param chromosomeSet the names of the contigs whose read mappings will be used. Every name must be present in
 *                      the PS4G file; small contigs or scaffolds are typically excluded.
 */
class MostLikelyPs4gParents(val ps4gReader: Ps4gFileReader, val chromosomeSet: Set<String>) {

    /**
     * Returns the gamete indices of the [numberOfParents] most likely parents, using the read mappings for the
     * contigs in [chromosomeSet].
     *
     * @param numberOfParents the number of parents to choose. Must be greater than 1.
     * @throws IllegalArgumentException if [numberOfParents] is not greater than 1.
     * @throws IllegalStateException if a contig in [chromosomeSet] is not in the PS4G file, or if the file does not
     * contain mappings to enough distinct gametes to choose [numberOfParents] parents.
     */
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

    /**
     * Returns the gamete index hit by the most gamete sets in [gameteSets], after discarding every gamete set that
     * also hits [excludedParent]. The result is the gamete that best complements [excludedParent], in the sense that
     * it accounts for the most reads that [excludedParent] does not explain.
     *
     * @return the best complementary gamete index, or null if every gamete set hits [excludedParent].
     */
    fun bestParentFromFilteredGameteSets(gameteSets: List<Ps4gFileReader.Ps4gGameteSet>, excludedParent: Int): Int? {
        val filteredGameteSets = gameteSets.filter { !it.gameteIndices.contains(excludedParent) }
        val gameteCounts = gameteCountsFromGameteSets(filteredGameteSets)
        return gameteCounts.maxByOrNull { it.value }?.key
    }

    /**
     * Tallies the read mapping hits to each gamete in [gameteSets].
     *
     * @return a map of gamete index to the number of gamete sets in [gameteSets] that contain that index.
     */
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