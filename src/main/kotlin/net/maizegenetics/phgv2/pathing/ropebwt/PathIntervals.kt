package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.utils.Position

/**
 * A run of consecutive bins that share the same call, expressed in 0-based half-open
 * [start, end) reference coordinates per the BED specification.
 *
 * @param T the type of the call. Haploid paths use the parent name; diploid paths use an
 *   ordered (parent1, parent2) pair.
 */
data class PathInterval<T>(val contig: String, val start: Int, val end: Int, val call: T)

/**
 * Converts an ordered, per-bin path into reference-coordinate intervals.
 *
 * The cut between bin i and bin i+1 is placed at the midpoint of the two bins,
 * `((pos[i] + pos[i+1]) / 2) * binSize`. Each cut is the exclusive end of the interval before it
 * and the inclusive start of the interval after it, so the returned intervals are gapless and
 * non-overlapping. The first interval starts at 0 and the last ends at `pos[last] * binSize`.
 *
 * When [mergeAdjacent] is true, consecutive bins whose calls are equal are emitted as a single
 * interval and the cuts between them are dropped. Equality is whatever `==` means for [T], so for
 * diploid paths (lineA:0, lineB:0) and (lineB:0, lineA:0) are different calls and are not merged,
 * preserving phase.
 *
 * @param path (position, call) pairs for a single contig, ordered by bin position.
 * @param binSize the bin size used to create the ps4g file.
 * @param mergeAdjacent if false, every bin is emitted as its own interval.
 * @return the intervals in path order, or an empty list if [path] is empty.
 */
fun <T> pathToIntervals(
    path: List<Pair<Position, T>>,
    binSize: Int,
    mergeAdjacent: Boolean = true
): List<PathInterval<T>> {

    if (path.isEmpty()) return emptyList()

    //the cut between bin index and bin index + 1: exclusive end of the former, inclusive start of the latter
    fun cut(index: Int) = (path[index].first.position + path[index + 1].first.position) / 2 * binSize

    val intervals = mutableListOf<PathInterval<T>>()
    var runStart = 0

    for (index in path.indices) {
        val isLast = index == path.lastIndex
        val endsRun = isLast || !mergeAdjacent || path[index].second != path[index + 1].second
        if (!endsRun) continue

        intervals.add(
            PathInterval(
                path[runStart].first.contig,
                if (runStart == 0) 0 else cut(runStart - 1),
                if (isLast) path[index].first.position * binSize else cut(index),
                path[runStart].second
            )
        )
        runStart = index + 1
    }

    return intervals
}
